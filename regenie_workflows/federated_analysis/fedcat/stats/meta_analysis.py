"""Inverse-variance weighted (IVW) fixed-effects meta-analysis of REGENIE results.

[`meta_analyse_regenie`][fedcat.stats.meta_analysis.meta_analyse_regenie] reads
`.regenie.gz` files directly and meta-analyses them in the `fedcat_rs` extension, with
Cochran's Q and I² heterogeneity statistics. P-values are kept as -log10(P), so
extreme values do not underflow.
"""

from __future__ import annotations

import asyncio
import math
from pathlib import Path

import numpy as np
import polars as pl
from loguru import logger

from fedcat.core.logging import NCPU, Permits
from fedcat_rs import meta_analyse_regenie_merge_rs as _meta_analyse_regenie_merge_rs
from fedcat_rs import meta_analyse_regenie_streaming_rs as _meta_analyse_regenie_streaming_rs
from fedcat_rs import z_to_log10p_rs as _z_to_log10p_rs

# Meta-analysis implementations: streaming (default) or merge, which needs sorted input
_BACKENDS = {
    "streaming": _meta_analyse_regenie_streaming_rs,
    "merge": _meta_analyse_regenie_merge_rs,
}

_LN10 = math.log(10.0)


def z_to_log10p(z_array: np.ndarray, use_rust: bool = True) -> np.ndarray:
    """Convert Z-scores to -log10 of their two-sided p-values.

    Args:
        z_array: 1-D array of Z-scores.
        use_rust: Convert every value in Rust; otherwise |Z| <= 37 is converted with
            scipy, which loses precision beyond that.

    Returns:
        1-D float64 array of -log10(P), all >= 0.
    """
    z_array = np.asarray(z_array, dtype=np.float64)
    if use_rust:
        return np.asarray(_z_to_log10p_rs(z_array))
    az = np.abs(z_array)
    result = np.empty_like(z_array)

    # scipy is precise for moderate Z
    moderate = az <= 37.0
    if moderate.any():
        from scipy.stats import norm

        log_sf = norm.logsf(az[moderate])
        result[moderate] = -(log_sf + math.log(2.0)) / _LN10

    extreme = ~moderate
    if extreme.any():
        result[extreme] = np.asarray(_z_to_log10p_rs(z_array[extreme]))

    np.maximum(result, 0.0, out=result)
    return result


def p_to_log10p(p: float) -> float:
    """Convert a p-value, e.g. a config threshold, to -log10(P)."""
    return -math.log10(p)


def log10p_to_p(log10p: np.ndarray | float) -> np.ndarray:
    """Convert -log10(P) to P, clipped to [0, 1]; values above 308 underflow to 0."""
    arr = np.asarray(log10p, dtype=np.float64)
    p = np.power(10.0, -arr)
    return np.clip(p, 0.0, 1.0)


# Rows read to judge a file's variant ID form
_ID_SNIFF_ROWS = 200


def _variant_id_forms(path: Path, id_col: str, separator: str | None) -> set[str]:
    """Return the forms (`"chr"`, `"bare"`) of the `chrom:pos:...` IDs in a file's first rows.

    IDs without a `:` (aggregates, rsIDs) are ignored.
    """
    for sep in [separator] if separator else ["\t", " "]:
        try:
            ids = pl.read_csv(
                path,
                separator=sep,
                n_rows=_ID_SNIFF_ROWS,
                columns=[id_col],
                comment_prefix="##",
                infer_schema=False,
            )[id_col]
        except (pl.exceptions.PolarsError, OSError):
            continue
        variant_ids = ids.drop_nulls().filter(ids.drop_nulls().str.contains(":"))
        n_chr = variant_ids.str.contains(r"(?i)^chr[^:]*:").sum()
        return ({"chr"} if n_chr else set()) | ({"bare"} if n_chr < variant_ids.len() else set())
    return set()


def _warn_on_mixed_chr_ids(
    cohort_list: list[tuple[str, list[str]]], column_config: dict[str, str] | None
) -> None:
    """Warn when cohorts differ in a leading `chr` on variant IDs, which are matched exactly.

    Each cohort is judged by the first rows of its first file.
    """
    cfg = column_config or {}
    by_form: dict[str, list[str]] = {"chr": [], "bare": []}
    for label, files in cohort_list:
        for form in _variant_id_forms(Path(files[0]), cfg.get("id", "ID"), cfg.get("separator")):
            by_form[form].append(label)
    if by_form["chr"] and by_form["bare"]:
        logger.warning(
            "Variant IDs mix 'chr'-prefixed and bare forms, so the same variant is not combined across "
            f"them: 'chr' in {', '.join(by_form['chr'])}; bare in {', '.join(by_form['bare'])}"
        )


def meta_analyse_regenie(
    per_cohort_files: dict[str, list[Path]],
    output_file: Path,
    res_type: str,
    column_config: dict[str, str] | None = None,
    compression: str | None = None,
    compression_level: int | None = None,
    backend: str = "streaming",
    aaf_bins: list[str] | None = None,
) -> None:
    """Meta-analyse REGENIE results across cohorts by IVW fixed effects.

    Results are combined on exact variant ID, without checking alleles.

    Args:
        per_cohort_files: `.regenie.gz` files per cohort label, e.g. `"eur_ukb"`.
        output_file: Tab-separated output file: the `META_OUTPUT_DTYPES` columns, then
            `DIRECTIONS` (and `AAF_BIN` with bins).
        res_type: Result type, e.g. `"single_variants"`, `"burden"` or `"mac10"`; only logged.
        column_config: Input column names where they differ from REGENIE's, keyed by
            `chrom`, `genpos`, `id`, `allele0`, `allele1`, `a1freq`, `n`, `test`, `beta`,
            `se`; and `separator`.
        compression: Output compression, `"zst"` or `"gz"`, matching the file
            suffix; None to infer it from the suffix.
        compression_level: Compression level; None for the default.
        backend: `"streaming"` accumulates each cohort in turn; `"merge"` merges
            cohorts in sorted order, using least memory, and needs sorted input.
        aaf_bins: Allele-frequency bins, e.g. `["0.001", "0.01"]`, for single-variant
            input to ACAT. A cohort contributes to bin `t` only where its A1FREQ is
            `<= t`, and the output has a row per variant and bin, with an `AAF_BIN` column.

    Raises:
        ValueError: If there are no input files, `backend` or `compression` is invalid,
            or a cohort gives one variant different BETA or SE values.
    """
    if not per_cohort_files:
        raise ValueError("per_cohort_files must not be empty")

    backend_fn = _BACKENDS.get(backend)
    if backend_fn is None:
        raise ValueError(f"Unknown backend {backend!r}; choose from {sorted(_BACKENDS)}")

    cohort_list = [(label, [str(f) for f in files]) for label, files in per_cohort_files.items() if files]
    if not cohort_list:
        raise ValueError("No cohort frames produced — check input files")
    _warn_on_mixed_chr_ids(cohort_list, column_config)

    if compression is not None:
        if compression not in {"zst", "gz"}:
            raise ValueError(f"Invalid compression '{compression}'; if passed must be one of {'zst', 'gz'}")
        if compression != output_file.suffix[1:]:
            raise ValueError(
                f"Output file suffix '{output_file.suffix}' does not match specified compression '{compression}'. "
            )
    else:
        match output_file.suffix:
            case ".zst":
                compression = "zst"
            case ".gz":
                compression = "gz"
            case _:
                compression = "none"
        logger.info(f"Inferred compression='{compression}' from output file suffix '{output_file.suffix}'")

    n_variants, n_cohorts = backend_fn(
        cohort_list,
        str(output_file),
        res_type,
        column_config=column_config,
        compression=compression,
        compression_level=compression_level,
        aaf_bins=aaf_bins,
    )
    logger.info(
        f"Meta-analysis complete ({backend}): {n_variants} variants from {n_cohorts} cohorts → {output_file.name}"
    )


async def meta_analyse_regenie_async(
    per_cohort_files: dict[str, list[Path]],
    output_file: Path,
    res_type: str,
    column_config: dict[str, str] | None = None,
    compression: str | None = None,
    compression_level: int | None = None,
    backend: str = "streaming",
    aaf_bins: list[str] | None = None,
) -> None:
    """Run [`meta_analyse_regenie`][fedcat.stats.meta_analysis.meta_analyse_regenie] in a thread."""
    async with Permits(NCPU - 1):
        await asyncio.to_thread(
            meta_analyse_regenie,
            per_cohort_files,
            output_file,
            res_type,
            column_config=column_config,
            compression=compression,
            compression_level=compression_level,
            backend=backend,
            aaf_bins=aaf_bins,
        )
