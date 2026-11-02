"""Type stubs and documentation of the `fedcat_rs` Rust extension.

Meta-analysis, ACAT and annotation conversion, used by
[`fedcat.stats`][fedcat.stats]. Inputs that cannot be read raise
`pyo3_runtime.PanicException` from the ACAT and conversion functions.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

def acat_log10p_rs(
    log10p_values: NDArray[np.float64],
    weights: NDArray[np.float64] | None = None,
) -> float:
    """Combine -log10(P) values by ACAT, as REGENIE's `get_acat_robust` does.

    Args:
        log10p_values: 1-D array of -log10(P) values.
        weights: Weights of the same length, normalised to sum to 1; None for equal weights.

    Returns:
        -log10(P) of the combination, or NaN for no values.

    Raises:
        ValueError: If `weights` is not the length of `log10p_values`.
    """
    ...

def z_to_log10p_rs(
    z_array: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Convert Z-scores to -log10 of their two-sided p-values, precise at any |Z|.

    Args:
        z_array: 1-D array of Z-scores.

    Returns:
        1-D array of -log10(P), all >= 0.
    """
    ...

def meta_analyse_regenie_streaming_rs(
    cohort_files: list[tuple[str, list[str]]],
    output_path: str,
    res_type: str,
    column_config: dict | None = None,
    compression: str = "zst",
    compression_level: int | None = None,
    aaf_bins: list[str] | None = None,
) -> tuple[int, int]:
    """Meta-analyse REGENIE results by IVW fixed effects, accumulating one cohort at a time.

    Rows with SE <= 0, a missing or non-finite value, or EXTRA `TEST_FAIL` are dropped,
    and chromosome 23 becomes X. Results are combined on (CHROM, position, ID, TEST); a
    variant in one cohort more than once, in one file or several, counts once and must have
    the same BETA and SE each time. Output
    is tab-separated after a `## Cohorts:` line, with columns CHROM, POS, ID, A0, A1,
    A1FREQ (N-weighted), N, TEST, BETA, SE, Z, LOG10P, Q, Q_P (`NA` from one study), I2,
    N_STUDIES and DIRECTIONS (`+`, `-` or `?` per cohort), sorted by chromosome, position, ID and test.

    Args:
        cohort_files: `(label, files)` per cohort.
        output_path: Output file.
        res_type: Result type, for the log only.
        column_config: Input column names where they differ from REGENIE's, keyed by
            `chrom`, `genpos`, `id`, `allele0`, `allele1`, `a1freq`, `n`, `test`, `beta`,
            `se` and `extra` (`""` for none); and `separator` (detected if not given).
            A `log10p` key is accepted but unused.
        compression: `"zst"`, `"gz"`, or anything else for none.
        compression_level: Compression level; None for 3 (zst) or 6 (gz).
        aaf_bins: Allele-frequency bins, e.g. `["0.001", "0.01"]`. A cohort contributes
            to bin `t` only where its own A1FREQ is `<= t`, and the output has a row per
            variant and bin, with an `AAF_BIN` column.

    Returns:
        `(n_rows_written, n_cohorts)`.

    Raises:
        ValueError: If there are no input files, a bin is not a number in (0, 1], or a
            cohort gives one variant different results.
        RuntimeError: If Polars cannot read or process the input.
        OSError: If the output cannot be written.
    """
    ...

def meta_analyse_regenie_merge_rs(
    cohort_files: list[tuple[str, list[str]]],
    output_path: str,
    res_type: str,
    column_config: dict | None = None,
    compression: str = "zst",
    compression_level: int | None = None,
    aaf_bins: list[str] | None = None,
) -> tuple[int, int]:
    """Meta-analyse sorted REGENIE results by IVW fixed effects, merging all files at once.

    Uses memory proportional to the number of files. Input must be sorted by chromosome,
    position, ID and TEST, as REGENIE writes it; out-of-order rows are warned about. Output,
    arguments and the handling of repeated variants are as for
    [`meta_analyse_regenie_streaming_rs`][fedcat_rs.meta_analyse_regenie_streaming_rs].

    Returns:
        `(n_rows_written, n_cohorts)`.

    Raises:
        ValueError: If there are no input files, a column is missing, or a cohort gives
            one variant different results.
        OSError: If a file cannot be read or written.
    """
    ...

def compute_acat_streaming_rs(
    singles_tsv: str,
    burden_tsv: str,
    schema_dir: str,
    aggregate_freqs: list[str],
    weighting_scheme: str = "standard",
    output_path: str = "",
    column_config: dict | None = None,
    compression: str = "zst",
    compression_level: int | None = None,
    backend: str = "sweep",
) -> dict[str, int]:
    """Compute the ACAT of every unit, mask and frequency bin, and write them as a TSV.

    Each ACAT combines a mask's single-variant results in the bin with its burden
    result, if any (test `ADD`, ID `UNIT.MASK.BIN`). Output columns are those of
    `ACAT_OUTPUT_DTYPES` in [`fedcat.core.io`][fedcat.core.io].

    Args:
        singles_tsv: Single-variant meta-analysis results. With an `AAF_BIN` column each
            row is used only in its own bin; otherwise bins filter on A1FREQ.
        burden_tsv: Burden (MAC-10) meta-analysis results, or `""` for none.
        schema_dir: Annotation schema directory from `fedcat convert-annots`.
        aggregate_freqs: Frequency bin labels, e.g. `["0.001", "0.01"]`, written as
            REGENIE labels the bins.
        weighting_scheme: `"standard"` (REGENIE and GCTA weights, maf × (1 - maf)^49)
            or `"equal"`.
        output_path: Output TSV, or `""` to write nothing.
        column_config: Input column names, keyed by `chrom`, `id`, `a1freq`, `test` and
            `log10p`; and `separator` (default tab).
        compression: `"zst"`, `"gz"`, or anything else for none.
        compression_level: Compression level; None for 3 (zst) or 6 (gz).
        backend: `"sweep"`, for schemas sorted by unit start as `convert_tsv_to_parquet_rs`
            writes them, or `"arrow"`, which needs no sort order but allows at most 64
            annotations and masks.

    Returns:
        Counts `_n_singles`, `_n_burden`, `_n_annots` and `_n_results`.

    Raises:
        ValueError: If `aggregate_freqs` is empty or not numbers.
    """
    ...

def convert_tsv_to_parquet_rs(
    annotation_tsvs: list[str],
    masks_tsvs: list[str],
    output_dir: str,
    valid_variant_ids: list[str] | None = None,
    full_bounds: bool = True,
) -> dict[str, int]:
    """Convert REGENIE annotation and mask files to the annotation parquet schema.

    The schema's files are described in [`fedcat.stats.acat`][fedcat.stats.acat].
    Rows whose variant ID is not `chr:pos:ref:alt` are skipped, and
    `annotations.parquet` is sorted by chromosome, unit start, unit and position.

    Args:
        annotation_tsvs: REGENIE annotation files (`variant<tab>unit<tab>annotation`),
            plain, gzip or zstd.
        masks_tsvs: REGENIE mask files (`MASK<tab>annot1,annot2,...`).
        output_dir: Schema directory, created if needed.
        valid_variant_ids: Variant IDs to keep; None for all.
        full_bounds: With `valid_variant_ids`, take each mask's positions from all its
            annotated variants rather than only those kept.

    Returns:
        Counts `n_rows`, `n_variants`, `n_units` and `n_annots`, and `full_bounds`:
        whether all variants' positions were used.
    """
    ...
