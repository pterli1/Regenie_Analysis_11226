"""Schemas of fedcat's intermediate files, and readers and converters for REGENIE-style results."""

from __future__ import annotations

import asyncio
import tarfile
from collections import OrderedDict
from pathlib import Path

import polars as pl
from loguru import logger

from fedcat.core.logging import NCPU, Permits

# ── File schemas ──────────────────────────────────────────────────────────────

# Conditioning step input (tophits_input), written without a header
TOPHITS_DTYPES = pl.Schema(
    {
        "ID": pl.Utf8,
        "LOG10P_IN": pl.Float64,
        "A1FREQ_IN": pl.Float64,
        "BETA_IN": pl.Float64,
        "SE_IN": pl.Float64,
        "TEST": pl.Utf8,
        "CHROM": pl.Utf8,
        "START": pl.Int32,
        "END": pl.Int32,
    }
)
TOPHITS_HEADER = list(TOPHITS_DTYPES.keys())

# IVW meta-analysis output
META_OUTPUT_DTYPES = pl.Schema(
    {
        "CHROM": pl.Utf8,
        "POS": pl.Int64,
        "ID": pl.Utf8,
        "A0": pl.Utf8,
        "A1": pl.Utf8,
        "A1FREQ": pl.Float64,
        "N": pl.Int32,
        "TEST": pl.Utf8,
        "BETA": pl.Float64,
        "SE": pl.Float64,
        "Z": pl.Float64,
        "LOG10P": pl.Float64,
        "Q": pl.Float64,
        "Q_P": pl.Float64,
        "I2": pl.Float64,
        "N_STUDIES": pl.Int32,
    }
)

# ACAT output
ACAT_OUTPUT_DTYPES = pl.Schema(
    {
        "CHR": pl.Utf8,
        "START": pl.Int32,
        "END": pl.Int32,
        "ID": pl.Utf8,
        "UNIT": pl.Utf8,
        "MASK": pl.Utf8,
        "THRESHOLD": pl.Utf8,
        "N_AGG_COMBINED": pl.Int32,
        "N_SINGLEVAR_COMBINED": pl.Int32,
        "LOG10P_ACAT": pl.Float64,
    }
)

# Conditioning step summary_output: each unit before (_IN) and after (_OUT) the step
SUMMARY_OUTPUT_DTYPES = pl.Schema(
    {
        "ID": pl.Utf8,
        "LOG10P_IN": pl.Float64,
        "A1FREQ_IN": pl.Float64,
        "BETA_IN": pl.Float64,
        "SE_IN": pl.Float64,
        "TEST": pl.Utf8,
        "CHROM": pl.Utf8,
        "START": pl.Int32,
        "END": pl.Int32,
        "A1FREQ_OUT": pl.Float64,
        "BETA_OUT": pl.Float64,
        "SE_OUT": pl.Float64,
        "LOG10P_OUT": pl.Float64,
        "STEP": pl.Int32,
        "JUMPER": pl.Boolean,
    }
)

# Per-variant collinearity with the selected set
COLLINEARITY_DTYPES = pl.Schema(
    {
        "ID": pl.Utf8,
        "cand_rsq": pl.Float64,
        "max_rsq": pl.Float64,
    }
)
COLLINEARITY_HEADER = list(COLLINEARITY_DTYPES.keys())


def parse_cohort_file(path: Path) -> dict[str, list[Path]]:
    """Parse a cohort file listing each cohort's result files.

    Each `# <cohort>` line starts a cohort, and the following non-empty lines are its
    files:

    ```text
    # cohort_name_1
    /path/to/file1.regenie.gz
    /path/to/file2.regenie.gz
    # cohort_name_2
    /path/to/file3.regenie.gz
    ```

    Args:
        path: Cohort file.

    Returns:
        Files per cohort label, in file order.

    Raises:
        ValueError: If a file comes before any cohort, a cohort name is empty, or there
            are no cohorts.
    """
    cohorts: OrderedDict[str, list[Path]] = OrderedDict()
    current_cohort: str | None = None

    with open(path) as fh:
        for lineno, raw_line in enumerate(fh, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                current_cohort = line.lstrip("#").strip()
                if not current_cohort:
                    raise ValueError(f"{path}:{lineno}: empty cohort name after '#'")
                if current_cohort not in cohorts:
                    cohorts[current_cohort] = []
                continue
            if current_cohort is None:
                raise ValueError(f"{path}:{lineno}: file path before any '# cohort_name' header")
            cohorts[current_cohort].append(Path(line))

    if not cohorts:
        raise ValueError(f"{path}: no cohorts found (file is empty or has no '# cohort_name' headers)")
    return dict(cohorts)


# ── Tarball helpers ───────────────────────────────────────────────────────────


def extract_tarball(tar_path: Path, destination_folder: Path) -> None:
    """Extract a tarball into a new folder.

    Args:
        tar_path: Tarball, compressed or not.
        destination_folder: Folder to create and extract into.

    Raises:
        FileExistsError: If `destination_folder` exists.
    """
    if destination_folder.exists():
        raise FileExistsError(f"Cannot extract tarball to {destination_folder}: already exists")
    with tarfile.open(tar_path, "r:*") as tar:
        tar.extractall(path=destination_folder)


async def extract_tarball_async(tar_path: Path, destination_folder: Path) -> None:
    """Run [`extract_tarball`][fedcat.core.io.extract_tarball] in a thread."""
    async with Permits(1):
        await asyncio.to_thread(extract_tarball, tar_path, destination_folder)


# ── Manc-COJO output reformatter ──────────────────────────────────────────────

# Columns of a single-cohort .cma.cojo, in file order. manc-cojo runs once per stratum,
# so its multi-cohort form is not read.
MANCOJO_CMA_DTYPES = pl.Schema(
    {
        "Chr": pl.Utf8,
        "SNP": pl.Utf8,
        "bp": pl.Int64,
        "A1": pl.Utf8,
        "A2": pl.Utf8,
        "freq": pl.Float64,
        "b": pl.Float64,
        "se": pl.Float64,
        "p": pl.Float64,
        # Effective sample size, which is fractional
        "n": pl.Float64,
        "bC": pl.Float64,
        "bC_se": pl.Float64,
        "pC": pl.Float64,
    }
)


# Columns of a converted .cma.cojo, as in REGENIE output
MANCOJO_OUTPUT_SCHEMA = pl.Schema(
    {
        "CHROM": pl.Utf8,
        "GENPOS": pl.Int64,
        "ID": pl.Utf8,
        "ALLELE0": pl.Utf8,
        "ALLELE1": pl.Utf8,
        "A1FREQ": pl.Float64,
        "N": pl.Int64,
        "TEST": pl.Utf8,
        "BETA": pl.Float64,
        "SE": pl.Float64,
        "LOG10P": pl.Float64,
    }
)


def reformat_mancojo_output(files: list[Path], output_path: Path, test_label: str = "ADD") -> int:
    """Convert Manc-COJO `.cma.cojo` conditional results to REGENIE's column layout.

    The output can then be meta-analysed like REGENIE results, by
    [`meta_analyse_regenie`][fedcat.stats.meta_analysis.meta_analyse_regenie].
    `LOG10P` is computed from `bC / bC_se`, as manc-cojo's `pC` underflows to 0 below
    about 1e-308. Rows without a usable conditional estimate are dropped.

    Args:
        files: Single-cohort `.cma.cojo` files.
        output_path: Output file, space-separated with the `MANCOJO_OUTPUT_SCHEMA` columns.
        test_label: Value of the `TEST` column.

    Returns:
        Number of rows written; 0 means the output only has a header and should not be
        meta-analysed.
    """
    from fedcat.stats.meta_analysis import (
        z_to_log10p,
    )

    df = pl.concat(
        [
            pl.read_csv(
                f,
                separator="\t",
                has_header=True,
                schema_overrides=MANCOJO_CMA_DTYPES,
                null_values=["NA", "nan", "-nan", ""],
            )
            for f in files
        ]
    )
    # A zero SE would give an infinite Z
    df = df.filter(pl.col("bC").is_not_null() & pl.col("bC_se").is_not_null() & (pl.col("bC_se") > 0))

    if df.height == 0:
        logger.warning(f"No usable conditional results in {[f.name for f in files]}")
        pl.DataFrame(schema=MANCOJO_OUTPUT_SCHEMA).write_csv(output_path, separator=" ", null_value="NA")
        return 0

    log10p = z_to_log10p((df["bC"] / df["bC_se"]).to_numpy())

    df = df.with_columns(
        # manc-cojo, like REGENIE, writes the X chromosome as 23
        pl.when(pl.col("Chr") == "23").then(pl.lit("X")).otherwise(pl.col("Chr")).alias("CHROM"),
        pl.lit(test_label).alias("TEST"),
        pl.Series("LOG10P", log10p),
    ).rename(
        {
            "SNP": "ID",
            "bp": "GENPOS",
            "A2": "ALLELE0",
            "A1": "ALLELE1",
            "freq": "A1FREQ",
            "bC": "BETA",
            "bC_se": "SE",
        }
    )
    df = df.with_columns(pl.col("n").round().cast(pl.Int64).alias("N"))

    df.select(list(MANCOJO_OUTPUT_SCHEMA)).write_csv(output_path, separator=" ", null_value="NA")
    return df.height


async def reformat_mancojo_output_async(files: list[Path], output_path: Path, test_label: str = "ADD") -> int:
    """Run [`reformat_mancojo_output`][fedcat.core.io.reformat_mancojo_output] in a thread."""
    async with Permits(NCPU - 1):
        return await asyncio.to_thread(reformat_mancojo_output, files, output_path, test_label)


# ── Collinearity reformatter ──────────────────────────────────────────────────


def reformat_collinearity_output(files: list[Path], output_path: Path) -> None:
    """Concatenate collinearity files (`COLLINEARITY_DTYPES` columns) into one TSV.

    Args:
        files: Collinearity files, e.g. one per chromosome.
        output_path: Output TSV.
    """
    df = pl.concat(
        [pl.read_csv(f, separator="\t", has_header=True, schema_overrides=COLLINEARITY_DTYPES) for f in files]
    )
    df.write_csv(output_path, separator="\t")


async def reformat_collinearity_output_async(files: list[Path], output_path: Path) -> None:
    """Run [`reformat_collinearity_output`][fedcat.core.io.reformat_collinearity_output] in a thread."""
    async with Permits(1):
        await asyncio.to_thread(reformat_collinearity_output, files, output_path)
