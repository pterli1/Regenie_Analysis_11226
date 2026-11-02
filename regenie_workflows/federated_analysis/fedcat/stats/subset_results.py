"""Step-0 subsets: the significant meta-analysis results conditioning starts from.

For a phenotype, the significant burden, ACAT and unconditional single-variant results
are written, as `tophits_input`, to an aggregate and a singles tarball. The aggregate
tarball also holds the annotations of the significant units: an `annotations_parquet/`
schema directory, or REGENIE `annotations` and `masks` files.
"""

from __future__ import annotations

import re
import shutil
import tarfile
from pathlib import Path

import polars as pl
from loguru import logger

from fedcat.core.io import TOPHITS_HEADER

DEFAULT_LOG10P_THRESHOLD = 6.0
# Frequency bin suffix of an aggregate ID, stripped to get its UNIT.MASK group
_AGG_BIN_SUFFIX = re.compile(r"(\.singleton|\.0\.\d+)$")

# Files of an annotation parquet schema directory
_PARQUET_SCHEMA_FILES = (
    "unit_mask_bounds.parquet",
    "dim_variants.parquet",
    "dim_units.parquet",
    "dim_annots.parquet",
    "mask_annotations.parquet",
    "annotations.parquet",
)


def _optional_str_col(df: pl.DataFrame, col: str) -> pl.Expr:
    """Return `col` as a string column, or a null string literal if the frame lacks it."""
    return pl.col(col).cast(pl.Utf8) if col in df.columns else pl.lit(None).cast(pl.Utf8)


def _unit_from_agg_id(agg_id: str) -> str:
    """Return the unit of a `UNIT.MASK.BIN` aggregate ID; units may contain dots, masks cannot."""
    return _AGG_BIN_SUFFIX.sub("", agg_id).rsplit(".", 1)[0]


def _test_priority_expr() -> pl.Expr:
    """Return a sort key where BURDEN outranks ACAT when LOG10P ties."""
    return (
        pl.when(pl.col("TEST") == "BURDEN")
        .then(pl.lit(0))
        .when(pl.col("TEST") == "ACAT")
        .then(pl.lit(1))
        .otherwise(pl.lit(2))
        .alias("_test_priority")
    )


def _make_tarball(tar_path: Path, *members: Path) -> None:
    """Write a gzipped tarball of the existing members, each at the top level under its own name."""
    with tarfile.open(tar_path, "w:gz") as tar:
        for m in members:
            if m.exists():
                tar.add(m, arcname=m.name)


def _subset_annotations_parquet(schema_dir: Path, gene_set: set[str], output_dir: Path) -> None:
    """Write the annotation schema restricted to the units in `gene_set` to `output_dir`."""
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_list = list(gene_set)

    dim_units = pl.read_parquet(schema_dir / "dim_units.parquet")
    dim_variants = pl.read_parquet(schema_dir / "dim_variants.parquet")
    dim_annots = pl.read_parquet(schema_dir / "dim_annots.parquet")

    dim_units_sub = dim_units.filter(pl.col("UNIT").is_in(gene_list))
    kept_unit_codes = dim_units_sub["UNIT_CODE"].to_list()

    annotations = pl.read_parquet(schema_dir / "annotations.parquet")
    annotations_sub = annotations.filter(pl.col("UNIT_CODE").is_in(kept_unit_codes))
    kept_variant_codes = annotations_sub["VARIANT_CODE"].unique().to_list()
    kept_annot_codes = annotations_sub["ANNOT_CODE"].unique().to_list()

    dim_variants_sub = dim_variants.filter(pl.col("VARIANT_CODE").is_in(kept_variant_codes))
    dim_annots_sub = dim_annots.filter(pl.col("ANNOT_CODE").is_in(kept_annot_codes))

    mask_annotations = pl.read_parquet(schema_dir / "mask_annotations.parquet")
    mask_annotations_sub = mask_annotations.filter(pl.col("ANNOT_CODE").is_in(kept_annot_codes))

    unit_mask_bounds = pl.read_parquet(schema_dir / "unit_mask_bounds.parquet")
    unit_mask_bounds_sub = unit_mask_bounds.filter(pl.col("UNIT").is_in(gene_list))

    annotations_sub.write_parquet(output_dir / "annotations.parquet")
    dim_units_sub.write_parquet(output_dir / "dim_units.parquet")
    dim_variants_sub.write_parquet(output_dir / "dim_variants.parquet")
    dim_annots_sub.write_parquet(output_dir / "dim_annots.parquet")
    mask_annotations_sub.write_parquet(output_dir / "mask_annotations.parquet")
    unit_mask_bounds_sub.write_parquet(output_dir / "unit_mask_bounds.parquet")


def _merge_masks_tsv(masks_paths: list[Path], output_path: Path) -> None:
    """Merge REGENIE mask files, combining the annotations of masks with the same name.

    Missing files are skipped with a warning.

    Raises:
        ValueError: If a line is not `MASK<tab>annotations`.
    """
    masks: dict[str, set[str]] = {}
    for masks_path in masks_paths:
        if not masks_path.exists():
            logger.warning(f"Masks file not found, skipping: {masks_path}")
            continue
        with open(masks_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t", 1)
                if len(parts) != 2:
                    raise ValueError(f"Malformed mask line: {line}")
                mask_name, annot_str = parts
                if mask_name not in masks:
                    masks[mask_name] = set()
                masks[mask_name].update(annot_str.split(","))
    with open(output_path, "w") as out:
        for mask_name, annots in sorted(masks.items()):
            out.write(f"{mask_name}\t{','.join(sorted(annots))}\n")


def _subset_annotations_tsv(annotation_paths: list[Path], gene_set: set[str], output_path: Path) -> None:
    """Write the rows of REGENIE annotation files for units in `gene_set`; missing files are skipped."""
    lfs: list[pl.LazyFrame] = []
    for ann_path in annotation_paths:
        if not ann_path.exists():
            logger.warning(f"Annotation file not found, skipping: {ann_path}")
            continue
        lfs.append(
            pl.scan_csv(
                ann_path, separator="\t", has_header=False, new_columns=["VAR_ID", "GENE", "ANNOT"]
            ).filter(pl.col("GENE").is_in(list(gene_set)))
        )
    if lfs:
        pl.concat(lfs).sink_csv(output_path, separator="\t", include_header=False)
    else:
        output_path.touch()


def subset_results_for_phenotype(
    pheno: str,
    *,
    acat_path: Path | None = None,
    burden_path: Path | None = None,
    singles_path: Path | None = None,
    annotation_paths: list[Path] | None = None,
    masks_paths: list[Path] | None = None,
    schema_dir: Path | None = None,
    agg_threshold: float = DEFAULT_LOG10P_THRESHOLD,
    singles_threshold: float = DEFAULT_LOG10P_THRESHOLD,
    chr_label: str | None = None,
) -> None:
    """Write a phenotype's step-0 subset tarballs to the working directory.

    Writes `<pheno>_agg_cond_step0_subset[_chr<chr>].tar.gz` from the ACAT and burden
    results, and `<pheno>_singles_cond_step0_subset[_chr<chr>].tar.gz` from the singles,
    for the result files given. Burden results take their positions from the ACAT
    results of their unit, so need the ACAT file.

    Args:
        pheno: Phenotype name.
        acat_path: ACAT output (TSV).
        burden_path: Burden meta-analysis results (TSV).
        singles_path: Unconditional single-variant meta-analysis results (TSV).
        annotation_paths: REGENIE annotation files; the alternative to `schema_dir`.
        masks_paths: REGENIE mask files, with `annotation_paths`.
        schema_dir: Annotation schema directory; the alternative to `annotation_paths`.
        agg_threshold: Least -log10(P) of the ACAT and burden results kept.
        singles_threshold: Least -log10(P) of the single-variant results kept.
        chr_label: Chromosome added to the tarball names.

    Raises:
        FileNotFoundError: If a given file is missing.
        ValueError: If both annotation sources are given, or burden results are given
            without ACAT results.
    """
    for filepath in [acat_path, burden_path, singles_path]:
        if filepath is not None and not filepath.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")
    if annotation_paths and schema_dir:
        raise ValueError("Cannot provide both annotation_paths and schema_dir.")
    if annotation_paths is not None:
        for ann_path in annotation_paths:
            if not ann_path.exists():
                raise FileNotFoundError(f"Annotation file not found: {ann_path}")
    if schema_dir is not None:
        for schema_file in _PARQUET_SCHEMA_FILES:
            if not (schema_dir / schema_file).exists():
                raise FileNotFoundError(f"Parquet schema file not found: {schema_dir / schema_file}")

    chr_suffix = f"_chr{chr_label}" if chr_label else ""

    # ── Scan burden ───────────────────────────────────────────────────
    sig_burden_df: pl.DataFrame | None = None
    if burden_path is not None:
        logger.info(f"[{pheno}] Scanning burden: {burden_path}")
        sig_burden_df = (
            pl.scan_csv(
                burden_path,
                separator="\t",
                comment_prefix="##",
                null_values=["NA"],
                schema_overrides={"CHROM": pl.Utf8, "POS": pl.Int64},
            )
            .filter(pl.col("LOG10P").is_not_null() & (pl.col("LOG10P") >= agg_threshold))
            .collect()
        )
        logger.info(f"[{pheno}]   {len(sig_burden_df)} significant burden hits (LOG10P >= {agg_threshold:.2f})")

    needed_units: set[str] = set()
    if sig_burden_df is not None:
        for gene_id in sig_burden_df["ID"].to_list():
            needed_units.add(_unit_from_agg_id(gene_id))

    # ── Scan ACAT: significant results, and the units of significant burden results ──
    sig_acat_df: pl.DataFrame | None = None
    unit_bounds: dict[str, tuple[int, int]] = {}
    gene_set: set[str] = set()

    if needed_units and acat_path is None:
        raise ValueError("ACAT scan required to get coordinates for burden hits, but no ACAT file provided.")
    elif acat_path:
        logger.info(f"[{pheno}] Scanning ACAT: {acat_path}")

        if acat_path and needed_units:
            acat_filter = (
                pl.col("LOG10P_ACAT").is_not_null() & (pl.col("LOG10P_ACAT") >= agg_threshold)
            ) | pl.col("UNIT").is_in(list(needed_units))
        else:
            acat_filter = pl.col("LOG10P_ACAT").is_not_null() & (pl.col("LOG10P_ACAT") >= agg_threshold)

        acat_kept = (
            pl.scan_csv(
                acat_path,
                separator="\t",
                null_values=["NA"],
                schema_overrides={"CHR": pl.Utf8, "START": pl.Int64, "END": pl.Int64},
            )
            .filter(acat_filter)
            .collect()
        )

        if len(acat_kept):
            for row in (
                acat_kept.group_by("UNIT")
                .agg(pl.col("START").min().alias("s"), pl.col("END").max().alias("e"))
                .iter_rows(named=True)
            ):
                unit_bounds[row["UNIT"]] = (row["s"], row["e"])

        if acat_path:
            sig_acat_df = acat_kept.filter(
                pl.col("LOG10P_ACAT").is_not_null() & (pl.col("LOG10P_ACAT") >= agg_threshold)
            )
            gene_set.update(sig_acat_df["UNIT"].to_list())
            logger.info(f"[{pheno}]   {len(sig_acat_df)} significant ACAT hits")

    if sig_burden_df is not None:
        gene_set.update(needed_units)

    # ── Aggregate conditioning (ACAT + burden) ────────────────────────
    if acat_path or burden_path:
        agg_dir = Path(f"{pheno}_agg_cond_step0_subset")
        agg_dir.mkdir(exist_ok=True)

        tophits_frames: list[pl.DataFrame] = []

        if sig_acat_df is not None and len(sig_acat_df):
            tophits_frames.append(
                sig_acat_df.select(
                    pl.col("ID"),
                    pl.col("LOG10P_ACAT").alias("LOG10P_IN"),
                    pl.lit(None).cast(pl.Utf8).alias("A1FREQ_IN"),
                    pl.lit(None).cast(pl.Utf8).alias("BETA_IN"),
                    pl.lit(None).cast(pl.Utf8).alias("SE_IN"),
                    pl.lit("ACAT").alias("TEST"),
                    pl.col("CHR").alias("CHROM"),
                    pl.col("START"),
                    pl.col("END"),
                )
            )

        if sig_burden_df is not None and len(sig_burden_df):
            starts, ends = [], []
            for gene_id in sig_burden_df["ID"].to_list():
                coords = unit_bounds.get(_unit_from_agg_id(gene_id))
                starts.append(coords[0] if coords else None)
                ends.append(coords[1] if coords else None)

            burden_with_coords = sig_burden_df.with_columns(
                pl.Series("START", starts, dtype=pl.Int64),
                pl.Series("END", ends, dtype=pl.Int64),
            )
            tophits_frames.append(
                burden_with_coords.select(
                    pl.col("ID"),
                    pl.col("LOG10P").alias("LOG10P_IN"),
                    _optional_str_col(burden_with_coords, "A1FREQ").alias("A1FREQ_IN"),
                    _optional_str_col(burden_with_coords, "BETA").alias("BETA_IN"),
                    _optional_str_col(burden_with_coords, "SE").alias("SE_IN"),
                    pl.lit("BURDEN").alias("TEST"),
                    pl.col("CHROM"),
                    pl.col("START"),
                    pl.col("END"),
                )
            )

        tophits_path = agg_dir / "tophits_input"
        if tophits_frames:
            tophits_df = (
                pl.concat(tophits_frames, how="diagonal_relaxed")
                .with_columns(_test_priority_expr())
                .sort(["LOG10P_IN", "_test_priority"], descending=[True, False], nulls_last=True)
                # ACAT and burden share UNIT.MASK.BIN IDs; keep both tests for reference
                .unique(subset=["ID", "TEST"], keep="first")
                .sort(["LOG10P_IN", "_test_priority"], descending=[True, False], nulls_last=True)
                .select(TOPHITS_HEADER)
            )
            tophits_df.write_csv(tophits_path, separator="\t", include_header=False, null_value="NA")
            logger.info(f"[{pheno}]   Wrote {len(tophits_df)} aggregate tophits")
        else:
            tophits_path.touch()

        tarball_members: list[Path] = [tophits_path]
        masks_out = agg_dir / "masks"
        if schema_dir is not None and gene_set:
            ann_pq_dir = agg_dir / "annotations_parquet"
            _subset_annotations_parquet(schema_dir, gene_set, ann_pq_dir)
            tarball_members.append(ann_pq_dir)
            logger.info(f"[{pheno}]   Wrote subsetted parquet schema to {ann_pq_dir}")
        elif annotation_paths and gene_set:
            ann_out = agg_dir / "annotations"
            _subset_annotations_tsv(annotation_paths, gene_set, ann_out)
            tarball_members.append(ann_out)
            logger.info(f"[{pheno}]   Wrote annotations to {ann_out}")
            if masks_paths:
                _merge_masks_tsv(masks_paths, masks_out)
                tarball_members.append(masks_out)
                logger.info(f"[{pheno}]   Wrote masks to {masks_out}")
        else:
            ann_out = agg_dir / "annotations"
            ann_out.touch()
            tarball_members.append(ann_out)

        tar_path = Path(f"{pheno}_agg_cond_step0_subset{chr_suffix}.tar.gz")
        _make_tarball(tar_path, *tarball_members)
        logger.info(f"[{pheno}]   Created tarball: {tar_path}")

        shutil.rmtree(agg_dir)

    # ── Single-variant conditioning ────────────────────────────────────
    if singles_path is not None and singles_path.exists():
        logger.info(f"[{pheno}] Scanning unconditional singles: {singles_path}")
        sig_sng_df = (
            pl.scan_csv(
                singles_path,
                separator="\t",
                comment_prefix="##",
                null_values=["NA"],
                schema_overrides={"CHROM": pl.Utf8, "POS": pl.Int64},
            )
            .filter(pl.col("LOG10P").is_not_null() & (pl.col("LOG10P") >= singles_threshold))
            .collect()
        )
        logger.info(f"[{pheno}]   {len(sig_sng_df)} significant singles hits (LOG10P >= {singles_threshold:.2f})")

        sng_dir = Path(f"{pheno}_singles_cond_step0_subset")
        sng_dir.mkdir(exist_ok=True)
        tophits_path = sng_dir / "tophits_input"

        if len(sig_sng_df):
            sig_sng_df.select(
                pl.col("ID"),
                pl.col("LOG10P").alias("LOG10P_IN"),
                _optional_str_col(sig_sng_df, "A1FREQ").alias("A1FREQ_IN"),
                _optional_str_col(sig_sng_df, "BETA").alias("BETA_IN"),
                _optional_str_col(sig_sng_df, "SE").alias("SE_IN"),
                pl.lit("ADD").alias("TEST"),
                pl.col("CHROM"),
                pl.col("POS").alias("START"),
                pl.col("POS").alias("END"),
            ).write_csv(tophits_path, separator="\t", include_header=False, null_value="NA")
            logger.info(f"[{pheno}]   Wrote {len(sig_sng_df)} singles tophits")
        else:
            tophits_path.touch()

        tar_path = Path(f"{pheno}_singles_cond_step0_subset{chr_suffix}.tar.gz")
        _make_tarball(tar_path, tophits_path)
        logger.info(f"[{pheno}]   Created tarball: {tar_path}")

        shutil.rmtree(sng_dir)
