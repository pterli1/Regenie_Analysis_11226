"""ACAT tests of aggregates, combining burden and single-variant results, and annotation conversion.

ACAT reads annotations as a parquet schema directory, written by `fedcat convert-annots`:

- `dim_variants.parquet`, `dim_units.parquet`, `dim_annots.parquet`: integer codes of
  variant IDs, units (genes or regions) and annotations.
- `annotations.parquet`: each variant's unit and annotation, with its chromosome and
  position, sorted by unit start.
- `mask_annotations.parquet`: the annotations in each mask.
- `unit_mask_bounds.parquet`: each unit and mask's first and last position.
"""

import os
import tempfile
from collections import defaultdict
from pathlib import Path

import polars as pl
from loguru import logger

from fedcat_rs import compute_acat_streaming_rs as _compute_acat_streaming_rs
from fedcat_rs import convert_tsv_to_parquet_rs as _convert_tsv_to_parquet_rs


def compute_acat(
    single_variants_file: Path,
    burden_file: Path | None = None,
    schema_dir: Path | None = None,
    masks_tsvs: list[Path] | None = None,
    annotations_tsvs: list[Path] | None = None,
    output_file: Path = Path("acat_output"),
    aggregate_freqs: list[str] | None = None,
    id_col: str = "ID",
    allele_freq_col: str = "A1FREQ",
    log10p_col: str = "LOG10P",
    n_threads: int | None = None,
    weighting_scheme: str = "standard",
    compression: str | None = None,
    compression_level: int | None = None,
    backend: str = "sweep",
    full_bounds: bool = True,
) -> None:
    """Compute ACAT p-values of every unit, mask and frequency bin.

    Each ACAT combines the mask's burden result (if any) with its rare single-variant
    results. Annotations are read as a parquet schema directory, or converted from
    REGENIE TSVs, which is slower and uses more memory.

    Args:
        single_variants_file: Single-variant meta-analysis results (TSV). With an
            `AAF_BIN` column (`fedcat meta --aaf-bins`), each row is used only in its
            own bin; otherwise bins filter on A1FREQ.
        burden_file: Burden (MAC-10) meta-analysis results (TSV), if any.
        schema_dir: Annotation schema directory from `fedcat convert-annots`; the
            alternative to `masks_tsvs` and `annotations_tsvs`.
        masks_tsvs: REGENIE mask files.
        annotations_tsvs: REGENIE annotation files.
        output_file: Output TSV (`ACAT_OUTPUT_DTYPES` columns).
        aggregate_freqs: Upper allele-frequency bounds of the bins, e.g.
            `["0.001", "0.01"]`, used verbatim as bin labels; defaults to `["0.001"]`.
        id_col: ID column of the results files.
        allele_freq_col: Allele-frequency column of the results files.
        log10p_col: -log10(P) column of the results files.
        n_threads: Threads; None for one fewer than the CPUs.
        weighting_scheme: `"standard"` (GCTA/REGENIE weights) or `"equal"`.
        compression: Output compression, `"zst"` or `"gz"`, matching the file
            suffix; None to infer it from the suffix.
        compression_level: Compression level; None for the default.
        backend: `"sweep"`, for schemas sorted by unit start, or `"arrow"`.
        full_bounds: When converting TSVs, take each mask's START/END from all its
            annotated variants rather than only those with single-variant results,
            so masks with only a burden result get positions.

    Raises:
        ValueError: If the annotations are not given one way, files are missing, or
            an option is invalid.
        FileNotFoundError: If a schema file is missing.
    """
    if aggregate_freqs is None:
        aggregate_freqs = ["0.001"]
    if not aggregate_freqs:
        raise ValueError("aggregate_freqs must contain at least one threshold")
    max_freq = max(float(f) for f in aggregate_freqs)
    if n_threads is None:
        n_cpu = os.cpu_count()
        n_threads = n_cpu - 1 if n_cpu and n_cpu > 1 else 1
    if weighting_scheme not in {"standard", "equal"}:
        raise ValueError(f"Invalid weighting_scheme '{weighting_scheme}'; must be one of {'standard', 'equal'}")
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

    if schema_dir is not None:
        if not schema_dir.is_dir():
            raise ValueError(f"Schema directory does not exist or is not a directory: {schema_dir}")
        if annotations_tsvs is not None or masks_tsvs is not None:
            raise ValueError(
                "When schema_dir is provided, annotations_tsvs and masks_tsvs should not be provided"
            )
        required_files = [
            "unit_mask_bounds.parquet",
            "dim_variants.parquet",
            "dim_units.parquet",
            "dim_annots.parquet",
            "mask_annotations.parquet",
            "annotations.parquet",
        ]
        for filename in required_files:
            if not (schema_dir / filename).exists():
                raise FileNotFoundError(f"Missing schema file: {schema_dir / filename}")
    else:
        if annotations_tsvs is None or masks_tsvs is None:
            raise ValueError("When schema_dir is not provided, annotations_tsvs and masks_tsvs must be provided")
        elif not all(f.exists() for f in annotations_tsvs):
            missing = [str(f) for f in annotations_tsvs if not f.exists()]
            raise ValueError(f"One or more annotation files do not exist: {missing}")
        elif not all(f.exists() for f in masks_tsvs):
            missing = [str(f) for f in masks_tsvs if not f.exists()]
            raise ValueError(f"One or more masks files do not exist: {missing}")
        logger.warning(
            "TSV annotation input detected: converting to parquet on the fly. "
            "Intended usage is pre-converted parquet input "
            "(fedcat convert-annots --masks <masks1.tsv> [<masks2.tsv> ...] --annots <annot1.tsv> [<annot2.tsv> ...] --output-dir <dir>)."
        )

        present_singles = (
            pl.scan_csv(
                single_variants_file,
                separator="\t",
                comment_prefix="##",
                schema_overrides={id_col: pl.Utf8, allele_freq_col: pl.Float64},
            )
            .filter(pl.col(allele_freq_col) <= max_freq)
            .select(pl.col(id_col))
            .unique()
            .collect(engine="streaming")[id_col]
            .to_list()
        )

        schema_dir = Path(tempfile.mkdtemp(prefix="fedcat_acat_schema_"))
        conv_meta = _convert_tsv_to_parquet_rs(
            annotation_tsvs=[str(f) for f in annotations_tsvs],
            masks_tsvs=[str(f) for f in masks_tsvs],
            output_dir=str(schema_dir),
            valid_variant_ids=present_singles,
            full_bounds=full_bounds,
        )
        logger.info(
            f"[TSV -> parquet] Converted {conv_meta['n_rows']} rows, "
            f"{conv_meta['n_variants']} variants, "
            f"{conv_meta['n_units']} units, "
            f"{conv_meta['n_annots']} annots"
        )

    logger.info("Using Rust streaming ACAT processor (all I/O in Rust)")
    col_cfg = {
        "id": id_col,
        "a1freq": allele_freq_col,
        "log10p": log10p_col,
    }
    meta = _compute_acat_streaming_rs(
        singles_tsv=str(single_variants_file),
        burden_tsv=str(burden_file) if burden_file is not None else "",
        schema_dir=str(schema_dir),
        aggregate_freqs=aggregate_freqs,
        weighting_scheme=weighting_scheme,
        output_path=str(output_file) if output_file else "",
        column_config=col_cfg,
        compression=compression,
        compression_level=compression_level,
        backend=backend,
    )

    logger.info(
        f"ACAT complete: combined {meta['_n_singles']} single variants, {meta['_n_burden']} aggregate tests, "
        f"{meta['_n_annots']} annotations -> {meta['_n_results']} ACAT aggregates"
    )


def _load_masks(masks_tsvs: list[Path]) -> tuple[list[str], list[list[str]]]:
    """Read REGENIE mask files (`MASK<tab>annot1,annot2,...`), merging masks of the same name.

    Returns:
        `(mask_names, mask_annots)`: parallel lists of names and their annotations.

    Raises:
        ValueError: If a mask name contains a `.`.
    """
    mask_order: list[str] = []
    mask_annots: dict[str, list[str]] = {}
    for path in masks_tsvs:
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) < 2:
                    continue
                name = parts[0]
                if "." in name:
                    raise ValueError(f"mask name '{name}' contains a dot, which is not allowed")
                annots = [a.strip() for a in parts[1].split(",") if a.strip()]
                if name not in mask_annots:
                    mask_order.append(name)
                    mask_annots[name] = annots
                else:
                    existing = mask_annots[name]
                    for a in annots:
                        if a not in existing:
                            existing.append(a)
    return mask_order, [mask_annots[n] for n in mask_order]


def convert_tsv_to_parquet_polars(
    annotation_tsvs: list[Path],
    masks_tsvs: list[Path],
    output_dir: Path,
    valid_variant_ids: list[str] | None = None,
) -> dict[str, int]:
    """Convert REGENIE annotation files to the parquet schema in Python, as a reference for the Rust conversion.

    Args:
        annotation_tsvs: Uncompressed REGENIE annotation files (`variant unit annot`).
        masks_tsvs: REGENIE mask files.
        output_dir: Schema directory, created if needed.
        valid_variant_ids: Variant IDs to keep; None for all.

    Returns:
        Counts `n_rows`, `n_variants`, `n_units` and `n_annots`.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lfs = [
        pl.scan_csv(
            str(p),
            has_header=False,
            separator="\t",
            new_columns=["VARIANT_ID", "UNIT", "ANNOT"],
            schema_overrides={"VARIANT_ID": pl.Utf8, "UNIT": pl.Utf8, "ANNOT": pl.Utf8},
        )
        for p in annotation_tsvs
    ]
    lf = pl.concat(lfs) if len(lfs) > 1 else lfs[0]

    # Keep chr:pos:ref:alt IDs
    lf = lf.filter(
        (pl.col("VARIANT_ID").str.count_matches(":") == 3) & (pl.col("VARIANT_ID").str.len_chars() >= 7)
    )

    if valid_variant_ids is not None:
        valid_series = pl.Series("VARIANT_ID", valid_variant_ids, dtype=pl.Utf8)
        lf = lf.filter(pl.col("VARIANT_ID").is_in(valid_series))

    lf = lf.with_columns(
        pl.col("VARIANT_ID")
        .str.split(":")
        .list.get(0)
        .str.replace(r"^(?i)chr", "")
        .str.to_uppercase()
        .alias("CHR"),
        pl.col("VARIANT_ID").str.split(":").list.get(1).cast(pl.Int64).alias("POS"),
    )

    df = lf.collect()
    logger.info(f"Loaded {df.height} valid annotation rows from {len(annotation_tsvs)} TSV files")

    if df.height == 0:
        for name, cols in [
            ("dim_variants", {"VARIANT_CODE": pl.Int64, "VARIANT_ID": pl.Utf8}),
            ("dim_units", {"UNIT_CODE": pl.Int64, "UNIT": pl.Utf8}),
            ("dim_annots", {"ANNOT_CODE": pl.Int64, "ANNOT": pl.Utf8}),
            (
                "annotations",
                {
                    "VARIANT_CODE": pl.Int64,
                    "UNIT_CODE": pl.Int64,
                    "ANNOT_CODE": pl.Int64,
                    "CHR": pl.Utf8,
                    "POS": pl.Int64,
                },
            ),
            (
                "unit_mask_bounds",
                {"UNIT": pl.Utf8, "CHR": pl.Utf8, "MASK": pl.Utf8, "MIN_POS": pl.Int64, "MAX_POS": pl.Int64},
            ),
            ("mask_annotations", {"MASK": pl.Utf8, "ANNOT_CODE": pl.Int64, "ANNOT": pl.Utf8}),
        ]:
            pl.DataFrame(schema=cols).write_parquet(output_dir / f"{name}.parquet")
        return {"n_rows": 0, "n_variants": 0, "n_units": 0, "n_annots": 0}

    # Integer codes, in sorted order of the values
    dim_variants = (
        df.select("VARIANT_ID")
        .unique()
        .sort("VARIANT_ID")
        .with_row_index("VARIANT_CODE")
        .cast({"VARIANT_CODE": pl.Int64})
    )
    dim_units = df.select("UNIT").unique().sort("UNIT").with_row_index("UNIT_CODE").cast({"UNIT_CODE": pl.Int64})
    dim_annots = (
        df.select("ANNOT").unique().sort("ANNOT").with_row_index("ANNOT_CODE").cast({"ANNOT_CODE": pl.Int64})
    )

    df = (
        df.join(dim_variants, on="VARIANT_ID", how="left")
        .join(dim_units, on="UNIT", how="left")
        .join(dim_annots, on="ANNOT", how="left")
    )

    # Chromosome sort order: autosomes, X, Y, MT, then others
    chr_rank_map = {"X": 23, "23": 23, "Y": 24, "24": 24, "MT": 25, "M": 25, "25": 25}
    df = df.with_columns(
        pl.col("CHR")
        .replace(chr_rank_map, default=None)
        .cast(pl.UInt32, strict=False)
        .fill_null(pl.col("CHR").cast(pl.UInt32, strict=False).fill_null(26))
        .alias("CHR_RANK")
    )

    # Sort by unit start
    df = df.with_columns(pl.col("POS").min().over(["UNIT_CODE", "CHR_RANK"]).alias("UNIT_MIN_POS")).sort(
        ["CHR_RANK", "UNIT_MIN_POS", "UNIT_CODE", "POS"]
    )

    n_rows = df.height
    n_variants = dim_variants.height
    n_units = dim_units.height
    n_annots = dim_annots.height
    logger.info(f"Found {n_variants} unique variants, {n_units} units, {n_annots} annotations")

    dim_variants.write_parquet(output_dir / "dim_variants.parquet")
    dim_units.write_parquet(output_dir / "dim_units.parquet")
    dim_annots.write_parquet(output_dir / "dim_annots.parquet")

    df.select(["VARIANT_CODE", "UNIT_CODE", "ANNOT_CODE", "CHR", "POS"]).write_parquet(
        output_dir / "annotations.parquet"
    )

    # Each unit and mask's position bounds
    mask_names, mask_annot_lists = _load_masks(masks_tsvs)

    annot_to_masks: dict[str, list[int]] = defaultdict(list)
    for mi, annot_list in enumerate(mask_annot_lists):
        for a in annot_list:
            annot_to_masks[a].append(mi)

    annot_code_map: dict[str, int] = dict(
        zip(
            dim_annots["ANNOT"].to_list(),
            dim_annots["ANNOT_CODE"].to_list(),
        )
    )

    bounds_rows: list[dict] = []
    for (unit, chr_val), grp in df.group_by(["UNIT", "CHR"]):
        annots_in_grp = grp["ANNOT"].to_list()
        positions = grp["POS"].to_list()
        mask_positions: dict[int, list[int]] = defaultdict(list)
        for a, p in zip(annots_in_grp, positions):
            for mi in annot_to_masks.get(a, []):
                mask_positions[mi].append(p)
        for mi, pos_list in mask_positions.items():
            bounds_rows.append(
                {
                    "UNIT": unit,
                    "CHR": chr_val,
                    "MASK": mask_names[mi],
                    "MIN_POS": min(pos_list),
                    "MAX_POS": max(pos_list),
                }
            )

    if bounds_rows:
        pl.DataFrame(
            bounds_rows,
            schema={"UNIT": pl.Utf8, "CHR": pl.Utf8, "MASK": pl.Utf8, "MIN_POS": pl.Int64, "MAX_POS": pl.Int64},
        ).write_parquet(output_dir / "unit_mask_bounds.parquet")
    else:
        pl.DataFrame(
            schema={"UNIT": pl.Utf8, "CHR": pl.Utf8, "MASK": pl.Utf8, "MIN_POS": pl.Int64, "MAX_POS": pl.Int64}
        ).write_parquet(output_dir / "unit_mask_bounds.parquet")

    ma_rows = [
        {"MASK": mask_names[mi], "ANNOT_CODE": annot_code_map[a], "ANNOT": a}
        for mi, annot_list in enumerate(mask_annot_lists)
        for a in annot_list
        if a in annot_code_map
    ]
    pl.DataFrame(
        ma_rows,
        schema={"MASK": pl.Utf8, "ANNOT_CODE": pl.Int64, "ANNOT": pl.Utf8},
    ).write_parquet(output_dir / "mask_annotations.parquet")

    logger.info(
        f"Wrote optimised schema to {output_dir}: "
        f"{n_rows} rows, {n_variants} variants, {n_units} units, {n_annots} annots"
    )
    return {"n_rows": n_rows, "n_variants": n_variants, "n_units": n_units, "n_annots": n_annots}
