# fedcat convert-annots

Convert REGENIE annotation and mask files to the parquet schema that [`fedcat acat`](acat.md) and the meta-analysis read fastest. Run it once per annotation set and point `--schema-dir` (or `meta_analysis_config.schema_directory`) at the output.

```bash
fedcat convert-annots --masks masks.tsv --annots annot_chr*.tsv.gz --output-dir schema/
```

## Inputs

- **Mask files** (`--masks`): REGENIE `--mask-def` files, `MASK<tab>annot1,annot2,...`.
- **Annotation files** (`--annots`): REGENIE `--anno-file` files, `VARIANT<tab>UNIT<tab>ANNOT`, plain, gzip or zstd.

Rows whose variant ID is not `chr:pos:ref:alt` are skipped.

## The schema

| File | Contents |
|---|---|
| `dim_variants.parquet`, `dim_units.parquet`, `dim_annots.parquet` | Integer codes of variant IDs, units (genes or regions) and annotations |
| `annotations.parquet` | Each variant's unit and annotation, with its chromosome and position, sorted by unit start |
| `mask_annotations.parquet` | The annotations in each mask |
| `unit_mask_bounds.parquet` | Each unit and mask's first and last position |

## Options

--8<-- "options/convert-annots.md"
