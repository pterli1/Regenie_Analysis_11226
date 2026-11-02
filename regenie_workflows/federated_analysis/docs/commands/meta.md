# fedcat meta

Meta-analyse REGENIE results from several cohorts by inverse-variance weighted fixed effects, as `fedcat discover` does in its meta-analysis and conditioning jobs. See [Meta-analysis and ACAT](../guide/meta-analysis.md#meta-analysis) for the method and output columns.

```bash
fedcat meta --cohort-file cohorts.txt --output-file results.tsv.zst --res-type single_variants
```

## The cohort file

Each `# <name>` line starts a cohort, followed by its REGENIE result files (e.g. one per chromosome):

```text
# eur_aou
/path/to/eur_aou_chr1.regenie.gz
/path/to/eur_aou_chr2.regenie.gz
# eur_ukb
/path/to/eur_ukb_chr1.regenie.gz
/path/to/eur_ukb_chr2.regenie.gz
```

`.regenie.gz` files are read directly. If no files are listed, no output is written.

## Notes

- **Backends.** `streaming` (the default) adds one cohort at a time. `merge` merges all files at once; it needs input sorted as REGENIE writes it, and memory grows with the number of files.
- **Frequency bins.** With `--aaf-bins`, a cohort contributes to bin `F` only where its own A1FREQ ≤ `F`, and the output has a row per variant and bin (`AAF_BIN` column): the single-variant input [`fedcat acat`](acat.md) expects.
- **Other input formats.** The `--*-col` options name non-REGENIE columns.
- **Compression** is inferred from the output file name (`.zst` or `.gz`) unless `--compression` is given.

## Options

--8<-- "options/meta.md"
