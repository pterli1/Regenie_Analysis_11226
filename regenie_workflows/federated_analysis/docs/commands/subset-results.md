# fedcat subset-results

Build the step-0 subsets that conditioning starts from: a phenotype's significant meta-analysis results, packaged as tarballs. `fedcat discover` does this in the cloud; this command does it on local files. See [The step-0 subset](../guide/meta-analysis.md#the-step-0-subset).

```bash
fedcat subset-results --phenotypes BMI HEIGHT --chr 1 \
    --singles-file singles.tsv.zst \
    --acat-file acat.tsv.zst --burden-file burden.tsv.zst --schema-dir schema/ \
    --singles-p-threshold 1e-6 --agg-p-threshold 1e-6
```

## Outputs

For each phenotype, in the working directory:

- `<pheno>_singles_cond_step0_subset[_chr<chr>].tar.gz`, from `--singles-file`.
- `<pheno>_agg_cond_step0_subset[_chr<chr>].tar.gz`, from `--acat-file` and `--burden-file`, with the significant units' annotations (a parquet schema, or REGENIE annotation and mask files).

Start conditioning from one with `fedcat discover --cond-singles-tarball` or `--cond-agg-tarball`.

## Notes

- **Thresholds.** `--singles-p-threshold` and `--agg-p-threshold` correspond to `cond_p_singles_retest` and `cond_p_agg_retest`; without them, `--log10p-threshold` applies to both.
- **Annotations** (`--schema-dir`, or `--annotations` with `--masks-tsvs`) are needed only for the aggregate subset.
- **Burden results** take their positions from their unit's ACAT results, so `--burden-file` needs `--acat-file`.

## Options

--8<-- "options/subset-results.md"
