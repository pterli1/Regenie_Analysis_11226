# fedcat acat

Compute ACAT p-values for every unit, mask and allele-frequency bin from single-variant and, optionally, burden meta-analysis results. See [Meta-analysis and ACAT](../guide/meta-analysis.md#acat) for the method.

```bash
# With a parquet annotation schema from `fedcat convert-annots` (fastest)
fedcat acat --single-variants-file sv.tsv.zst --burden-file burden.tsv.zst \
    --output-file acat.tsv.zst --schema-dir schema/

# With REGENIE annotation and mask files, converted on the fly
fedcat acat --single-variants-file sv.tsv.zst --output-file acat.tsv.zst \
    --masks-tsvs masks.tsv --annotations-tsvs annot_chr1.tsv.gz \
    --aggregate-freqs 0.001,0.01
```

## Notes

- **Annotations** come from `--schema-dir`, or `--masks-tsvs` with `--annotations-tsvs`. TSVs are converted on every run, which is slower and uses more memory; to run ACAT repeatedly, convert them once with [`fedcat convert-annots`](convert-annots.md).
- **Frequency bins.** `--aggregate-freqs` sets the bins, used verbatim as labels, so match REGENIE's (`0.001`, not `1e-3`). With an `AAF_BIN` column (from `fedcat meta --aaf-bins`), each row is used only in its own bin; otherwise rows are binned on A1FREQ.
- **Weights.** `standard` uses REGENIE's allele-frequency weights; `equal` weights variants equally.
- **Mask positions.** With TSV annotations, mask START/END come from all annotated variants; `--no-full-bounds` uses only those with single-variant results, which is faster but leaves burden-only masks without positions.
- **Backends.** `sweep` (the default) needs a schema sorted by unit start, as `convert-annots` writes; `arrow` needs no sorting but allows at most 64 annotations and masks.

## Options

--8<-- "options/acat.md"
