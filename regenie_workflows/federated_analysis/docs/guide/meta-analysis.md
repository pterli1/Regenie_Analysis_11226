# Meta-analysis and ACAT

The `meta_analysis` phase combines the REGENIE results of every cohort and stratum, then tests each aggregate by ACAT. Both are also standalone commands, [`fedcat meta`](../commands/meta.md) and [`fedcat acat`](../commands/acat.md).

## Where it runs

`meta_analysis_config.location` sets where the meta-analysis job runs:

- `aou`: a workflow in the AoU workspace (`meta_workflow`). UKB summary statistics are first copied into the AoU bucket.
- `local`: a DNAnexus applet (`meta_workflow` is then an applet ID).

Results go to `meta_analysis/all_biobank_meta/` in that biobank's [results directory](outputs.md).

## Meta-analysis

Results are combined by **inverse-variance weighted (IVW) fixed-effects meta-analysis**, with Cochran's Q, its p-value and I² for heterogeneity. P-values are kept as -log10(P), so extreme values do not underflow.

- Variants are matched on chromosome, position, ID and test. Alleles are not checked, so cohorts must use the same variant IDs.
- Rows with SE ≤ 0, a missing or non-finite value, or a REGENIE `TEST_FAIL` are dropped; chromosome 23 becomes X.
- `A1FREQ` is N-weighted across cohorts, and `DIRECTIONS` shows each cohort's effect direction (`+`, `-`, or `?` where it has no result).

The output is a tab-separated table after a `## Cohorts:` line, with columns `CHROM POS ID A0 A1 A1FREQ N TEST BETA SE Z LOG10P Q Q_P I2 N_STUDIES DIRECTIONS`.

### What is meta-analysed

`--meta-mode` on `fedcat discover` chooses the result types:

| Mode | Meta-analyses |
|---|---|
| `auto` (default) | `singles_only` if only singles phases (step 2 or conditioning) are in `--phases-to-run`, `aggregates_only` if only aggregate phases are, otherwise `all` |
| `all` | Single variants, burden tests, MAC-10 singletons and ACAT |
| `singles_only` | Single variants (unconditional, and conditional on the `cond_file` if set) |
| `cond_singles_only` | Only single variants conditional on the `cond_file` |
| `uncond_singles_only` | Only unconditional single variants |
| `aggregates_only` | Burden tests, MAC-10 singletons and ACAT |

The burden and ACAT stages need annotations: either `schema_directory` (a parquet schema from [`fedcat convert-annots`](../commands/convert-annots.md)) or `regenie_annotations` in `meta_analysis_config`.

## ACAT

The aggregated Cauchy association test (ACAT) gives each aggregate (unit, e.g. a gene; mask; and allele-frequency bin) one p-value, combining its single-variant results in the bin with its burden result. Variants are weighted by allele frequency, as in REGENIE, or equally with `fedcat acat --weighting-scheme equal`. Results have test `ADD` and ID `UNIT.MASK.BIN`.

Single-variant results are meta-analysed per frequency bin for ACAT: a cohort contributes to bin `F` only where its own `A1FREQ` ≤ `F`.

### ACAT options

`acat_full_bounds` (default `true`)
:   Take each mask's START/END from all its annotated variants, not only those with single-variant results, so burden-only masks get positions too. Makes ACAT 30–40% slower.

`acat_singles_stratum_rare` (default `true`)
:   A stratum's single-variant result enters a bin's ACAT only where the variant is rare (`A1FREQ` ≤ bin) in that stratum, rather than filtering on the meta-analysed `A1FREQ`. Fixed for a run: changing it means rerunning from the meta-analysis.

## The step-0 subset

Results significant at `cond_p_singles_retest` and `cond_p_agg_retest` are packaged as the starting point for [conditioning](conditioning.md): per chromosome, a singles tarball and an aggregate tarball with the significant units' annotations. The meta-analysis job builds them while its outputs are on disk; the `step0_subset` phase reuses them and builds only new thresholds. Subsets are stored per threshold (e.g. `step0_subsets/singles_p1e-06`), so thresholds never mix.
