# fedcat: deferred work

Points raised in code review and left for later, with the reason. Decided items are at the end so they aren't reopened.

## Analysis

- **Allele harmonisation in the meta-analysis.** Cohorts are combined on exact variant ID; alleles are not checked. Wanted as options (check, warn, drop) rather than automatic flipping.
- **ACAT of non-additive burden tests.** ACAT combines only burden results with TEST `ADD`, the only test the MAC-10 burden runs. A non-additive MAC-10 burden (e.g. a dominant model) would need ACAT to take the test as an option.
- **Cohorts with mixed `chr`-prefixed IDs.** The meta-analysis only warns, as such variants are never combined. The alternative is normalising IDs when reading results, while conditioning keeps IDs as in each pvar so jobs can find them.

## Platforms

- **DNAnexus job polling in batches.** Each job is polled on its own. Fine for now with DNAnexus's 100 concurrent jobs and the polling backoff; revisit once DNAnexus is back.
- **Re-staging a dependencies version mid-run.** `setup-workflows` overwrites that version's scripts in place, so tasks that start later in a running analysis pick up the new scripts. Bump the dependencies version when changing scripts during a run.

## Code health

- **Remaining bugbear rules.** `ruff check --extend-select B` reports 22 `zip()` calls without `strict=` (B905; most should be `strict=True`, each needs a look), 3 unused loop variables (B007) and 1 `raise` without `from` (B904, benchmark). Enable with `extend-select = ["B"]` in `[tool.ruff.lint]`.
- **Ruff's default rules change between versions.** Ruff is pinned; when upgrading, consider pinning `select` to the current set.
- **Dead code.** `convert_tsv_to_parquet_polars` and `_load_masks` in `stats/acat.py` are only used by `tests/benchmark_pipeline.py`; `z_to_log10p(use_rust=False)` and `log10p_to_p` in `stats/meta_analysis.py` only by tests.

## Testing

- **Untested:** the CLI glue in `__main__.py` (about 60% covered), the DNAnexus and Cromwell client wrappers (need the real services), and the reset flags on resuming conditioning (`--cond-reset-step`, `--reset-jobs`, `--force-reset` in `DiscoveryRun._run_conditioning`).
- **Coverage tooling.** `coverage`/`pytest-cov` are not in the dev extras. CI with test and coverage badges was left out while the repo is private.

## Decided, not changing

- The conditioning set for a selected aggregate has no frequency cap: conditioning on all its signal, including common causal variants, is intended.
- A burden result counts as empty when it has 2 lines: regenie writes all masks on one `##MASKS` line.
- Ruff's timezone rules (DTZ) are ignored: fedcat uses local time on purpose.
