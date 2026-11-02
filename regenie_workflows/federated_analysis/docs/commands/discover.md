# fedcat discover

Run a federated analysis: REGENIE step 1 and step 2 in every cohort, the meta-analysis, and conditioning. The analysis is set by the [config file](../guide/configuration.md); these options control which phases run and how a run is resumed, inspected and repaired.

## Examples

```bash
# Run all phases from the start, or resume from the latest phase
fedcat discover --config config.json

# Run a subset of phases
fedcat discover --config config.json --phases-to-run step2_singles step2_aggregates meta_analysis

# Start at a specific phase
fedcat discover --config config.json --start-phase step2_singles

# Start conditioning from a meta-analysis in the cloud, or after changing cond_p_*_retest
fedcat discover --config config.json --phases-to-run step0_subset conditioning_singles conditioning_aggregates

# Build jobs and print a summary without submitting anything
fedcat discover --config config.json --dry-run

# Restrict REGENIE step 2 to some regions (whitespace-separated CHROM START END)
fedcat discover --config config.json --run-blockfile regions.txt

# Resubmit failed jobs after raising their resources in the config
fedcat discover --config config.json --overwrite-config --force-job-retry

# Job states, task counts and costs of the whole run
fedcat discover --config config.json --report
```

## Guides

- [How a run works](../guide/how-it-works.md): the phases and choosing which to run.
- [Running and resuming](../guide/running.md): resuming, retrying failed jobs, inspecting and aborting a run.
- [Conditioning](../guide/conditioning.md): the conditioning phases and their options.

## Options

--8<-- "options/discover.md"

## State codes

`--reset-jobs` takes a `JobState` value, and `--cond-reset-step` a `StepState` value:

--8<-- "state-codes.md"
