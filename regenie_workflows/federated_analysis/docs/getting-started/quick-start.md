# Quick start

A first run: a single-variant GWAS of two phenotypes in AoU, from REGENIE step 1 to the meta-analysis across strata.

## 1. Write a config

Start from a template in [`config_templates/`](https://github.com/ExeterGenetics/biobank_workflows/tree/main/regenie_workflows/federated_analysis/config_templates). A minimal AoU config (paths shortened):

```json
{
  "local_env": "aou",
  "cohorts": {
    "aou": {
      "workflows": {
        "regenie_step1": "regenie_step1_scatter_aou",
        "regenie_step2_singlevar": "regenie_step2_single_variants_aou",
        "regenie_step2_aggregate": "regenie_step2_aggregates_aou",
        "conditioning_singlevar": "regenie_conditioning_singles_aou",
        "conditioning_aggregate": "regenie_conditioning_agg_aou"
      },
      "pheno_file": "gs://workspace-bucket/path/to/pheno.tsv.gz",
      "cov_file": "gs://workspace-bucket/path/to/covs.tsv.gz",
      "strata": ["afr", "amr", "eas", "eur", "mid", "sas"],
      "cat_vars": "seq_center",
      "step1_directory": "regenie_step1_output",
      "root_results_directory": "wgs_analysis",
      "pfile_template": "gs://wgs-bucket/pgens/chr{chr}_pgen_name",
      "step1_input_template": "gs://wgs-bucket/step1_input/regenie_step1_{stratum}"
    }
  },
  "meta_analysis_config": {
    "location": "aou",
    "meta_workflow": "meta_analysis_discovery_aou",
    "subset_workflow": "step0_subsetting_only_aou"
  },
  "pheno_names": ["pheno1", "pheno2"],
  "trait_type": "Continuous",
  "trait_transform": "RINT",
  "chromosomes_to_run": ["1", "2", "X"],
  "local_run_path": "/home/jupyter/fedcat_runs/my_first_run",
  "run_prefix": "my_first_run"
}
```

Only the options a run uses are required: with no aggregate tests this needs no annotations, and with no conditioning no `cond_suffix`. [Configuration](../guide/configuration.md) explains every section.

## 2. Check it with a dry run

```bash
fedcat discover --config config.json --phases-to-run step1 step2_singles meta_analysis --dry-run
```

This checks the config, your authentication, the workflows and the input files, then builds the first phase's jobs and prints a summary without submitting anything.

## 3. Run it

```bash
fedcat discover --config config.json --phases-to-run step1 step2_singles meta_analysis
```

FedCAT submits each phase's jobs, polls them, and moves on when they finish. Progress is reported as `x/y tasks completed (a running, b failed)`.

Stop it at any time with ++ctrl+c++: launched jobs keep running, and the same command resumes the run. See [Running and resuming](../guide/running.md).

!!! tip "Long runs"

    Run FedCAT inside `tmux` or `screen` so it survives a disconnected terminal.

## 4. Find the results

- Local state files and logs: `<local_run_path>/<run_prefix>_discovery_run/`.
- Per-cohort REGENIE results and the meta-analysis: under `root_results_directory` in the workspace bucket (AoU) or DNAnexus project (UKB).

See [Outputs](../guide/outputs.md) for the layout.

## Next steps

- Aggregate tests: add `regenie_annotations` to the cohort and the `step2_aggregates` phase ([Configuration](../guide/configuration.md#annotations)).
- Independent signals: add `cond_suffix` and the `step0_subset` and conditioning phases ([Conditioning](../guide/conditioning.md)).
- UKB: add a `ukb` cohort with DNAnexus applet IDs ([Configuration](../guide/configuration.md#cohorts)).
