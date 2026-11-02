# Configuration

A run is configured by one JSON file, passed with `fedcat discover --config`. This page goes through it section by section; the [config file reference](../reference/config.md) lists every field with its type and default.

Only the options a run uses are required: e.g. a singles-only run needs no annotations, and a run without the meta-analysis phases no `meta_analysis_config`.

!!! tip "Starting points"

    - [`config_templates/`](https://github.com/ExeterGenetics/biobank_workflows/tree/main/regenie_workflows/federated_analysis/config_templates) has templates, including [`v9_all_options_with_defaults.json`](https://github.com/ExeterGenetics/biobank_workflows/blob/main/regenie_workflows/federated_analysis/config_templates/v9_all_options_with_defaults.json) with every option at its default.
    - [`configs/`](https://github.com/ExeterGenetics/biobank_workflows/tree/main/regenie_workflows/federated_analysis/configs) has configs from real runs.

The examples below have comments; a real config file is plain JSON and cannot.

## The run

```json
{
  // Environment fedcat itself runs in: "aou" on an AoU Jupyter instance, "local" on
  // another server or laptop. A cohort that forbids data egress (AoU) needs "aou".
  "local_env": "aou",

  // Phenotypes to run (up to 20 at once; see "Large phenotype sets" for more).
  // Every phenotype and covariate file must have a column for each name.
  "pheno_names": ["pheno1", "pheno2"],
  // "Binary" or "Continuous". Binary phenotypes are not currently tested or
  // recommended for the conditioning phases.
  "trait_type": "Continuous",
  // "RINT", "MCC" or "RAW". MCC is not recommended with the current REGENIE
  // version; binary phenotypes should use RAW.
  "trait_transform": "RINT",
  // Chromosomes to run
  "chromosomes_to_run": ["1", "2", "X"],
  // State files, logs and conditioning results go in subdirectories of this path,
  // which is created if it does not exist
  "local_run_path": "/home/jupyter/looping/run_folder",
  // Name of the run, used in its directory and output names
  "run_prefix": "run_identifier",

  // One entry per biobank; see "Cohorts" below
  "cohorts": {}
}
```

## Cohorts

`cohorts` has an entry per biobank to run: `aou`, `ukb` or both. They share a structure but differ in how workflows, resources and paths are given.

=== "AoU"

    ```json
    "aou": {
      // Workflow names in the Verily workspace (or WDL file names with local_cromwell)
      "workflows": {
        "regenie_step1": "regenie_step1_scatter_aou",
        "regenie_step2_singlevar": "regenie_step2_single_variants_aou",
        "regenie_step2_aggregate": "regenie_step2_aggregates_aou",
        "conditioning_singlevar": "regenie_conditioning_singles_aou",
        "conditioning_aggregate": "regenie_conditioning_agg_aou",
        // Optional machine resources; see "Resources" below
        "resources": {
          "conditioning_aggregate": {"n_cpu": 8, "mem_gb": 128},
          "regenie_step2_singlevar": {"disk": "200 SSD"},
          "regenie_step1": {"preemptible": 0, "max_retries": 1}
        }
      },
      // Google bucket paths; the gs:// prefix is optional
      "pheno_file": "gs://workspace-bucket/path/to/pheno.tsv.gz",
      "cov_file": "gs://workspace-bucket/path/to/covs.tsv.gz",
      // Optional file of variants to ALWAYS condition on
      "cond_file": "gs://workspace-bucket/path/to/cond_variants",
      "strata": ["afr", "amr", "eas", "eur", "mid", "sas"],
      "cat_vars": "seq_center",
      "step1_directory": "regenie_step1_output",
      "root_results_directory": "wgs_analysis2/vYYYY.MM",
      "pfile_template": "gs://wgs-bucket/path_to_pgens/chr{chr}_pgen_name",
      "step1_input_template": "gs://wgs-bucket/path_to_regenie_step1_input/regenie_step1_{stratum}",
      "regenie_annotations": [
        {"name": "coding", "mask_file": "gs://...", "annot_file": "gs://...", "setlist_file": "gs://..."}
      ],
      "cohort_sizes": {"afr": 84000, "amr": 79000, "eas": 10000, "eur": 234000, "mid": 1500, "sas": 5500},
      // Version of the workflows staged by `fedcat setup-workflows` (default 1.7.1)
      "dependencies_version": "1.7.1",
      // "wb_cli" (default) or "local_cromwell"
      "workflow_backend": "wb_cli"
    }
    ```

=== "UKB"

    ```json
    "ukb": {
      // DNAnexus applet IDs
      "workflows": {
        "regenie_step1": "applet-xxxx",
        "regenie_step2_singlevar": "applet-xxxx",
        "regenie_step2_aggregate": "applet-xxxx",
        "conditioning_singlevar": "applet-xxxx",
        "conditioning_aggregate": "applet-xxxx",
        // Optional; see "Resources" below
        "resources": {
          "conditioning_aggregate": {"instance_type": "mem3_ssd1_v2_x32"},
          "regenie_step1": {"priority": "high"}
        }
      },
      // DNAnexus file IDs, or paths including the project (project-xxx:/path/to/file)
      "pheno_file": "file-xxxx",
      "cov_file": "file-xxxx",
      "cond_file": "file-xxxx",
      "strata": ["afr", "eur", "sas"],
      "cat_vars": "centre,sequencing_centre",
      "step1_directory": "regenie_step1_output",
      "root_results_directory": "wgs_analysis2/v2026.2",
      "pfile_template": "project-xxx:/500k_wgs_pgens/chr{chr}_trimmed_merged_nofail",
      "step1_input_template": "project-xxx:/regenie_step1_input/regenie_step1_{stratum}",
      // File IDs or paths; they are passed to the applets as paths, so every file
      // must have a unique path even if given as an ID
      "regenie_annotations": [
        {"name": "coding", "mask_file": "project-xxx:/...", "annot_file": "project-xxx:/...", "setlist_file": "project-xxx:/..."}
      ],
      "cohort_sizes": {"afr": 10000, "eur": 460000, "sas": 10000},
      // Optional: project this cohort's jobs run and write in (defaults to the
      // top-level dx_project_id), and the project holding the applets
      "dx_project_id": "project-xxxx",
      "applet_project_id": "project-xxxx"
    }
    ```

### Input files

| Field | |
|---|---|
| `pheno_file`, `cov_file` | The cohort's phenotype and covariate files, holding the individuals of every stratum |
| `cat_vars` | Comma-separated categorical covariates, or `NONE` |
| `cond_file` | Optional variants to always condition on, set for every cohort or none. Step 2 then also tests single variants conditional on them |
| `pfile_template` | PLINK2 pgen file set per chromosome, without the `.pgen`/`.pvar`/`.psam` extension, with a `{chr}` placeholder |
| `step1_input_template` | REGENIE step 1 input file set per stratum, without extension, with a `{stratum}` placeholder |

### Strata

`strata` are the sample groups run and meta-analysed separately, e.g. genetic ancestry groups. Each stratum's samples are those in its step 1 input file. Names may contain only letters, digits and `-`.

`cohort_sizes` gives the individuals per stratum, used to size step 2 jobs. Accurate values reduce costs; missing strata default to 100,000.

### Output locations

`step1_directory` and `root_results_directory` are relative to the workspace bucket (AoU) or DNAnexus project (UKB). Results go in a directory under `root_results_directory` named after `run_prefix` and `trait_transform`. See [Outputs](outputs.md).

### Annotations

`regenie_annotations` lists REGENIE annotation file triplets for aggregate tests (`--mask-def`, `--anno-file` and `--set-list`), each named to label its results. Paths may contain a `{chr}` placeholder and a `.zst` extension. Only the `step2_aggregates` and `conditioning_aggregates` phases need them; omit the field for a singles-only run.

### Resources

`workflows.resources` sets the machine for a kind of job, keyed by the workflow names (`regenie_step1`, ..., `conditioning_aggregate`) plus `meta_workflow` and `subset_workflow`. Others keep their defaults.

| Biobank | Fields |
|---|---|
| AoU | `n_cpu`, `mem_gb`, `disk` (`"<GB> SSD"` or `"<GB> HDD"`), `max_retries` (retries after a task fails; default 0), `preemptible` (attempts on a spot VM before a standard one; 0 never uses spot) |
| UKB | `instance_type` (a DNAnexus instance type), `priority` (`"high"` for on-demand machines; `"normal"` or `"low"` for cheaper spot machines, where a job fails if its machine is reclaimed) |

If jobs run out of memory, raise their resources and resume with `--overwrite-config --force-job-retry` to resubmit them (see [Troubleshooting](troubleshooting.md)).

### Workflow versions and backends (AoU)

`dependencies_version` must match a version staged by [`fedcat setup-workflows`](../commands/setup-workflows.md). With the `wb_cli` backend, fedcat finds the registered workflow running that version's WDL, whatever its name.

`workflow_backend` is `wb_cli` (the default: submit with the Workbench CLI) or `local_cromwell` (submit to a Cromwell server on this machine; see [`fedcat cromwell`](../commands/cromwell.md)).

## Meta-analysis

`meta_analysis_config` is needed by the `meta_analysis` and `step0_subset` phases.

```json
"meta_analysis_config": {
  // Where the meta-analysis and step-0 subset jobs run: "aou" for AoU, "local"
  // for a DNAnexus applet
  "location": "aou",
  // Workspace workflow name / WDL name, or DNAnexus applet ID (needed for meta_analysis)
  "meta_workflow": "meta_analysis_discovery_aou",
  // Workflow name / WDL name, or the regenie_step0_subset applet ID (needed for step0_subset)
  "subset_workflow": "step0_subsetting_only_aou",

  // ONE of the two annotation inputs below, needed only when the meta-analysis runs
  // its burden/ACAT stages (--meta-mode all or aggregates_only)
  // A: a parquet schema directory written by `fedcat convert-annots`
  "schema_directory": "gs://workspace-bucket/wgs_annotations/v2026.2/parquet_all_chr{chr}",
  // B: REGENIE-style annotation triplets, as for the cohorts
  "regenie_annotations": [
    {"name": "coding", "mask_file": "...", "annot_file": "...", "setlist_file": "..."}
  ]
}
```

See [Meta-analysis and ACAT](meta-analysis.md).

## Covariates

### Excluding covariates per chromosome

All optional, and combined per chromosome into REGENIE's `--covarExcludeList`.

```json
// LOCO PGS covariates: column {pgs_prefix}{chrom} is excluded on each chromosome,
// so step 2 for chr1 drops XX_SCORE_CHR1. Omit if you have no LOCO PGS covariates.
"pgs_prefix": "XX_SCORE_CHR",
// Chromosomes with no PGS covariate to exclude
"chromosomes_no_pgs": ["X"],
// Covariates excluded on every chromosome
"exclude_covars": ["irrelevant_covar1", "irrelevant_covar2"],
// Further exclusions per chromosome
"exclude_covars_dict": {"X": ["X_PGS"], "6": ["HLA_score2"]}
```

### Stratifying by a covariate

With `stratify_by_covar`, `fedcat discover` runs a complete, independent analysis for each listed value of a covariate (e.g. sex or a disease state) instead of one over everyone.

```json
"stratify_by_covar": {
  // Column of each cohort's cov_file to stratify on
  "covar_name": "sex",
  // Values to analyse, compared as strings exactly as they appear in the file,
  // each mapped to a covariate stratum name. Any number of values.
  "covar_values_names": {"2": "female", "1": "male"},
  // Also run the unstratified analysis (default true)
  "also_run_combined": true
}
```

Each covariate stratum gets its own run prefix, `{run_prefix}_{name}` (e.g. `{run_prefix}_female`), and so its own step 1 fit, outputs, run directory and state files. The combined run keeps `run_prefix`. Within a covariate stratum the stratifying covariate is dropped from the REGENIE covariates.

All covariate strata run in parallel by default; `--covar-stratum` restricts a run, e.g. `--covar-stratum female combined`.

## Large phenotype sets

For many phenotypes (e.g. proteomics), split up the work:

```json
// Phenotypes fitted together per step 1 job (default 0: all in one fit). Each set
// writes to "<step1_run_prefix>_set<N>" when there is more than one.
"step1_pheno_group_size": 20,
// Phenotypes tested together per step 2 task (default 0: all together). A group reads
// the predictions of every step 1 set holding its phenotypes, so it can be larger
// than a step 1 set; outputs are tagged "_g<N>" when grouped.
"step2_pheno_group_size": 250
```

and in the AoU cohort, split phases over several workflows and cap the load:

```json
// Split each step1/step2 phase over workflows of at most this many tasks (default 0: one workflow)
"max_tasks_per_workflow": 5000,
// Most jobs (AoU workflows / UKB applet jobs) running at once (default 0: no limit)
"max_concurrent_jobs": 10
```

## Testing options

```json
// Prefix of step 1 outputs (default: run_prefix). Set it to reuse a step 1 fit for
// several step 2 runs with different run_prefix (e.g. annotation versions).
"step1_run_prefix": "generic_phenotype_group",
// REGENIE version every step1/step2/conditioning job runs, on both biobanks
// (default "4.1.3"; also "4.1")
"regenie_version": "4.1.3",
// Minimum minor allele count for testing single variants (default 11)
"singles_mac_threshold": 11,
// Aggregate allele-frequency bins, passed verbatim to REGENIE --aaf-bins (default ["0.001"])
"aggregate_freqs": ["0.001"],
// Ceilings on the step 2 block length in bp, which otherwise grows as cohorts shrink
// (default: none for singles; 20 Mb for WDL aggregates; 30 Mb for DNAnexus aggregates)
"max_step2_singles_block_length": 50000000,
"max_step2_aggregates_block_length": 40000000,
// Restrict aggregate tests to sets starting with "coding_" (default false). Only
// recommended with old annotation files where set types are mixed.
"coding_only": false,
// GxE interaction testing (REGENIE --interaction; default off). Names a covariate
// column; use "VAR[BASE]" to set the baseline of a categorical variable, and list it
// in cat_vars. Only implemented for step 2 single variants.
"interaction": "bmi"
```

The ACAT options `acat_full_bounds` and `acat_singles_stratum_rare` are in [Meta-analysis and ACAT](meta-analysis.md#acat-options).

## Conditioning

The conditioning phases need a `cond_suffix`, added to the conditioning run directory names, and take their thresholds from the config:

```json
"cond_suffix": "run1",
"cond_p_singles_retest": 1e-6,
"cond_p_singles_select": 5e-8,
"cond_p_agg_retest": 1e-6,
"cond_p_agg_select": 1e-7,
"cond_rsq_threshold": 0.9,
"cond_check_coll": true,
"cond_max_steps": -1,
"cond_engine": "regenie"
```

See [Conditioning](conditioning.md) for what each does.

## Projects and buckets

```json
// DNAnexus project where UKB jobs run and results are written
// (default: the main UKB_500k_WGS project)
"dx_project_id": "project-xxxx",
// GCS bucket AoU results are written to. If unset, fedcat lists the workspace's
// buckets for you to choose from and saves the choice in the state files.
"aou_bucket": "workspace-bucket",
// Workbench resource id of that bucket (for `wb --bucket-id`). Resolved
// automatically; only set it to override.
"aou_bucket_resource_id": "workspace-bucket-resource-id"
```

## Changing the config of a run

The config is saved in the run's state files. To resume after editing it, pass `--overwrite-config` to replace the saved values that changed.
