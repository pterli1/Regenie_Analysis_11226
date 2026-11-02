# Outputs

## Discovery results in the cloud

Step 1, step 2 and meta-analysis results stay in the cloud, under each cohort's `root_results_directory`:

- AoU: `gs://<aou_bucket>/<root_results_directory>/<run_prefix>_<transform>_wgs/`
- UKB: `<dx_project_id>:/<root_results_directory>/<run_prefix>_<transform>_wgs/`

`<transform>` is `rint`, `mcc`, `raw` (continuous traits) or `binary`. The meta-analysis is in `meta_analysis/all_biobank_meta/` in the results directory of the biobank it ran in. Step 1 predictions go to `step1_directory`.

## The local run directory

`fedcat discover` keeps its state locally, in `<local_run_path>/<run_prefix>_discovery_run/`:

| Path | Contents |
|---|---|
| `discovery_state_<phase>.json` | Each phase's state, reloaded on resume |
| `logs/` | Log files per phase |

With [`stratify_by_covar`](configuration.md#stratifying-by-a-covariate), each covariate stratum has its own run directory, `<run_prefix>_<name>_discovery_run/`.

## Conditioning

Conditioning results are written locally; the per-cohort jobs' intermediate output stays in the cloud. Each phenotype and conditioning type has a directory in `local_run_path`:

```
<run_prefix>_<pheno>_<singlevar|aggregate>_conditioning_<cond_suffix>/
├── singles_cond_step0/        # one directory per step (agg_cond_step<N> for aggregates)
│   ├── conditioning_state_singles_step0.json
│   ├── tophits_input
│   └── summary_output
├── singles_cond_step1/
│   └── ...
└── singles_cond_results/      # written when conditioning finishes
    ├── all_sig_conditioned_singles_<pheno>.tsv
    ├── all_top_sig_conditioned_singles_<pheno>.tsv
    └── trace_chr<chrom>_<id>.tsv
```

Each step directory holds the step's intermediate files and state file (also uploaded to the bucket for AoU runs). The results directory has:

`all_sig_conditioned_<type>_<pheno>.tsv`
:   Every test of each selected signal, at the last step it was tested.

`all_top_sig_conditioned_<type>_<pheno>.tsv`
:   The most significant test of each selected signal.

`trace_chr<chrom>_<id>.tsv`
:   One per selected signal, tracing its statistics across the steps.

If conditioning stopped at `cond_max_steps` with signals left, results go to `<type>_cond_partial_results_step<N>/` instead. Results are also uploaded to the cloud, and [`fedcat summarise-indep`](../commands/summarise-indep.md) can rewrite them for a finished run.

!!! note "Runs from fedcat 0.6.0 and earlier"

    These named conditioning directories without the leading `<run_prefix>_`; add it to resume such a run.
