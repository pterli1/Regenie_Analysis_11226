# fedcat setup-workflows

Stage a biobank's WDL workflows and their scripts into the current Verily Workbench workspace, once per workspace and again for each new workflow version. Only AoU is implemented; UKB uses DNAnexus applets built from [`ukb_regenie_apps`](https://github.com/ExeterGenetics/biobank_workflows/tree/main/regenie_workflows/ukb_regenie_apps).

```bash
# Stage the default version, choosing the workspace bucket interactively
fedcat setup-workflows

# A given version and bucket, replacing existing workflows without asking
fedcat setup-workflows --version 1.7.1 --bucket gs://my-workspace-bucket --force

# See what it would do
fedcat setup-workflows --dry-run
```

## What it does

1. Uploads the scripts the WDL tasks copy at run time to `gs://<bucket>/workflows/<set>/v<version>/dependencies/`.
2. Sets each WDL's `dependencies_path` default to that folder, uploads the WDL to `.../v<version>/wdls/`, and registers it with `wb workflow create` as `<wdl>-v<version>`, with dots as dashes (e.g. `regenie_step1_scatter_aou-v1-7-1`).

Each version has its own names, so versions sit side by side and a run uses the one matching its cohort's `dependencies_version`. With the `wb_cli` backend, fedcat finds the registered workflow running that version's WDL, whatever its name.

!!! warning "Changing scripts during a run"

    Restaging a version overwrites its scripts in place, so tasks starting later in a running analysis get the new scripts. Use a new version number to change scripts mid-analysis.

## Options

--8<-- "options/setup-workflows.md"
