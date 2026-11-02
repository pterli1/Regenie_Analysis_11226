# Installation

FedCAT runs from the terminal of an AoU Verily Workbench Jupyter instance. It submits AoU jobs as workspace workflows and UKB jobs as DNAnexus applets.

## Requirements

- Python 3.10 or later (already present on AoU).
- **AoU authentication:** `fedcat` uses Application Default Credentials and the `wb` command-line tool, both normally set up already on AoU.
- **DNAnexus authentication** (for UKB cohorts): run `dx login`.

## 1. Clone the repository

Clone this repository onto the Jupyter instance, most easily with the [Verily GitHub integration](https://support.workbench.verily.com/docs/guides/cloud_apps/git_integrations/).

??? note "Cloning with a personal access token instead"

    1. Go to <https://github.com/settings/personal-access-tokens> and create a token with `ExeterGenetics` as resource owner.
    2. Choose *All repositories* (or select just this one) and add the *Contents* permission (leave it read-only).
    3. Run `git clone https://github.com/ExeterGenetics/biobank_workflows.git` and paste the token when asked for a password.

## 2. Install FedCAT

The install script creates a working directory (`fedcat_runs`, or its first argument) with a virtual environment holding FedCAT and its Rust extension, `fedcat_rs`.

```bash
bash repos/biobank_workflows/regenie_workflows/create_working_directory.sh
cd fedcat_runs
source venv/bin/activate
```

!!! note "The fedcat_rs extension"

    The script installs a prebuilt `fedcat_rs` wheel (Linux x86_64), or builds it from source with Rust. If both fail, run `pip install ./fedcat_rs` in `regenie_workflows/federated_analysis`, then install FedCAT. After pulling, reinstall `fedcat_rs` if FedCAT needs a newer version.

## 3. Stage the workflows

Add the AoU workflows and their scripts to your Workbench workspace, once per workspace and again for each new workflow version:

```bash
fedcat setup-workflows
```

Each workflow is registered under a versioned name, so versions sit side by side. The version must match your config's `dependencies_version`. See [`fedcat setup-workflows`](../commands/setup-workflows.md).

UKB runs use DNAnexus applets, built from [`ukb_regenie_apps`](https://github.com/ExeterGenetics/biobank_workflows/tree/main/regenie_workflows/ukb_regenie_apps); their IDs go in the config.

## Installing without the script

If you install FedCAT some other way:

- Put [`variant_density_blocks.tsv`](https://github.com/ExeterGenetics/biobank_workflows/blob/main/regenie_workflows/federated_analysis/setup_files/variant_density_blocks.tsv) (or your own file of roughly 100k-variant blocks, passed with `--variant-density-blockfile`) in the working directory.
- Add the WDL workflows to your workspace (with `fedcat setup-workflows`, or see [`setup_helpers`](https://github.com/ExeterGenetics/biobank_workflows/tree/main/setup_helpers)). To use a Cromwell server on the instance instead, see [Running Cromwell yourself](../commands/cromwell.md).

Next: [Quick start](quick-start.md).
