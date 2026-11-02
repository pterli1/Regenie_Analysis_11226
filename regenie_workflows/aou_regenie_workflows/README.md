# AoU REGENIE workflows

## Contents
- `wdls/` — WDL workflow files run using Cromwell
- `scripts/` — Scripts that are called by the WDL workflows
- `docker/` - Docker with regenie and gcloud that can run some of the wdls (prefer using docker under [../federated_analysis](../federated_analysis/))

## Usage
- The wdl workflows are intended to be added to your workspace and run through FedCAT (see [../federated_analysis](../federated_analysis/)). To run them independently follow steps descibed in [setup_helpers](../../setup_helpers/)
- Inputs are as defined in the WDLs

## Notes
- To make use of these workflows you will need read permissions in the AoU genetic data workspace, and write permissions in the workspace where you are running the workflow and uploading results.
- The scripts and any other dependencies will be downloaded from the `dependencies_path` google bucket folder at runtime, the FedCAT setup-workflows command will upload the dependencies and insert the path as a default for you.

Contact h.wright4@exeter.ac.uk for help
