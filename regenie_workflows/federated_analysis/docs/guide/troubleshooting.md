# Troubleshooting

**Checking a config without spending anything**
:   `--dry-run` builds a phase's jobs and prints a summary without submitting them; it needs cloud authentication to check the input files. `--list-jobs` prints the current jobs and their states, without authentication.

**Jobs fail for lack of memory**
:   Raise the job's resources under the cohort's `workflows.resources` (`mem_gb` and `n_cpu` for AoU, `instance_type` for UKB), then resume with `--overwrite-config --force-job-retry`.

**Transient Google Batch failures**
:   Tasks are not retried by default, so failures surface rather than being paid for twice. If you expect only transient Batch failures, set `max_retries` (AoU) for that job in `workflows.resources`.

**The config was changed but the run ignores it**
:   The config is saved in the state files; pass `--overwrite-config` to replace the saved values that changed.

**Reading the progress lines**
:   `x/y tasks completed (a running, b failed)`, where `y` is the number of tasks in the running jobs, counting each shard of a nested scatter. `--logging-level DEBUG` logs how long each metadata read took.

**What did it cost?**
:   `--report` prices every phase and conditioning step. A `local_cromwell` run is priced from Cromwell's task metadata, counting every attempt, preempted and failed ones included, at its size and billable time. Once the server is stopped it falls back to successful tasks' logs, which miss those attempts, so report a run before stopping its server. Both are us-central1 list prices, not billed amounts.

**More detail in the logs**
:   Use `--logging-level DEBUG` (or `TRACE`).

**Large phenotype sets (e.g. proteomics)**
:   Set `step1_pheno_group_size` and `step2_pheno_group_size`, split AoU phases with `max_tasks_per_workflow`, and cap the load with `max_concurrent_jobs`. See [Large phenotype sets](configuration.md#large-phenotype-sets).

**Errors from `fedcat_rs` after pulling a new version**
:   Reinstall it (`pip install ./fedcat_rs`); FedCAT may need a newer version.

**Getting help**
:   Contact Harry Wright (<h.wright4@exeter.ac.uk>).
