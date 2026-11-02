# Running and resuming

A run can take days and meet transient cloud failures, so it can be stopped, inspected, fixed and resumed. The options here are all on [`fedcat discover`](../commands/discover.md).

## Resuming

Progress is saved in state files in `<local_run_path>/<run_prefix>_discovery_run/`. To resume, run the same command again:

```bash
fedcat discover --config config.json
```

With the default `--start-phase latest`, it resumes from the latest phase with a state file, or starts at the first of `--phases-to-run`.

++ctrl+c++ pauses the run: launched jobs keep running in the cloud, and resuming picks them up.

To restart from an earlier phase, pass it to `--start-phase`; fedcat asks first, as later phases will be overwritten. `--one-phase` stops after the current phase.

## When jobs fail

A phase stops with an error when any of its jobs fails. The next run retries the failed jobs from their state before the error.

| Option | Effect |
|---|---|
| `--force-job-retry` | Resubmit failed jobs with inputs rebuilt from the current config, leaving running and completed jobs alone. The usual fix after changing the config. |
| `--reset-jobs N` | Reset failed jobs to [state](../commands/discover.md#state-codes) `N` instead of their previous state. |
| `--reset-jobs-list "desc1" "desc2"` | Reset only the jobs with these exact descriptions (see `--list-jobs`). Not available in conditioning phases. |
| `--force-reset` | Reset even if the previous run did not stop with an error, by forcing the phase into `ERROR` first. |

Tasks are not retried by default, so failures surface rather than being paid for twice. If you expect only transient Google Batch failures, set `max_retries` in the job's [resources](configuration.md#resources).

## Changing the config

The config is saved in the state files, so edits to it apply on resume only with `--overwrite-config`, which replaces the saved values that changed. For example, after jobs run out of memory:

```bash
# After raising mem_gb / instance_type under workflows.resources
fedcat discover --config config.json --overwrite-config --force-job-retry
```

## Inspecting a run

`--dry-run`
:   Run the checks, build the jobs and print a summary without submitting anything. Needs cloud authentication to check the input files.

`--list-jobs`
:   List the current jobs and their states, then exit. Needs no authentication.

`--report`
:   Report job states, task counts and costs for every phase and conditioning step, then exit. Costs come from DNAnexus, or from Cromwell's task metadata, which includes preempted and failed attempts. Without either they are estimated from successful tasks' logs, which miss those attempts, so report a `local_cromwell` run before stopping its server. AoU costs are us-central1 list prices, not billed amounts.

`--skip-validation`
:   Skip checking that the workflows and remote input files exist, for a faster start.

## Stopping all jobs

```bash
fedcat discover --config config.json --abort-all
```

aborts every running job in the run's state, across all phases and conditioning steps, after you confirm. Stop any fedcat process driving the run first. State files are unchanged: on resume the aborted jobs count as failed, so rerun them with `--force-job-retry`.

## Starting later

`--delay-start` waits before launching jobs, e.g. to start a large phase overnight. It takes a delay or a time:

```bash
fedcat discover --config config.json --delay-start in 2 hours
fedcat discover --config config.json --delay-start tomorrow 9am
```

The checks run first, so problems show up straight away.

## Conditioning steps

Conditioning runs per phenotype in numbered steps, each with its own state file:

| Option | Effect |
|---|---|
| `--cond-step N` | Start at step `N` (default: detected from the step directories on disk). |
| `--cond-reset-step N` | Redo a failed step from [`StepState`](../commands/discover.md#state-codes) `N`. |
| `--cond-one-step` | Run a single conditioning step, then exit. |

In conditioning phases `--reset-jobs` resets the step's jobs, and `--force-reset` forces the step into `ERROR`.

If a run stopped at `cond_max_steps` with signals left, raise `cond_max_steps` and resume with `--overwrite-config`.

## Logging

Logs go to the terminal and to `logs/` in the run directory. `--logging-level DEBUG` adds detail, such as how long each Cromwell metadata read took; `NOFILE` logs to the terminal only.
