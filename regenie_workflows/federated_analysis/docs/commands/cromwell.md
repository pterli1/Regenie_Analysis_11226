# fedcat cromwell

Manage a Cromwell server on the Workbench app for AoU cohorts with `"workflow_backend": "local_cromwell"`, which fedcat then submits workflows to instead of the Verily workflow runner.

```bash
fedcat cromwell setup     # write cromwell.conf via 'wb cromwell generate-config' plus fedcat's settings
fedcat cromwell start     # start the server in the background and wait until it is healthy
fedcat cromwell status    # is it healthy, and on the Google Batch backend?
fedcat cromwell stop      # stop it; running workflows resume when it restarts
```

## How it works

`setup` writes a Google Batch configuration rooted in a workspace bucket, plus fedcat's settings:

- **The Batch helper image.** Cromwell 92's default (`cloud-sdk:461.0.0-alpine`) was deleted by Google's one-year image retention, and its alpine images lack the Python path Batch expects, so tasks fail at startup ([cromwell#7899](https://github.com/broadinstitute/cromwell/issues/7899)).
- **A longer HTTP request timeout** than Cromwell's 40 s, which metadata reads for large scatters can exceed.

Keep `--concurrent-job-limit` within the project's Batch CPU quota.

Files live in `~/cromwell` (`--dir`). The database location and port are passed at startup, not written to `cromwell.conf`. The database is a file, so a restarted server still knows the workflows fedcat is polling.

## Running fedcat against it

`fedcat discover` stops early if no healthy server is running, and reads `CROMWELL_URL` for a non-default port.

It submits the WDLs from the repository with fedcat's default workflow options (no call caching, `ContinueWhilePossible`). A `<workflow>.wdl` or `cromwell_options.json` in the working directory takes priority, to try a change without editing the repository.

`fedcat discover --report` prices the run from Cromwell's task metadata, including failed and preempted attempts, so report a run before stopping its server.

## Options

--8<-- "options/cromwell.md"
