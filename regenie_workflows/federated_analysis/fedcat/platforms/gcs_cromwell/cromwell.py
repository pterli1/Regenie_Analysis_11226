"""Cromwell workflow submission, polling, task counts and costs, through one shared client.

The client is set by `setup_cromwell_client` for the cohort's `workflow_backend`:

- `wb_cli`: the Verily Workbench `wb workflow job` CLI.
- `local_cromwell`: a Cromwell server on this machine, without authentication.
"""

from __future__ import annotations

import asyncio
import json
import os
import time
from collections import Counter
from pathlib import Path

from loguru import logger

from fedcat.core.structures import WorkflowBackend, wdl_workflow_name
from fedcat.platforms.base import TaskStatus

from .cromwell_client import CromwellClient, WbCliClient
from .cromwell_server import server_url
from .pricing import COST_METADATA_KEYS, MetadataCost, cost_from_metadata

# ── Module-level state ───────────────────────────────────────────────────────
_CROMWELL_CLIENT: CromwellClient | WbCliClient | None = None
_CROMWELL_CLIENT_CONFIGURED = False
_CROMWELL_API_SEMAPHORE = asyncio.BoundedSemaphore(1)
# Metadata reads grow with the size of a workflow, so queue separately from status
# polls: a slow progress report then delays only the next report
_CROMWELL_METADATA_SEMAPHORE = asyncio.BoundedSemaphore(1)
_CROMWELL_SETUP_LOCK = asyncio.Lock()


def cromwell_url() -> str:
    """Return the local Cromwell server's URL: `CROMWELL_URL`, or the port `fedcat cromwell start` uses."""
    return os.environ.get("CROMWELL_URL", server_url())


class CromwellWorkflowError(Exception):
    """Raised when a Cromwell workflow submission or execution fails."""


# Task status of each `wb workflow job task list` status (upper case) and Cromwell call
# status (camel case); others are not counted. Only a task's latest attempt is read, so
# RetryableFailure counts as failed only when nothing followed it; preempted attempts
# are reported as Preempted (see leaf_task_statuses)
_TASK_STATUS: dict[str, TaskStatus] = {
    "PENDING": TaskStatus.QUEUED,
    "QUEUE_IN_ENGINE": TaskStatus.QUEUED,
    "STARTING": TaskStatus.RUNNING,
    "RUNNING": TaskStatus.RUNNING,
    "STOPPING": TaskStatus.RUNNING,
    "COMPLETED": TaskStatus.SUCCEEDED,
    "SKIPPED": TaskStatus.SUCCEEDED,
    "FAILED": TaskStatus.FAILED,
    "CANCELLED": TaskStatus.FAILED,
    "DELETED": TaskStatus.FAILED,
    "NotStarted": TaskStatus.QUEUED,
    "WaitingForQueueSpace": TaskStatus.QUEUED,
    "QueuedInCromwell": TaskStatus.QUEUED,
    "Starting": TaskStatus.RUNNING,
    "Running": TaskStatus.RUNNING,
    "Aborting": TaskStatus.RUNNING,
    "Done": TaskStatus.SUCCEEDED,
    "Bypassed": TaskStatus.SUCCEEDED,
    "RetryableFailure": TaskStatus.FAILED,
    "Failed": TaskStatus.FAILED,
    "Aborted": TaskStatus.FAILED,
    "Unstartable": TaskStatus.FAILED,
    "Preempted": TaskStatus.PREEMPTED,
}


# ── Client setup ─────────────────────────────────────────────────────────────


async def setup_cromwell_client(
    workflow_backend: WorkflowBackend,
    output_bucket_id: str = "",
) -> None:
    """Create the shared Cromwell client for a backend; later calls do nothing.

    Args:
        workflow_backend: Submission backend. `local_cromwell` connects to
            [`cromwell_url`][fedcat.platforms.gcs_cromwell.cromwell.cromwell_url].
        output_bucket_id: Workbench resource id of the bucket workflows write to;
            required for `wb_cli`.

    Raises:
        ValueError: If `output_bucket_id` is missing for `wb_cli`, or the backend is unknown.
    """
    global _CROMWELL_CLIENT_CONFIGURED, _CROMWELL_CLIENT

    if _CROMWELL_CLIENT_CONFIGURED:
        logger.trace("Cromwell client already configured, skipping setup.")
        return

    async with _CROMWELL_SETUP_LOCK:
        if _CROMWELL_CLIENT_CONFIGURED:
            return

        match workflow_backend:
            case WorkflowBackend.WB_CLI:
                if not output_bucket_id:
                    raise ValueError("output_bucket_id is required for the wb_cli backend.")
                _CROMWELL_CLIENT = WbCliClient(output_bucket_id=output_bucket_id)
                logger.info(f"Cromwell backend: wb CLI (output_bucket_id={output_bucket_id!r}).")

            case WorkflowBackend.LOCAL_CROMWELL:
                url = cromwell_url()
                _CROMWELL_CLIENT = CromwellClient(
                    cromwell_server=url,
                    # The time the generated cromwell.conf gives a metadata request
                    timeout_seconds=300,
                    verify_ssl=False,
                )
                logger.info(f"Cromwell backend: local server at {url}.")

            case _:
                raise ValueError(f"Unsupported workflow backend {workflow_backend!r}")

        _CROMWELL_CLIENT_CONFIGURED = True


async def check_cromwell_server() -> dict:
    """Return the local Cromwell server's `/backends` response; empty with `wb_cli`.

    Raises:
        CromwellWorkflowError: If the server cannot be reached.
    """
    if not isinstance(_CROMWELL_CLIENT, CromwellClient):
        return {}
    try:
        async with _CROMWELL_API_SEMAPHORE:
            return await _CROMWELL_CLIENT.check_backends()
    except Exception as e:
        raise CromwellWorkflowError(
            f"No Cromwell server reachable at {cromwell_url()} ({e}); start one with 'fedcat cromwell start'"
        ) from e


# ── Job submission and polling ───────────────────────────────────────────────


async def submit_cromwell_job(
    input_dict: dict,
    workflow_name: str,
    input_path: Path,
    description: str,
    delay_seconds: int = 0,
) -> str:
    """Write a workflow's inputs JSON and submit it.

    With `local_cromwell`, the WDL and workflow options are sent with the submission,
    taken from the working directory if there, else from the repo checkout and
    fedcat's default options.

    Args:
        input_dict: WDL inputs.
        workflow_name: Registered workflow id (`wb_cli`), or a name whose WDL is
            submitted (`local_cromwell`).
        input_path: Local path to write the inputs JSON to.
        description: Job name, used as the workflow label.
        delay_seconds: Seconds to wait before submitting.

    Returns:
        The workflow ID (Cromwell UUID or `wb` job ID).

    Raises:
        RuntimeError: If the client is not set up.
        CromwellWorkflowError: If the WDL is not found, or the submission is not accepted.
    """
    if delay_seconds:
        await asyncio.sleep(delay_seconds)

    input_path.write_text(json.dumps(input_dict, indent=4))
    logger.debug(f"Submitting job {description!r} ({workflow_name})")

    if _CROMWELL_CLIENT is None:
        raise RuntimeError("Cromwell client not configured. Call setup_cromwell_client() first.")

    if isinstance(_CROMWELL_CLIENT, WbCliClient):
        async with _CROMWELL_API_SEMAPHORE:
            job_id = await _CROMWELL_CLIENT.submit_workflow(
                workflow_name=workflow_name,
                inputs_json_path=input_path,
                label=description,
            )
        return job_id

    # local_cromwell: the WDL and options are sent with every submission
    from .workflow_setup import local_workflow_options, resolve_local_wdl

    stem = wdl_workflow_name(workflow_name)
    wdl_path = resolve_local_wdl(stem)
    if wdl_path is None:
        raise CromwellWorkflowError(
            f"{stem}.wdl not found in the working directory {Path.cwd()} or the biobank_workflows repo checkout"
        )

    async with _CROMWELL_API_SEMAPHORE:
        res = await _CROMWELL_CLIENT.submit_workflow(
            wdl_path=wdl_path,
            inputs_json_path=input_path,
            options_json=local_workflow_options(),
            labels={"description": description} if description else None,
        )

    if not res or "id" not in res:
        raise CromwellWorkflowError(f"Cannot parse Cromwell submit response for {description!r}: {res}")
    status = res.get("status", "")
    if status and status not in ("Submitted", "Running"):
        raise CromwellWorkflowError(f"Unexpected initial Cromwell status '{status}' for {description!r}")
    return str(res["id"])


async def abort_cromwell_job(workflow_id: str) -> None:
    """Abort a running workflow, with `wb workflow job cancel` or Cromwell's `/abort`.

    Args:
        workflow_id: Workflow ID.

    Raises:
        RuntimeError: If the client is not set up.
        CromwellWorkflowError: If the abort is rejected.
    """
    if _CROMWELL_CLIENT is None:
        raise RuntimeError("Cromwell client not configured. Call setup_cromwell_client() first.")
    try:
        async with _CROMWELL_API_SEMAPHORE:
            await _CROMWELL_CLIENT.abort(workflow_id)
    except Exception as e:
        raise CromwellWorkflowError(str(e)) from e


async def cromwell_task_counts(workflow_ids: list[str]) -> Counter[TaskStatus]:
    """Count workflows' tasks by status, from the latest attempt of each.

    Reads each workflow's metadata, which grows with its history, so this is for
    progress reports rather than every poll.

    Args:
        workflow_ids: Workflow IDs.

    Returns:
        Tasks per status.

    Raises:
        CromwellWorkflowError: If there is no client or a read fails.
    """
    if _CROMWELL_CLIENT is None:
        raise CromwellWorkflowError("No Cromwell client to read task statuses from")
    counts: Counter[TaskStatus] = Counter()
    for workflow_id in workflow_ids:
        started = time.monotonic()
        try:
            async with _CROMWELL_METADATA_SEMAPHORE:
                statuses = await _CROMWELL_CLIENT.task_statuses(workflow_id)
        except Exception as e:
            raise CromwellWorkflowError(str(e)) from e
        logger.debug(f"Read {len(statuses)} task statuses for {workflow_id} in {time.monotonic() - started:.1f}s")
        unknown = {s for s in statuses if s not in _TASK_STATUS}
        if unknown:
            logger.debug(f"Ignoring unrecognised task statuses for {workflow_id}: {sorted(unknown)}")
        counts.update(_TASK_STATUS[s] for s in statuses if s in _TASK_STATUS)
    return counts


async def cromwell_metadata_cost(workflow_ids: list[str]) -> MetadataCost:
    """Price every task attempt of workflows from the local Cromwell server's metadata.

    Workflows whose metadata cannot be read are left out with a warning.

    Args:
        workflow_ids: Workflow IDs.

    Returns:
        The summed cost.

    Raises:
        CromwellWorkflowError: If the client is not a local server, or no workflow can be read.
    """
    if not isinstance(_CROMWELL_CLIENT, CromwellClient):
        raise CromwellWorkflowError("Task metadata is only available from a reachable local Cromwell server")
    total = MetadataCost()
    problems: list[str] = []
    for workflow_id in workflow_ids:
        started = time.monotonic()
        try:
            async with _CROMWELL_METADATA_SEMAPHORE:
                metadata = await _CROMWELL_CLIENT.metadata(
                    workflow_id, include_keys=COST_METADATA_KEYS, expand_subworkflows=True
                )
        except Exception as e:  # noqa: BLE001 — any read failure falls back to the estimate
            problems.append(f"{workflow_id}: {e}")
            continue
        cost = cost_from_metadata(metadata)
        logger.debug(
            f"Priced {cost.priced} task attempts of {workflow_id} at ${cost.usd:,.2f} "
            f"in {time.monotonic() - started:.1f}s ({cost.unpriced} unpriced)"
        )
        total += cost
    if problems and len(problems) == len(workflow_ids):
        raise CromwellWorkflowError(f"No workflow metadata could be read ({problems[0]})")
    for problem in problems:
        logger.warning(f"Leaving a workflow out of the cost report: {problem}")
    return total


async def wait_for_cromwell_job(
    workflow_id: str,
    poll_interval: float = 180,
    max_errors: int = 5,
    description: str = "",
    initial_delay: float | None = None,
) -> tuple[str, bool]:
    """Poll a workflow until it is `Succeeded`, `Failed` or `Aborted`.

    `wb` statuses are mapped to these by the client. Errors reading the status are
    retried, as a local server may be restarting or busy while its workflows run.

    Args:
        workflow_id: Workflow ID.
        poll_interval: Seconds between polls.
        max_errors: Consecutive polling errors after which the workflow is treated as failed.
        description: Job name, for logs.
        initial_delay: Seconds before the first poll; None for one `poll_interval`.

    Returns:
        `(status, success)`; the status is `"polling_error"` after too many errors.

    Raises:
        RuntimeError: If the client is not set up.
    """
    logger.debug(f"Polling workflow {description!r} ({workflow_id}), interval {poll_interval}s")
    await asyncio.sleep(poll_interval if initial_delay is None else initial_delay)

    if _CROMWELL_CLIENT is None:
        raise RuntimeError("Cromwell client not configured. Call setup_cromwell_client() first.")

    consecutive_errors = 0
    while True:
        try:
            async with _CROMWELL_API_SEMAPHORE:
                status_obj = await _CROMWELL_CLIENT.status(workflow_id)
            status = status_obj.status
            if status in ("Succeeded", "Failed", "Aborted"):
                return status, status == "Succeeded"
            logger.debug(f"Workflow {description!r} ({workflow_id}) in state '{status}', waiting…")
            consecutive_errors = 0
        except Exception as e:  # noqa: BLE001 — counted towards max_errors before giving up
            consecutive_errors += 1
            logger.warning(
                f"Error {consecutive_errors}/{max_errors} polling workflow {description!r} ({workflow_id}): {e}"
            )
            if consecutive_errors >= max_errors:
                logger.error(
                    f"Giving up tracking {description!r} ({workflow_id}) after {consecutive_errors} polling "
                    "errors; the workflow itself may still be running — check 'fedcat cromwell status'"
                )
                return "polling_error", False
        await asyncio.sleep(poll_interval)
