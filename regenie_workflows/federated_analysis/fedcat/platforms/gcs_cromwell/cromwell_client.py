"""Clients submitting and polling Cromwell workflows.

- [`CromwellClient`][fedcat.platforms.gcs_cromwell.cromwell_client.CromwellClient]:
  the REST API of a Cromwell server.
- [`WbCliClient`][fedcat.platforms.gcs_cromwell.cromwell_client.WbCliClient]: the
  Verily Workbench `wb workflow job` CLI.
"""

from __future__ import annotations

import asyncio
import json
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import aiohttp

CROMWELL_WORKFLOWS_API = "/api/workflows/v1"

# wb final states, mapped to Cromwell's Succeeded and Failed
_WB_TERMINAL_SUCCESS = frozenset({"COMPLETED"})
_WB_TERMINAL_FAILURE = frozenset({"FAILED", "CANCELLED", "DELETED"})

# Page size of `wb workflow job task list` (its default)
_WB_TASK_PAGE_SIZE = 1000
# `wb workflow job run` takes the inputs as one argument, and Linux caps a single
# argument at 128 KiB (MAX_ARG_STRLEN); stay under it with room for the flag name.
WB_MAX_INPUTS_BYTES = 120_000


# Metadata keys needed to count tasks. Cromwell runs an inner scatter as a sub-workflow
# per outer shard, and only expands sub-workflows when subWorkflowId is asked for
TASK_STATUS_KEYS = ["executionStatus", "subWorkflowId", "failures"]


class CromwellClientError(RuntimeError):
    """A Cromwell or `wb` request failed or returned an unexpected response."""


class CromwellHttpError(CromwellClientError):
    """A Cromwell REST request returned an HTTP error."""

    def __init__(self, message: str, *, status: int | None = None, body: str | None = None):
        """Record the HTTP status and response body with the message."""
        super().__init__(message)
        self.status = status
        self.body = body


@dataclass(frozen=True)
class CromwellStatus:
    """A workflow's status."""

    id: str
    """Workflow ID."""
    status: str
    """Cromwell status, e.g. `Running` or `Succeeded`."""


def _preempted(attempt: dict[str, Any]) -> bool:
    """Return whether Cromwell retried a call attempt because its spot VM was preempted."""
    # Cromwell's wording for a preemption, absent from a maxRetries failure
    return attempt.get("executionStatus") == "RetryableFailure" and "was preempted" in json.dumps(
        attempt.get("failures") or []
    )


def leaf_task_statuses(metadata: dict[str, Any]) -> list[str]:
    """Return the `executionStatus` of the latest attempt of every task in workflow metadata.

    A call running a sub-workflow counts as the tasks inside it, or as itself before
    it has any. Every preempted attempt adds a `Preempted`, and a task whose latest
    attempt was preempted reports only that, so it does not count as failed while its
    retry starts.

    Args:
        metadata: Workflow metadata with sub-workflows expanded.

    Returns:
        One status per task, plus one `Preempted` per preempted attempt.
    """
    statuses: list[str] = []
    for attempts in (metadata.get("calls") or {}).values():
        latest: dict[int, dict[str, Any]] = {}
        for attempt in attempts:
            shard = attempt.get("shardIndex", -1)
            if attempt.get("attempt", 1) >= latest.get(shard, {}).get("attempt", 0):
                latest[shard] = attempt
            if _preempted(attempt):
                statuses.append("Preempted")
        for attempt in latest.values():
            if _preempted(attempt):
                continue
            sub = attempt.get("subWorkflowMetadata")
            inner = leaf_task_statuses(sub) if isinstance(sub, dict) else []
            statuses.extend(inner or [str(attempt.get("executionStatus", ""))])
    return statuses


class CromwellClient:
    """Client of a Cromwell server's REST API (`/api/workflows/v1`)."""

    def __init__(
        self,
        *,
        cromwell_server: str,
        timeout_seconds: int = 60,
        verify_ssl: bool = True,
    ) -> None:
        """Create the client.

        Args:
            cromwell_server: Server URL; the server is called without authentication.
            timeout_seconds: Request timeout.
            verify_ssl: Verify the server's TLS certificate.
        """
        self._server = cromwell_server.rstrip("/")
        self._timeout = aiohttp.ClientTimeout(total=timeout_seconds)
        self._verify_ssl = verify_ssl

    def _workflows_base(self) -> str:
        """Return the URL of the workflows API."""
        return f"{self._server}{CROMWELL_WORKFLOWS_API}"

    async def _request(self, method: str, path: str, *, what: str, **kwargs: Any) -> dict[str, Any]:
        """Make a request to the workflows API and decode the JSON response.

        Args:
            method: HTTP method.
            path: Path under the workflows API.
            what: Name of the operation, for errors.
            **kwargs: Passed to aiohttp, e.g. `params` or `data`.

        Raises:
            CromwellHttpError: On an HTTP error status.
            CromwellClientError: If the response is not JSON.
        """
        async with (
            aiohttp.ClientSession(timeout=self._timeout) as session,
            session.request(
                method,
                f"{self._workflows_base()}{path}",
                ssl=self._verify_ssl,
                **kwargs,
            ) as resp,
        ):
            body = await resp.text()
            if resp.status >= 400:
                raise CromwellHttpError(f"Cromwell {what} failed", status=resp.status, body=body)
        try:
            return json.loads(body)
        except json.JSONDecodeError as e:
            raise CromwellClientError(f"Cromwell {what} returned a non-JSON response: {body[:200]}") from e

    async def check_backends(self) -> dict[str, Any]:
        """Return the server's backends and its default backend."""
        return await self._request("GET", "/backends", what="backends check")

    async def submit_workflow(
        self,
        *,
        wdl_path: Path,
        inputs_json_path: Path,
        options_json: str | None = None,
        labels: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        """Submit a workflow.

        Args:
            wdl_path: WDL file.
            inputs_json_path: Inputs JSON file.
            options_json: Workflow options as JSON.
            labels: Workflow labels.

        Returns:
            The response, with the workflow `id` and `status`.
        """
        form = aiohttp.FormData()
        form.add_field("workflowSource", wdl_path.read_bytes(), filename=wdl_path.name, content_type="text/plain")
        form.add_field(
            "workflowInputs",
            inputs_json_path.read_bytes(),
            filename=inputs_json_path.name,
            content_type="application/json",
        )
        if options_json is not None:
            form.add_field(
                "workflowOptions",
                options_json.encode(),
                filename="options.json",
                content_type="application/json",
            )
        if labels:
            form.add_field("labels", json.dumps(labels), content_type="application/json")

        return await self._request("POST", "", what="submit", data=form)

    async def status(self, workflow_id: str) -> CromwellStatus:
        """Return a workflow's status.

        Raises:
            CromwellClientError: If the response has no status.
        """
        data = await self._request("GET", f"/{workflow_id}/status", what="status")
        workflow_status = data.get("status")
        if not workflow_status:
            raise CromwellClientError(f"Cromwell status response missing 'status': {data}")
        return CromwellStatus(id=str(data.get("id", workflow_id)), status=str(workflow_status))

    async def metadata(
        self,
        workflow_id: str,
        *,
        include_keys: list[str] | None = None,
        exclude_keys: list[str] | None = None,
        expand_subworkflows: bool | None = None,
    ) -> dict[str, Any]:
        """Return a workflow's metadata.

        Asking for only the keys needed keeps the response, and the server's query,
        small. Calls always include `shardIndex` and `attempt`.

        Args:
            workflow_id: Workflow ID.
            include_keys: Keys to include, matched as prefixes.
            exclude_keys: Keys to exclude.
            expand_subworkflows: Include sub-workflows' metadata; None for the server default.
        """
        params = [("includeKey", k) for k in include_keys or []]
        params += [("excludeKey", k) for k in exclude_keys or []]
        if expand_subworkflows is not None:
            params.append(("expandSubWorkflows", "true" if expand_subworkflows else "false"))
        return await self._request("GET", f"/{workflow_id}/metadata", what="metadata", params=params or None)

    async def abort(self, workflow_id: str) -> dict[str, Any]:
        """Abort a workflow."""
        return await self._request("POST", f"/{workflow_id}/abort", what="abort")

    async def task_statuses(self, workflow_id: str) -> list[str]:
        """Return the status of every task of a workflow; see [`leaf_task_statuses`][fedcat.platforms.gcs_cromwell.cromwell_client.leaf_task_statuses]."""
        return leaf_task_statuses(
            await self.metadata(workflow_id, include_keys=TASK_STATUS_KEYS, expand_subworkflows=True)
        )


async def _run_wb_command(*args: str) -> tuple[int, str, str]:
    """Run a `wb` command and return `(returncode, stdout, stderr)`."""
    proc = await asyncio.create_subprocess_exec(
        "wb",
        *args,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout_bytes, stderr_bytes = await proc.communicate()
    return proc.returncode or 0, stdout_bytes.decode(), stderr_bytes.decode()


class WbCliClient:
    """Client submitting and polling workflows registered in the workspace, with the `wb` CLI."""

    def __init__(self, *, output_bucket_id: str) -> None:
        """Create the client.

        Args:
            output_bucket_id: Workbench resource id (not GCS name) of the bucket
                workflows write their logs to.
        """
        self._output_bucket_id = output_bucket_id

    async def submit_workflow(
        self,
        *,
        workflow_name: str,
        inputs_json_path: Path,
        label: str | None = None,
    ) -> str:
        """Run a registered workflow with `wb workflow job run`.

        Args:
            workflow_name: Registered workflow id.
            inputs_json_path: Inputs JSON file, passed as a command-line argument.
            label: Job name; a timestamp is appended to make the `--job-id`.

        Returns:
            The `wb` job ID.

        Raises:
            CromwellClientError: If the inputs exceed `WB_MAX_INPUTS_BYTES`, or `wb` fails.
        """
        # Job id in the Verily style, e.g. plink-conversion_job_2026-07-08T095552534654704Z
        now_ns = time.time_ns()
        date_time_part = datetime.fromtimestamp(now_ns / 1e9, tz=timezone.utc).strftime("%Y-%m-%dT%H%M%S")
        nanos_fraction = now_ns % 1_000_000_000
        timestamp_string = f"{date_time_part}{nanos_fraction:09d}Z"
        job_id_label = f"{label}_{timestamp_string}" if label else f"{workflow_name}_job_{timestamp_string}"
        inputs_json = json.dumps(json.loads(inputs_json_path.read_text()), separators=(",", ":"))
        if len(inputs_json.encode()) > WB_MAX_INPUTS_BYTES:
            raise CromwellClientError(
                f"Inputs for {workflow_name!r} are {len(inputs_json.encode())} bytes, over the "
                f"{WB_MAX_INPUTS_BYTES} a single wb command-line argument allows; lower the cohort's "
                "max_tasks_per_workflow so each workflow carries fewer tasks"
            )
        cmd = [
            "workflow",
            "job",
            "run",
            f"--workflow={workflow_name}",
            f"--output-bucket-id={self._output_bucket_id}",
            f"--output-path=cromwell-logs/{workflow_name}",
            f"--job-id={job_id_label}",
            f"--inputs={inputs_json}",
            # No call caching, as with local Cromwell: tasks act by uploading to GCS, so a
            # cache hit could skip a re-staged script, or outputs deleted to force a re-run
            "--format=JSON",
        ]

        rc, stdout, stderr = await _run_wb_command(*cmd)
        if rc != 0:
            raise CromwellClientError(
                f"wb workflow job run failed for {workflow_name!r} (rc={rc}): {stderr.strip()}"
            )
        try:
            data = json.loads(stdout)
        except json.JSONDecodeError as e:
            raise CromwellClientError(f"wb workflow job run returned non-JSON output: {stdout[:200]}") from e

        job_id = data.get("runId")
        if not job_id:
            raise CromwellClientError(f"wb workflow job run response missing job ID field: {data}")
        return str(job_id)

    async def abort(self, job_id: str) -> None:
        """Cancel a job with `wb workflow job cancel`.

        Raises:
            CromwellClientError: If `wb` fails.
        """
        rc, _, stderr = await _run_wb_command("workflow", "job", "cancel", f"--job-id={job_id}")
        if rc != 0:
            raise CromwellClientError(f"wb workflow job cancel failed for {job_id!r} (rc={rc}): {stderr.strip()}")

    async def task_statuses(self, job_id: str) -> list[str]:
        """Return the status of every task of a job, paging through `wb workflow job task list`.

        The list has no shard index, so a retried shard appears once per attempt.

        Raises:
            CromwellClientError: If `wb` fails.
        """
        statuses: dict[str, str] = {}
        page: str | None = None
        while True:
            args = [
                "workflow",
                "job",
                "task",
                "list",
                f"--job-id={job_id}",
                f"--limit={_WB_TASK_PAGE_SIZE}",
                "--format=JSON",
            ]
            if page:
                args.append(f"--page={page}")
            rc, stdout, stderr = await _run_wb_command(*args)
            if rc != 0:
                raise CromwellClientError(
                    f"wb workflow job task list failed for {job_id!r} (rc={rc}): {stderr.strip()}"
                )
            try:
                tasks = json.loads(stdout)
            except json.JSONDecodeError as e:
                raise CromwellClientError(
                    f"wb workflow job task list returned non-JSON output: {stdout[:200]}"
                ) from e
            # --page names the task to start from, which may be included again
            new = [t for t in tasks if t.get("id") not in statuses]
            for task in new:
                statuses[task["id"]] = str(task.get("status", "")).upper()
            if len(tasks) < _WB_TASK_PAGE_SIZE or not new:
                return list(statuses.values())
            page = tasks[-1]["id"]

    async def status(self, job_id: str) -> CromwellStatus:
        """Return a job's status from `wb workflow job describe`.

        `COMPLETED` becomes `Succeeded` and `FAILED`, `CANCELLED` and `DELETED` become
        `Failed`; other statuses are returned unchanged.

        Raises:
            CromwellClientError: If `wb` fails or returns no status.
        """
        rc, stdout, stderr = await _run_wb_command(
            "workflow",
            "job",
            "describe",
            f"--job-id={job_id}",
            "--format=JSON",
        )
        if rc != 0:
            raise CromwellClientError(
                f"wb workflow job describe failed for {job_id!r} (rc={rc}): {stderr.strip()}"
            )
        try:
            data = json.loads(stdout)
        except json.JSONDecodeError as e:
            raise CromwellClientError(f"wb workflow job describe returned non-JSON output: {stdout[:200]}") from e

        wb_status = str(data.get("status", "")).upper()
        if not wb_status:
            raise CromwellClientError(f"wb workflow job describe response missing 'status': {data}")

        if wb_status in _WB_TERMINAL_SUCCESS:
            mapped = "Succeeded"
        elif wb_status in _WB_TERMINAL_FAILURE:
            mapped = "Failed"
        else:
            mapped = wb_status

        return CromwellStatus(id=job_id, status=mapped)
