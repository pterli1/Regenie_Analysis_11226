"""[`GcsCromwellPlatform`][fedcat.platforms.gcs_cromwell.platform.GcsCromwellPlatform]: the platform interface on GCS and Cromwell."""

from __future__ import annotations

import asyncio
from collections import Counter
from pathlib import Path
from typing import Any, ClassVar

from attrs import define, evolve
from loguru import logger

from fedcat.core.structures import (
    DEFAULT_DEPENDENCIES_VERSION,
    JobResources,
    WorkflowBackend,
    wdl_workflow_name,
)
from fedcat.platforms.base import (
    AbortError,
    Cost,
    JobSpec,
    PathSyntaxError,
    Platform,
    PlatformError,
    SubmitError,
    TaskStatus,
    TrackingError,
    TransferError,
)

from .cromwell import (
    CromwellWorkflowError,
    abort_cromwell_job,
    check_cromwell_server,
    cromwell_metadata_cost,
    cromwell_task_counts,
    setup_cromwell_client,
    submit_cromwell_job,
    wait_for_cromwell_job,
)
from .gcs import (
    blob_exists,
    download_files_async,
    list_files_async,
    upload_file_async,
)
from .pricing import WORKFLOW_RESOURCES, parse_disk
from .workflow_setup import (
    WorkflowSetupError,
    registered_workflows,
    resolve_registered_workflow,
    workspace_bucket_uuid,
)


@define(frozen=True, eq=True)
class GcsCromwellCfg:
    """Settings of a GCS/Cromwell platform; hashable, so platforms with equal settings are shared."""

    bucket: str
    """Workspace bucket name."""
    bucket_resource_id: str = ""
    """Workbench resource id of the bucket."""
    backend: WorkflowBackend = WorkflowBackend.WB_CLI
    """Workflow submission backend."""
    dependencies_version: str = DEFAULT_DEPENDENCIES_VERSION
    """Staged dependencies version; with `wb_cli`, jobs run the workflows registered for its WDLs."""


class GcsCromwellPlatform(Platform):
    """GCS storage with Cromwell workflows, submitted with the Workbench CLI or to a local server.

    Files are given as `bucket/path`, with or without `gs://`.
    """

    kind: ClassVar[str] = "gcs_cromwell"

    def __init__(self, cfg: GcsCromwellCfg) -> None:
        """Create the platform; nothing connects until `connect`."""
        self._cfg = cfg
        # wb_cli: registered workflow id per WDL name, from `wb workflow list`
        self._wb_workflow_ids: dict[str, str] = {}

    # ── connection ────────────────────────────────────────────────────────
    async def connect(self, *, need_submit: bool) -> None:
        """Set up the Cromwell client if `need_submit`, and check a local server runs on Google Batch.

        GCS reads use Application Default Credentials and need no set-up.
        """
        if need_submit:
            await setup_cromwell_client(
                self._cfg.backend,
                output_bucket_id=self._cfg.bucket_resource_id or self._cfg.bucket,
            )
            if self._cfg.backend is WorkflowBackend.LOCAL_CROMWELL:
                try:
                    backends = await check_cromwell_server()
                except CromwellWorkflowError as e:
                    raise PlatformError(str(e)) from e
                if backends.get("defaultBackend") == "Local":
                    raise PlatformError(
                        "Cromwell's default backend is 'Local', so tasks would run on this machine; "
                        "restart it with 'fedcat cromwell start' so it loads the generated cromwell.conf"
                    )

    # ── reference canonicalisation ────────────────────────────────────────
    def canonical_data_ref(self, ref: str) -> str:
        """Return the path without `gs://`, which job inputs add when built."""
        return ref.removeprefix("gs://")

    def canonical_path_ref(self, ref: str) -> str:
        """Return the `gs://` path."""
        return "gs://" + ref.removeprefix("gs://")

    def canonical_workflow_ref(self, workflow: str) -> str:
        """Return the WDL's workflow name, which prefixes its inputs.

        With `wb_cli`, the registered workflow running that WDL is found on submission.
        """
        return wdl_workflow_name(workflow)

    async def _lookup_wb_workflows(self, wdl_names: list[str]) -> dict[str, str | WorkflowSetupError]:
        """Return the registered workflow id of each WDL at this dependencies version, or why there is none.

        Workflows are matched by the WDL's staged path, so may have any name.
        """
        registered, bucket_uuid = await asyncio.gather(
            asyncio.to_thread(registered_workflows),
            asyncio.to_thread(workspace_bucket_uuid, self._cfg.bucket, self._cfg.bucket_resource_id),
        )
        out: dict[str, str | WorkflowSetupError] = {}
        for name in wdl_names:
            try:
                out[name] = resolve_registered_workflow(
                    name, self._cfg.dependencies_version, registered, bucket_uuid
                )
            except WorkflowSetupError as exc:
                out[name] = exc
        return out

    async def _wb_workflow_id(self, wdl_name: str) -> str:
        """Return the registered workflow id running a WDL, looked up once.

        Raises:
            WorkflowSetupError: If no workflow runs the WDL at this version.
        """
        if wdl_name not in self._wb_workflow_ids:
            found = (await self._lookup_wb_workflows([wdl_name]))[wdl_name]
            if isinstance(found, WorkflowSetupError):
                raise found
            self._wb_workflow_ids[wdl_name] = found
        return self._wb_workflow_ids[wdl_name]

    def resource_inputs(self, target: str, resources: JobResources) -> dict[str, Any]:
        """Return the WDL inputs of the resources the config sets; unset ones keep the WDL defaults."""
        prefix = wdl_workflow_name(target)
        return {
            f"{prefix}.{name}": value
            for name, value in (
                ("n_cpu", resources.n_cpu),
                ("mem_gb", resources.mem_gb),
                ("disk", resources.disk),
                ("max_retries", resources.max_retries),
                ("preemptible_tries", resources.preemptible),
            )
            if value is not None
        }

    @staticmethod
    def _as_gs(ref: str) -> str:
        """Return `ref` with a `gs://` prefix."""
        return ref if ref.startswith("gs://") else f"gs://{ref}"

    # ── storage ───────────────────────────────────────────────────────────
    async def exists(self, ref: str) -> bool:
        """Return whether the GCS object exists."""
        full = self._as_gs(ref)
        try:
            return await asyncio.to_thread(blob_exists, full)
        except ValueError as e:
            raise PathSyntaxError(str(e)) from e
        except Exception as e:  # pragma: no cover - defensive
            raise TransferError(f"exists({ref!r}) failed: {e}") from e

    async def list_files(self, folder: str, pattern: str | None = None) -> list[str]:
        """List the file names under a GCS folder, at any depth."""
        return await list_files_async(folder=folder, pattern=pattern)

    async def download(
        self,
        folder: str,
        pattern: str,
        dest: Path,
        max_workers: int | None = None,
    ) -> None:
        """Download matching files under a GCS folder, one at a time; `max_workers` is ignored."""
        await download_files_async(folder=folder, pattern=pattern, destination_folder=dest)

    async def upload(self, file: Path, destination: str, *, overwrite: bool = False) -> str:
        """Upload a file into a GCS folder, returning its `gs://` path; existing files are always replaced."""
        dest = f"{destination.removeprefix('gs://').rstrip('/')}/{file.name}"
        await upload_file_async(file, dest)
        return f"gs://{dest}"

    # ── workflow validation ───────────────────────────────────────────────
    async def check_targets(self, targets: list[str]) -> list[str]:
        """Check each WDL is registered at this version (`wb_cli`) or found locally (`local_cromwell`)."""
        deduped = list(dict.fromkeys(targets))
        if self._cfg.backend is WorkflowBackend.LOCAL_CROMWELL:
            from .workflow_setup import resolve_local_wdl

            return [
                f"{t}: {wdl_workflow_name(t)}.wdl not found in the working directory or the biobank_workflows repo checkout"
                for t in deduped
                if resolve_local_wdl(wdl_workflow_name(t)) is None
            ]
        # WB_CLI: a workflow running the WDL staged for this version must be registered
        try:
            found = await self._lookup_wb_workflows([wdl_workflow_name(t) for t in deduped])
        except WorkflowSetupError as exc:
            return [f"could not list registered Verily workflows: {exc}"]
        problems = []
        for t in deduped:
            match found[wdl_workflow_name(t)]:
                case WorkflowSetupError() as exc:
                    problems.append(f"{t}: {exc}")
                case workflow_id:
                    self._wb_workflow_ids[wdl_workflow_name(t)] = workflow_id
        return problems

    # ── jobs ──────────────────────────────────────────────────────────────
    async def submit(self, spec: JobSpec) -> str:
        """Submit the workflow, writing its inputs JSON to `spec.scratch_dir`."""
        workflow_name = spec.target
        if self._cfg.backend is WorkflowBackend.WB_CLI:
            try:
                workflow_name = await self._wb_workflow_id(spec.target)
            except WorkflowSetupError as e:
                raise SubmitError(f"{spec.target}: {e}") from e
        try:
            return await submit_cromwell_job(
                input_dict=spec.inputs,
                workflow_name=workflow_name,
                input_path=spec.scratch_dir / f"wdl_input_{spec.description}.json",
                description=spec.description,
                delay_seconds=0,
            )
        except CromwellWorkflowError as e:
            raise SubmitError(str(e)) from e

    async def wait(
        self,
        job_id: str,
        *,
        poll_interval: float,
        description: str,
        initial_delay: float | None = None,
    ) -> tuple[str, bool]:
        """Poll the workflow until it finishes."""
        return await wait_for_cromwell_job(
            job_id, poll_interval=poll_interval, description=description, initial_delay=initial_delay
        )

    async def abort(self, job_id: str) -> None:
        """Abort the workflow."""
        try:
            await abort_cromwell_job(job_id)
        except CromwellWorkflowError as e:
            raise AbortError(str(e)) from e

    async def task_counts(self, job_ids: list[str]) -> Counter[TaskStatus]:
        """Count the workflows' call attempts by status."""
        try:
            return await cromwell_task_counts(job_ids)
        except CromwellWorkflowError as e:
            raise TrackingError(str(e)) from e

    async def job_costs(self, job_ids: list[str]) -> Cost | None:
        """Return the workflows' cost priced from Cromwell's task metadata, for `local_cromwell` only.

        None, so the cost is estimated from task logs, with `wb_cli`, which exposes no
        metadata, or when the local server is not running.
        """
        if self._cfg.backend is not WorkflowBackend.LOCAL_CROMWELL:
            return None
        try:
            cost = await cromwell_metadata_cost(job_ids)
        except CromwellWorkflowError as e:
            logger.debug(f"No Cromwell metadata to price this run with, estimating from task logs: {e}")
            return None
        if not cost.priced:
            return None
        unpriced = f", {cost.unpriced} attempts unpriced" if cost.unpriced else ""
        return Cost(
            amount=cost.usd,
            currency="USD",
            basis=(
                f"priced from Cromwell metadata for {cost.priced} task attempts, "
                f"failed and preempted VM time included{unpriced}"
            ),
            estimated=True,
        )

    def task_hourly_cost(
        self, target: str, *, n_cpu: int | None = None, mem_gb: int | None = None, disk: str | None = None
    ) -> float | None:
        """Price one task from the WDL's default resources in `WORKFLOW_RESOURCES`, with any overrides."""
        resources = WORKFLOW_RESOURCES.get(wdl_workflow_name(target))
        if resources is None:
            return None
        disk_gb, disk_type = (disk and parse_disk(disk)) or (resources.disk_gb, resources.disk_type)
        return evolve(
            resources,
            cpus=n_cpu or resources.cpus,
            mem_gb=mem_gb or resources.mem_gb,
            disk_gb=disk_gb,
            disk_type=disk_type,
        ).hourly_cost

    # ── layout ────────────────────────────────────────────────────────────
    def results_root(self, root_results_directory: str) -> str:
        """Return `<bucket>/<root_results_directory>`.

        Raises:
            ValueError: If the bucket is unset.
        """
        if not self._cfg.bucket:
            raise ValueError("bucket is required for a GCS+Cromwell platform")
        return f"{self._cfg.bucket.removeprefix('gs://')}/{root_results_directory.lstrip('/')}"
