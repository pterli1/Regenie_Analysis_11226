"""[`DxPlatform`][fedcat.platforms.dnanexus.platform.DxPlatform]: the platform interface on DNAnexus."""

from __future__ import annotations

import asyncio
import re
from collections import Counter
from pathlib import Path
from typing import ClassVar

import dxpy
from attrs import define

from fedcat.core.structures import JobKind, JobResources
from fedcat.platforms.base import (
    AbortError,
    AuthError,
    Cost,
    JobSpec,
    PathSyntaxError,
    Platform,
    SubmitError,
    TaskStatus,
    TrackingError,
    TransferError,
)

from .dx import (
    UKB_APPLET_PROJECT,
    count_files_at_path,
    describe_jobs,
    download_files_async,
    get_file_id,
    get_file_path,
    list_files_async,
    setup_dxpy,
    submit_dx_job,
    terminate_dx_job,
    upload_file_async,
    wait_for_dx_job,
)

_FILE_ID_RE = re.compile(r"^file-[A-Za-z0-9]{24}$")

# Read-only describes run up to 16 at once, outside dx.py's one-at-a-time limit
_EXISTS_SEMAPHORE = asyncio.Semaphore(16)

# Instance type each kind of job runs on when its config names none
DEFAULT_INSTANCE_TYPES: dict[JobKind, str] = {
    JobKind.STEP1: "mem3_ssd1_v2_x16",
    JobKind.STEP2_SINGLES: "mem3_hdd2_v2_x8",
    JobKind.STEP2_AGG: "mem3_hdd2_v2_x8",
    JobKind.META: "mem3_ssd1_v2_x8",
    JobKind.STEP0_SUBSET: "mem3_ssd1_v2_x8",
    JobKind.COND_SINGLES: "mem3_ssd1_v2_x8",
    JobKind.COND_AGG: "mem3_ssd1_v2_x16",
}

# Task status of each DNAnexus job state
_TASK_STATUS: dict[str, TaskStatus] = {
    "idle": TaskStatus.QUEUED,
    "waiting_on_input": TaskStatus.QUEUED,
    "runnable": TaskStatus.QUEUED,
    "running": TaskStatus.RUNNING,
    "waiting_on_output": TaskStatus.RUNNING,
    "debug_hold": TaskStatus.RUNNING,
    "terminating": TaskStatus.RUNNING,
    "done": TaskStatus.SUCCEEDED,
    "failed": TaskStatus.FAILED,
    "terminated": TaskStatus.FAILED,
}


@define(frozen=True, eq=True)
class DxCfg:
    """Settings of a DNAnexus platform; hashable, so platforms with equal settings are shared."""

    project_id: str
    """Project jobs run in and results are written to."""
    applet_project_id: str = UKB_APPLET_PROJECT
    """Project holding the applets."""


class DxPlatform(Platform):
    """DNAnexus storage and applet jobs in one project.

    Files are given as `file-xxxx` IDs or `/folder/name` paths in the project.
    """

    kind: ClassVar[str] = "dnanexus"

    def __init__(self, cfg: DxCfg) -> None:
        """Create the platform; nothing connects until `connect`."""
        self._cfg = cfg

    @property
    def project_id(self) -> str:
        """Project jobs run in."""
        return self._cfg.project_id

    @property
    def applet_project_id(self) -> str:
        """Project holding the applets."""
        return self._cfg.applet_project_id

    # ── connection ────────────────────────────────────────────────────────
    async def connect(self, *, need_submit: bool) -> None:
        """Check the dxpy login, CONTRIBUTE access to the project and VIEW to the applet project.

        `need_submit` makes no difference, as reading files also needs the login.
        """
        requirements: dict[str, str] = {self._cfg.project_id: "CONTRIBUTE"}
        if self._cfg.applet_project_id != self._cfg.project_id:
            requirements[self._cfg.applet_project_id] = "VIEW"
        try:
            await setup_dxpy(requirements)
        except (dxpy.exceptions.DXError, PermissionError, ValueError) as e:
            raise AuthError(str(e)) from e

    # ── reference canonicalisation ────────────────────────────────────────
    def canonical_data_ref(self, ref: str) -> str:
        """Return the file's `file-xxxx` ID."""
        return get_file_id(ref, self._cfg.project_id)

    def canonical_path_ref(self, ref: str) -> str:
        """Return the file's `/folder/name` path."""
        return get_file_path(ref, self._cfg.project_id)

    def instance_type_for(self, kind: JobKind, resources: JobResources) -> str | None:
        """Return the configured instance type, or the kind's default in `DEFAULT_INSTANCE_TYPES`."""
        return resources.instance_type or DEFAULT_INSTANCE_TYPES[kind]

    def priority_for(self, kind: JobKind, resources: JobResources) -> str | None:
        """Return the configured priority, or None to keep the job's own."""
        return resources.priority

    # ── storage ───────────────────────────────────────────────────────────
    async def exists(self, ref: str) -> bool:
        """Return whether a file exists; a path naming several files raises TransferError."""
        if not (_FILE_ID_RE.match(ref) or ref.startswith("/")):
            raise PathSyntaxError(
                f"{ref!r} is not a DNAnexus reference (expected 'file-XXXX' id or '/folder/name' path)"
            )
        async with _EXISTS_SEMAPHORE:
            try:
                if ref.startswith("file-"):
                    await asyncio.to_thread(dxpy.describe, ref, fields={"id": True}, project=self._cfg.project_id)
                    return True
                n_files = await asyncio.to_thread(count_files_at_path, ref, self._cfg.project_id)
                if n_files > 1:
                    raise TransferError(f"several files are named {ref!r}; remove the stale copies")
                return n_files == 1
            except dxpy.exceptions.DXSearchError:
                return False
            except dxpy.exceptions.DXError as e:
                raise TransferError(f"exists({ref!r}) failed: {e}") from e

    async def list_files(self, folder: str, pattern: str | None = None) -> list[str]:
        """List the file names in a project folder and its subfolders."""
        return await list_files_async(project=self._cfg.project_id, folder=folder, pattern=pattern)

    async def download(
        self,
        folder: str,
        pattern: str,
        dest: Path,
        max_workers: int | None = None,
    ) -> None:
        """Download matching files of a project folder and its subfolders."""
        await download_files_async(
            project=self._cfg.project_id,
            folder=folder,
            pattern=pattern,
            destination_folder=dest,
            max_workers=max_workers,
        )

    async def upload(self, file: Path, destination: str, *, overwrite: bool = False) -> str:
        """Upload a file to a project folder, returning its file ID."""
        dxfile = await upload_file_async(file, destination, self._cfg.project_id, overwrite=overwrite)
        return dxfile.get_id()  # type: ignore

    # ── workflow validation ───────────────────────────────────────────────
    async def check_targets(self, targets: list[str]) -> list[str]:
        """Check each applet can be described in the applet project."""
        problems: list[str] = []
        for applet in dict.fromkeys(targets):
            async with _EXISTS_SEMAPHORE:
                try:
                    await asyncio.to_thread(
                        dxpy.describe, applet, fields={"id": True}, project=self._cfg.applet_project_id
                    )
                except (dxpy.exceptions.DXError, dxpy.exceptions.DXSearchError) as e:
                    problems.append(f"{applet}: not found in project {self._cfg.applet_project_id} ({e})")
        return problems

    # ── jobs ──────────────────────────────────────────────────────────────
    async def submit(self, spec: JobSpec) -> str:
        """Run the applet in the project, on `mem3_ssd1_v2_x8` if no instance type is given."""
        try:
            return await submit_dx_job(
                applet_id=spec.target,
                applet_project_id=self._cfg.applet_project_id,
                project_id=self._cfg.project_id,
                input_dict=spec.inputs,
                instance_type=spec.instance_type or "mem3_ssd1_v2_x8",
                priority=spec.priority,
                tags=spec.tags or [spec.description],
                delay_seconds=0,
            )
        except dxpy.exceptions.DXError as e:
            raise SubmitError(str(e)) from e

    async def wait(
        self,
        job_id: str,
        *,
        poll_interval: float,
        description: str,
        initial_delay: float | None = None,
    ) -> tuple[str, bool]:
        """Poll the job until it is done, failed or terminated."""
        return await wait_for_dx_job(
            job_id, poll_interval=poll_interval, description=description, initial_delay=initial_delay
        )

    async def abort(self, job_id: str) -> None:
        """Terminate the job."""
        try:
            await terminate_dx_job(job_id)
        except dxpy.exceptions.DXError as e:
            raise AbortError(str(e)) from e

    async def _describe(self, job_ids: list[str], fields: dict[str, bool]) -> list[dict]:
        """Describe jobs, raising TrackingError on failure."""
        try:
            return await describe_jobs(job_ids, fields)
        except dxpy.exceptions.DXError as e:
            raise TrackingError(str(e)) from e

    async def task_counts(self, job_ids: list[str]) -> Counter[TaskStatus]:
        """Count jobs by state, each job being one task."""
        descriptions = await self._describe(job_ids, {"state": True})
        return Counter(_TASK_STATUS[d["state"]] for d in descriptions if d.get("state") in _TASK_STATUS)

    async def job_costs(self, job_ids: list[str]) -> Cost | None:
        """Return the jobs' total price as DNAnexus reports it.

        None if no price is visible, which needs access to the project's billing.
        """
        descriptions = await self._describe(job_ids, {"totalPrice": True, "currency": True})
        priced = [d for d in descriptions if d.get("totalPrice") is not None]
        if not priced:
            return None
        currency = (priced[0].get("currency") or {}).get("code", "USD")
        return Cost(
            amount=sum(d["totalPrice"] for d in priced),
            currency=currency,
            basis=f"reported by DNAnexus for {len(priced)}/{len(job_ids)} jobs",
        )

    # ── layout ────────────────────────────────────────────────────────────
    def results_root(self, root_results_directory: str) -> str:
        """Return `/<root_results_directory>`."""
        return f"/{root_results_directory.lstrip('/')}"
