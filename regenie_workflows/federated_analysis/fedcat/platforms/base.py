"""The [`Platform`][fedcat.platforms.base.Platform] interface, its errors and job types.

A platform is the storage and job service a cohort runs on (GCS with Cromwell, or
DNAnexus), separate from the biobank: several biobanks can share one platform and
differ only in config. Every file and job operation fedcat performs goes through it.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import Counter
from enum import Enum
from pathlib import Path
from typing import Any, ClassVar

from attrs import define, field

from fedcat.core.structures import JobKind, JobResources

# ── Exception hierarchy ───────────────────────────────────────────────────────


class PlatformError(Exception):
    """Base class of platform errors."""


class PathSyntaxError(PlatformError):
    """A reference is not valid for this platform, e.g. a DNAnexus path given to GCS.

    Distinct from a valid reference to nothing, for which `Platform.exists` returns False.
    """


class AuthError(PlatformError):
    """Authentication or authorisation with the platform failed."""


class SubmitError(PlatformError):
    """A job could not be submitted."""


class TransferError(PlatformError):
    """A file upload, download, listing or existence check failed."""


class AbortError(PlatformError):
    """A job could not be aborted (e.g. it had already finished)."""


class TrackingError(PlatformError):
    """Task status or cost for a job could not be retrieved."""


# ── Tracking types ────────────────────────────────────────────────────────────


class TaskStatus(Enum):
    """Platform-neutral status of one task (a Cromwell call attempt or a DNAnexus job)."""

    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    PREEMPTED = "preempted"
    """Counts preempted attempts rather than tasks: the retry is counted separately."""


@define(frozen=True)
class Cost:
    """Cost of some jobs."""

    amount: float
    """Amount."""
    currency: str
    """Currency, e.g. `"USD"`."""
    basis: str
    """How the amount was obtained, for display."""
    estimated: bool = False
    """True when estimated by fedcat rather than reported by the platform."""


# ── Job specification ─────────────────────────────────────────────────────────


@define(frozen=True)
class JobSpec:
    """A job to submit, on any platform; platforms ignore the fields that do not apply to them."""

    description: str
    """Job name."""
    target: str
    """Workflow name (Cromwell) or applet ID (DNAnexus)."""
    inputs: dict[str, Any]
    """Workflow or applet inputs."""
    scratch_dir: Path
    """Local directory for transient files, e.g. the Cromwell inputs JSON."""
    instance_type: str | None = None
    """DNAnexus instance type; None for the applet default."""
    priority: str = "normal"
    """DNAnexus job priority."""
    tags: list[str] = field(factory=list)
    """DNAnexus job tags."""


# ── Platform interface ────────────────────────────────────────────────────────


class Platform(ABC):
    """Storage and job service of a cohort.

    Cloud paths are given without a scheme: `bucket/path` on GCS, `/path` or
    `project-xxx:/path` on DNAnexus.
    """

    kind: ClassVar[str]
    """Identifier matching `platform_kind` in `BIOBANK_PROFILES`: `"gcs_cromwell"` or `"dnanexus"`."""

    @abstractmethod
    async def connect(self, *, need_submit: bool) -> None:
        """Authenticate and set up the clients; repeat calls do nothing more.

        Args:
            need_submit: Also set up what is only needed to submit and poll jobs, such
                as a Cromwell client.

        Raises:
            AuthError: If authentication or an access check fails.
        """

    @abstractmethod
    def canonical_data_ref(self, ref: str) -> str:
        """Return a file reference in the form jobs take file inputs.

        Args:
            ref: File path or ID.

        Returns:
            A `file-xxxx` ID on DNAnexus, a `gs://` path on GCS.
        """

    @abstractmethod
    def canonical_path_ref(self, ref: str) -> str:
        """Return a file reference in the form jobs take path inputs.

        Args:
            ref: File path or ID.

        Returns:
            A `/folder/name` path on DNAnexus, a `gs://` path on GCS.
        """

    def canonical_workflow_ref(self, workflow: str) -> str:
        """Return a configured workflow in the form jobs are submitted to.

        Args:
            workflow: Workflow name or applet ID from the config.

        Returns:
            The applet ID unchanged on DNAnexus; the workflow name on Cromwell.
        """
        return workflow

    # ── resources ─────────────────────────────────────────────────────────
    # DNAnexus picks a machine and priority, Cromwell passes sizes as workflow inputs;
    # each platform overrides what applies to it

    def instance_type_for(self, kind: JobKind, resources: JobResources) -> str | None:
        """Return the machine to run a kind of job on.

        Args:
            kind: Job kind.
            resources: Its configured resources.

        Returns:
            The instance type, or None when the platform does not choose machines per job.
        """
        return None

    def resource_inputs(self, target: str, resources: JobResources) -> dict[str, Any]:
        """Return the workflow inputs setting a job's resources.

        Args:
            target: Workflow name, which prefixes the inputs.
            resources: Configured resources.

        Returns:
            The inputs; empty when the platform sizes jobs another way.
        """
        return {}

    def priority_for(self, kind: JobKind, resources: JobResources) -> str | None:
        """Return the priority to run a kind of job at.

        Args:
            kind: Job kind.
            resources: Its configured resources.

        Returns:
            The priority, or None when the platform has none.
        """
        return None

    @abstractmethod
    async def exists(self, ref: str) -> bool:
        """Return whether a file exists.

        Args:
            ref: File path or ID.

        Raises:
            PathSyntaxError: If `ref` is not valid for this platform.
            TransferError: If the check fails.
        """

    @abstractmethod
    async def list_files(self, folder: str, pattern: str | None = None) -> list[str]:
        """List the file names in a folder.

        Args:
            folder: Cloud folder.
            pattern: Regex the names must match; None for all files.

        Returns:
            File names, without the folder.
        """

    @abstractmethod
    async def download(
        self,
        folder: str,
        pattern: str,
        dest: Path,
        max_workers: int | None = None,
    ) -> None:
        """Download the files of a folder whose names match a regex.

        Args:
            folder: Cloud folder.
            pattern: Regex the file names must match.
            dest: Local directory, created if needed.
            max_workers: Most concurrent downloads; None for the platform default.
        """

    @abstractmethod
    async def upload(self, file: Path, destination: str, *, overwrite: bool = False) -> str:
        """Upload a file into a folder, keeping its name.

        Args:
            file: Local file.
            destination: Cloud folder.
            overwrite: Replace an existing file of the same name.

        Returns:
            A reference usable as a job input: a `gs://` path or a DNAnexus file ID.
        """

    @abstractmethod
    async def check_targets(self, targets: list[str]) -> list[str]:
        """Check workflows or applets exist and can be run.

        Args:
            targets: Workflow names or applet IDs.

        Returns:
            One problem description per target that cannot be run; empty if all can.
        """

    @abstractmethod
    async def submit(self, spec: JobSpec) -> str:
        """Submit a job.

        Args:
            spec: Job to submit.

        Returns:
            The platform job ID.

        Raises:
            SubmitError: If submission fails.
        """

    @abstractmethod
    async def wait(
        self,
        job_id: str,
        *,
        poll_interval: float,
        description: str,
        initial_delay: float | None = None,
    ) -> tuple[str, bool]:
        """Poll a job until it finishes.

        Args:
            job_id: Platform job ID.
            poll_interval: Seconds between polls.
            description: Job name, for logs.
            initial_delay: Seconds before the first poll; None for one `poll_interval`.

        Returns:
            `(final_state, success)`: the platform's final job state, and whether the
            job succeeded.
        """

    @abstractmethod
    async def abort(self, job_id: str) -> None:
        """Ask the platform to stop a job; returns once the request is accepted.

        Args:
            job_id: Platform job ID.

        Raises:
            AbortError: If the platform rejects the request.
        """

    @abstractmethod
    async def task_counts(self, job_ids: list[str]) -> Counter[TaskStatus]:
        """Count jobs' tasks by status.

        A Cromwell workflow has a task per call; a DNAnexus job is one task.

        Args:
            job_ids: Platform job IDs.

        Returns:
            Tasks per status.

        Raises:
            TrackingError: If the statuses cannot be retrieved.
        """

    async def job_costs(self, job_ids: list[str]) -> Cost | None:
        """Return the total cost the platform reports for jobs.

        Args:
            job_ids: Platform job IDs.

        Returns:
            The cost, or None if the platform reports none.
        """
        return None

    def task_hourly_cost(
        self, target: str, *, n_cpu: int | None = None, mem_gb: int | None = None, disk: str | None = None
    ) -> float | None:
        """Estimate the USD per hour of one task of a workflow, to price task logs.

        Args:
            target: Workflow name or applet ID.
            n_cpu: CPUs, overriding the workflow default.
            mem_gb: Memory in GB, overriding the workflow default.
            disk: Disk as `"<GB> SSD|HDD"`, overriding the workflow default.

        Returns:
            USD per task-hour, or None if the platform has no estimate.
        """
        return None

    @abstractmethod
    def results_root(self, root_results_directory: str) -> str:
        """Return the cloud folder of a directory relative to the cohort's storage root.

        Args:
            root_results_directory: Directory relative to the bucket or project root.

        Returns:
            `<bucket>/<dir>` on GCS, `/<dir>` on DNAnexus.
        """
