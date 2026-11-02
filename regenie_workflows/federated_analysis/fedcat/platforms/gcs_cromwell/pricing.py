"""Estimated GCP costs of AoU tasks.

Only compute and disk are priced, at us-central1 list prices for N2 machines (the WDLs
pin an Intel platform) and persistent disk, without discounts. Spot prices change
often; check them at https://cloud.google.com/compute/all-pricing before relying on an
estimate.

[`cost_from_metadata`][fedcat.platforms.gcs_cromwell.pricing.cost_from_metadata] prices
every task attempt in Cromwell's metadata, failed and preempted ones included, at the
size it ran and for the time its VM was billable. Without metadata (`wb_cli`, or a
stopped server), task logs of successful attempts are priced at the
`WORKFLOW_RESOURCES` sizes.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from datetime import datetime, timezone
from enum import Enum
from typing import Any

from attrs import define

# USD per vCPU-hour and GB-hour of memory, N2 on-demand and spot
ON_DEMAND_CPU_HOUR = 0.031611
ON_DEMAND_GB_HOUR = 0.004237
SPOT_CPU_HOUR = 0.00765
SPOT_GB_HOUR = 0.00103


class DiskType(Enum):
    """Persistent disk type, as in a WDL `disks: "local-disk N <type>"`."""

    SSD = "SSD"
    HDD = "HDD"

    @property
    def gb_hour(self) -> float:
        """USD per GB-hour: PD-SSD at $0.17 and PD-Standard at $0.04 per GB-month."""
        return (0.17 if self is DiskType.SSD else 0.04) / 730


def parse_disk(spec: str) -> tuple[int, DiskType] | None:
    """Parse a disk spec such as `"1000 HDD"` or `"local-disk 1000 HDD"`.

    Args:
        spec: Disk spec.

    Returns:
        `(size_gb, disk_type)`, e.g. `(1000, DiskType.HDD)`, or None if `spec` is not of
        that form.
    """
    parts = spec.replace("local-disk", "").split()
    if len(parts) != 2 or not parts[0].isdigit():
        return None
    try:
        return int(parts[0]), DiskType(parts[1].upper())
    except ValueError:
        return None


@define(frozen=True)
class TaskResources:
    """Machine a task runs on."""

    cpus: int
    """vCPUs."""
    mem_gb: float
    """Memory in GB."""
    disk_gb: int
    """Disk in GB."""
    preemptible: bool
    """Whether it runs on a spot VM."""
    disk_type: DiskType = DiskType.HDD
    """Disk type."""

    @property
    def hourly_cost(self) -> float:
        """Estimated USD per hour."""
        cpu, gb = (SPOT_CPU_HOUR, SPOT_GB_HOUR) if self.preemptible else (ON_DEMAND_CPU_HOUR, ON_DEMAND_GB_HOUR)
        return self.cpus * cpu + self.mem_gb * gb + self.disk_gb * self.disk_type.gb_hour


# Default task resources of each WDL (n_cpu, mem_gb, disk, preemptible_tries), which a
# test checks against the WDLs. Configured n_cpu, mem_gb and disk override these when
# pricing; a configured preemptible does not, so a job moved off spot is under-priced
WORKFLOW_RESOURCES: dict[str, TaskResources] = {
    "regenie_step1_scatter_aou": TaskResources(
        cpus=8, mem_gb=32, disk_gb=350, preemptible=True, disk_type=DiskType.SSD
    ),
    "regenie_step2_single_variants_aou": TaskResources(cpus=8, mem_gb=32, disk_gb=1000, preemptible=True),
    "regenie_step2_aggregates_aou": TaskResources(cpus=8, mem_gb=64, disk_gb=1000, preemptible=True),
    "meta_analysis_discovery_aou": TaskResources(
        cpus=16, mem_gb=128, disk_gb=200, preemptible=False, disk_type=DiskType.SSD
    ),
    "step0_subsetting_only_aou": TaskResources(cpus=8, mem_gb=32, disk_gb=400, preemptible=False),
    "regenie_conditioning_singles_aou": TaskResources(cpus=4, mem_gb=32, disk_gb=1000, preemptible=True),
    "regenie_conditioning_agg_aou": TaskResources(cpus=8, mem_gb=64, disk_gb=1000, preemptible=True),
}


# ── costing a run from Cromwell metadata ──────────────────────────────────────

# Metadata keys cost_from_metadata needs; subWorkflowId makes Cromwell expand sub-workflows
COST_METADATA_KEYS = [
    "executionStatus",
    "runtimeAttributes",
    "start",
    "end",
    "vmStartTime",
    "vmEndTime",
    "vmCostPerHour",
    "preemptible",
    "subWorkflowId",
]

# Statuses of attempts still running, billed up to now
_UNFINISHED = frozenset(
    {"NotStarted", "WaitingForQueueSpace", "QueuedInCromwell", "Starting", "Running", "Aborting"}
)


@define(frozen=True)
class MetadataCost:
    """What a set of task attempts cost, and how many could not be priced."""

    usd: float = 0.0
    """Cost of the priced attempts."""
    priced: int = 0
    """Attempts priced."""
    unpriced: int = 0
    """Attempts without usable timing or sizes, e.g. never started."""

    def __add__(self, other: MetadataCost) -> MetadataCost:
        """Sum two costs."""
        return MetadataCost(self.usd + other.usd, self.priced + other.priced, self.unpriced + other.unpriced)


def _parse_time(value: Any) -> datetime | None:
    """Parse a Cromwell ISO-8601 timestamp as an aware datetime, or None."""
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)


def _parse_memory_gb(value: Any) -> float | None:
    """Parse a memory size such as `"32 GB"` to GB, or None; GiB and GB are treated alike."""
    parts = str(value).split()
    try:
        amount = float(parts[0])
    except (IndexError, ValueError):
        return None
    unit = parts[1].upper() if len(parts) > 1 else "GB"
    if unit.startswith(("G", "GI")):
        return amount
    if unit.startswith("M"):
        return amount / 1024
    if unit.startswith("T"):
        return amount * 1024
    return None


def _attempt_preemptible(attempt: Mapping[str, Any], runtime: Mapping[str, Any]) -> bool:
    """Return whether an attempt ran on a spot VM.

    Taken from the attempt's own flag, or else from whether its attempt number is
    within the `preemptible` runtime attribute.
    """
    flag = attempt.get("preemptible")
    if isinstance(flag, bool):
        return flag
    try:
        return int(attempt.get("attempt", 1)) <= int(runtime.get("preemptible", 0))
    except (TypeError, ValueError):
        return False


def attempt_resources(attempt: Mapping[str, Any]) -> TaskResources | None:
    """Return the machine a task attempt ran on, from its `runtimeAttributes`.

    Args:
        attempt: Call attempt metadata.

    Returns:
        The resources, or None if the attributes do not give CPUs, memory and disk.
    """
    runtime = attempt.get("runtimeAttributes")
    if not isinstance(runtime, Mapping):
        return None
    disk = parse_disk(str(runtime.get("disks", "")))
    mem_gb = _parse_memory_gb(runtime.get("memory"))
    try:
        cpus = int(runtime["cpu"])
    except (KeyError, TypeError, ValueError):
        return None
    if disk is None or mem_gb is None:
        return None
    disk_gb, disk_type = disk
    return TaskResources(
        cpus=cpus,
        mem_gb=mem_gb,
        disk_gb=disk_gb,
        preemptible=_attempt_preemptible(attempt, runtime),
        disk_type=disk_type,
    )


def attempt_billable_hours(attempt: Mapping[str, Any], now: datetime) -> float | None:
    """Return the hours a task attempt's VM was billable.

    Uses `vmStartTime` and `vmEndTime`, falling back to the call's `start` and `end`.

    Args:
        attempt: Call attempt metadata.
        now: End time of an attempt still running.

    Returns:
        The hours, or None if they cannot be told.
    """
    start = _parse_time(attempt.get("vmStartTime")) or _parse_time(attempt.get("start"))
    if start is None:
        return None
    end = _parse_time(attempt.get("vmEndTime")) or _parse_time(attempt.get("end"))
    if end is None:
        if attempt.get("executionStatus") not in _UNFINISHED:
            return None
        end = now
    return max((end - start).total_seconds(), 0.0) / 3600


def attempt_cost(attempt: Mapping[str, Any], now: datetime) -> float | None:
    """Return the USD cost of a task attempt.

    Uses Cromwell's `vmCostPerHour` where its cost catalog service provides it, else the
    rate of the machine the attempt ran on.

    Args:
        attempt: Call attempt metadata.
        now: End time of an attempt still running.

    Returns:
        The cost, or None if the attempt cannot be priced.
    """
    hours = attempt_billable_hours(attempt, now)
    if hours is None:
        return None
    rate = attempt.get("vmCostPerHour")
    if isinstance(rate, (int, float)) and not isinstance(rate, bool) and rate >= 0:
        return hours * float(rate)
    resources = attempt_resources(attempt)
    return None if resources is None else hours * resources.hourly_cost


def _attempts(metadata: Mapping[str, Any]) -> Iterator[Mapping[str, Any]]:
    """Yield every task attempt in workflow metadata, in sub-workflows too, including failed ones."""
    for attempts in (metadata.get("calls") or {}).values():
        for attempt in attempts:
            sub = attempt.get("subWorkflowMetadata")
            if isinstance(sub, Mapping):
                yield from _attempts(sub)
            else:
                yield attempt


def cost_from_metadata(metadata: Mapping[str, Any], now: datetime | None = None) -> MetadataCost:
    """Price every task attempt in a workflow's metadata.

    Args:
        metadata: Workflow metadata, with sub-workflows expanded.
        now: End time of attempts still running; defaults to now.

    Returns:
        The total cost, and counts of attempts priced and not.
    """
    now = now or datetime.now(timezone.utc)
    total = MetadataCost()
    for attempt in _attempts(metadata):
        cost = attempt_cost(attempt, now)
        total += MetadataCost(0.0, 0, 1) if cost is None else MetadataCost(cost, 1, 0)
    return total
