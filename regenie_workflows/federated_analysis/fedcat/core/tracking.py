"""Task progress and costs of a run's jobs, for progress lines and `--report`.

A [`Job`][fedcat.core.job.Job] is one platform submission: a Cromwell workflow of many
tasks, or a DNAnexus applet job of one. Tasks are counted by status with
[`Platform.task_counts`][fedcat.platforms.base.Platform.task_counts].

Costs come from [`Platform.job_costs`][fedcat.platforms.base.Platform.job_costs] where
the platform reports them. Otherwise they are estimated from the logs task scripts
upload on success: first to last timestamp, at the platform's hourly rate. Failed and
preempted tasks upload no log and set-up time is not counted, so the estimate is a
lower bound.
"""

from __future__ import annotations

import json
import re
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path

from attrs import define
from loguru import logger

from fedcat.core.job import Job
from fedcat.core.structures import BioBank, CondType, DiscoveryConfig, JobState, Phase
from fedcat.platforms.base import Cost, Platform, PlatformError, TaskStatus

_LOG_TIMESTAMP_RE = re.compile(r"^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]", re.MULTILINE)
_TASK_LOG_PATTERN = r".*\.workflow\.log$"


# ── run units ─────────────────────────────────────────────────────────────────


@define(frozen=True)
class LogLocation:
    """Where some of a run unit's task logs are."""

    folder: str
    """Cloud folder."""
    pattern: str
    """Regex matched against file names in `folder`."""


@define(frozen=True)
class RunUnit:
    """The jobs of one state file: a discovery phase or a conditioning step."""

    label: str
    """Name in the report, e.g. `step1` or `BMI singles step 2`."""
    state_file: Path
    """State file recording the jobs."""
    logs: dict[BioBank, list[LogLocation]]
    """Where each biobank's tasks for this unit upload their logs."""


def run_units(config: DiscoveryConfig, run_dir: Path, platforms: dict[BioBank, Platform]) -> list[RunUnit]:
    """Return every discovery phase and conditioning step of a run that has a state file.

    Args:
        config: Run config.
        run_dir: Discovery run directory.
        platforms: Platform per biobank.

    Returns:
        The run units, phases first.
    """
    from fedcat.core.conditioning import ConditioningStep

    units: list[RunUnit] = []
    for phase in Phase:
        state_file = run_dir / f"discovery_state_{phase.value}.json"
        if state_file.exists():
            units.append(RunUnit(phase.value, state_file, _discovery_logs(config, phase, platforms)))

    if not config.cond_suffix:
        return units
    step_re = re.compile(r"_cond_step(\d+)$")
    for pheno in config.pheno_names:
        for cond_type in CondType:
            step = ConditioningStep.from_discovery_config(config, pheno, cond_type)
            step_dirs = [
                (int(m.group(1)), d) for d in step.run_path.glob("*_cond_step*") if (m := step_re.search(d.name))
            ]
            for n, step_dir in sorted(step_dirs):
                for state_file in sorted(step_dir.glob("conditioning_state_*.json")):
                    logs = {
                        bb: [
                            LogLocation(
                                f"{config.output_location(bb)}/meta_analysis/{step.dir_name}/{cond_type.value}_cond_step{n}",
                                _TASK_LOG_PATTERN,
                            )
                        ]
                        for bb in config.cohorts
                    }
                    units.append(RunUnit(f"{pheno} {cond_type.value} step {n}", state_file, logs))
    return units


def _discovery_logs(
    config: DiscoveryConfig, phase: Phase, platforms: dict[BioBank, Platform]
) -> dict[BioBank, list[LogLocation]]:
    """Return where each biobank's tasks of a discovery phase upload their logs."""
    logs: dict[BioBank, list[LogLocation]] = {}
    for bb, cohort in config.cohorts.items():
        meta_analysis = f"{config.output_location(bb)}/meta_analysis"
        match phase:
            case Phase.STEP1:
                logs[bb] = [
                    LogLocation(
                        platforms[bb].results_root(cohort.step1_directory),
                        rf"^{re.escape(config.effective_step1_run_prefix)}_.*_Step1_.*\.workflow\.log$",
                    )
                ]
            case Phase.STEP2_SINGLES:
                logs[bb] = [LogLocation(f"{meta_analysis}/single_variants", _TASK_LOG_PATTERN)]
            case Phase.STEP2_AGG:
                logs[bb] = [LogLocation(f"{meta_analysis}/aggregates", _TASK_LOG_PATTERN)]
            case Phase.META:
                logs[bb] = [
                    LogLocation(
                        f"{meta_analysis}/all_biobank_meta", rf"^{re.escape(config.run_prefix)}_meta_chr.*\.log$"
                    )
                ]
            case Phase.STEP0_SUBSET:
                # One folder per type. Logs the meta job leaves there have a single
                # timestamp, so add no runtime
                logs[bb] = [
                    LogLocation(
                        f"{meta_analysis}/all_biobank_meta/{config.step0_subset_dirname(cond_type)}",
                        rf"^{re.escape(config.run_prefix)}_step0_subset_chr.*\.log$",
                    )
                    for cond_type in CondType
                ]
    return logs


def load_jobs(state_file: Path, conv) -> list[Job]:
    """Return the jobs recorded in a discovery or conditioning state file.

    Args:
        state_file: State file.
        conv: Converter from [`make_conv`][fedcat.core.structures.make_conv].

    Returns:
        The jobs; empty, with a warning, if the file cannot be read.
    """
    try:
        raw_jobs = json.loads(state_file.read_text()).get("jobs", [])
    except (OSError, json.JSONDecodeError) as e:
        logger.warning(f"Skipping unreadable state file {state_file}: {e}")
        return []
    return [conv.structure(raw, Job) for raw in raw_jobs]


async def connect_platforms(
    platforms: dict[BioBank, Platform], biobanks: set[BioBank], *, storage_is_enough: bool = False
) -> dict[BioBank, Platform]:
    """Connect to biobanks' platforms with their job services.

    Args:
        platforms: Platform per biobank.
        biobanks: Biobanks to connect.
        storage_is_enough: Return a biobank whose job service is unreachable connected
            to storage alone, e.g. to report on a run whose Cromwell server has stopped.

    Returns:
        The platforms that connected; failures are logged.
    """
    connected: dict[BioBank, Platform] = {}
    for bb in sorted(biobanks, key=lambda b: b.value):
        try:
            await platforms[bb].connect(need_submit=True)
        except PlatformError as e:
            if not storage_is_enough:
                logger.error(f"Cannot connect to {bb.value}: {e}")
                continue
            logger.warning(f"Reading {bb.value} from storage only: {e}")
            try:
                await platforms[bb].connect(need_submit=False)
            except PlatformError as storage_error:
                logger.error(f"Cannot connect to {bb.value}: {storage_error}")
                continue
        connected[bb] = platforms[bb]
    return connected


# ── status ────────────────────────────────────────────────────────────────────


def format_job_states(jobs: list[Job]) -> str:
    """Return `"x/y jobs completed (a running, b failed)"` for a list of jobs."""
    counts = Counter(j.state for j in jobs)
    return f"{counts[JobState.COMPLETED]}/{len(jobs)} jobs completed ({counts[JobState.RUNNING]} running, {counts[JobState.ERROR]} failed)"


def format_task_counts(counts: Counter[TaskStatus], expected: int = 0) -> str:
    """Return `"x/y tasks completed (a running, b failed, c preempted)"`.

    Args:
        counts: Tasks per status.
        expected: Tasks the jobs were submitted with; the total when larger than the
            tasks counted, as tasks that have not started report no status.

    Returns:
        The progress text, or `"no tasks reported"`.
    """
    total = max(expected, sum(n for s, n in counts.items() if s is not TaskStatus.PREEMPTED))
    if not total:
        return "no tasks reported"
    preempted = f", {counts[TaskStatus.PREEMPTED]} preempted" if counts[TaskStatus.PREEMPTED] else ""
    return (
        f"{counts[TaskStatus.SUCCEEDED]}/{total} tasks completed "
        f"({counts[TaskStatus.RUNNING]} running, {counts[TaskStatus.FAILED]} failed{preempted})"
    )


async def task_counts_by_biobank(
    jobs: list[Job], platforms: dict[BioBank, Platform]
) -> dict[BioBank, Counter[TaskStatus]]:
    """Count the tasks of submitted jobs by status, per biobank.

    Args:
        jobs: Jobs; those without a platform ID are ignored.
        platforms: Platform per biobank.

    Returns:
        Task counts per biobank, omitting biobanks that fail to report.
    """
    ids: dict[BioBank, list[str]] = {}
    for job in jobs:
        if job.id:
            ids.setdefault(job.biobank, []).append(job.id)
    out: dict[BioBank, Counter[TaskStatus]] = {}
    for bb, job_ids in sorted(ids.items(), key=lambda kv: kv[0].value):
        if bb not in platforms:
            continue
        try:
            out[bb] = await platforms[bb].task_counts(job_ids)
        except PlatformError as e:
            logger.debug(f"Could not get {bb.value} task counts: {e}")
    return out


async def progress_line(jobs: list[Job], platforms: dict[BioBank, Platform]) -> str:
    """Return a progress line of job states and the running jobs' task counts.

    Args:
        jobs: Jobs of a phase or conditioning step.
        platforms: Platform per biobank.

    Returns:
        The line, or `""` if there are no jobs.
    """
    if not jobs:
        return ""
    line = f"In progress: {format_job_states(jobs)}"
    running = [j for j in jobs if j.state is JobState.RUNNING]
    for bb, counts in (await task_counts_by_biobank(running, platforms)).items():
        expected = sum(j.expected_tasks for j in running if j.biobank is bb)
        line += f" | {bb.value} {format_task_counts(counts, expected)}"
    return line


# ── cost ──────────────────────────────────────────────────────────────────────


def log_runtime_hours(text: str) -> float | None:
    """Return the hours between the first and last `[YYYY-mm-dd HH:MM:SS]` timestamps of a task log.

    Args:
        text: Log text.

    Returns:
        The hours, or None with fewer than two timestamps.
    """
    stamps = _LOG_TIMESTAMP_RE.findall(text)
    if len(stamps) < 2:
        return None
    first, last = (datetime.strptime(s, "%Y-%m-%d %H:%M:%S") for s in (stamps[0], stamps[-1]))
    return max((last - first).total_seconds(), 0.0) / 3600


async def unit_costs(
    unit: RunUnit, jobs: list[Job], platforms: dict[BioBank, Platform]
) -> dict[BioBank, Cost | None]:
    """Return the cost of a run unit's submitted jobs per biobank.

    Args:
        unit: Run unit.
        jobs: The unit's jobs.
        platforms: Connected platform per biobank.

    Returns:
        Cost per biobank, as reported by the platform or else estimated from task logs;
        None where neither is available.
    """
    out: dict[BioBank, Cost | None] = {}
    for bb in sorted({j.biobank for j in jobs if j.id}, key=lambda b: b.value):
        platform = platforms.get(bb)
        if platform is None:
            out[bb] = None
            continue
        bb_jobs = [j for j in jobs if j.biobank is bb and j.id]
        try:
            cost = await platform.job_costs([j.id for j in bb_jobs])
        except PlatformError as e:
            logger.debug(f"Could not get reported {bb.value} costs for {unit.label}: {e}")
            cost = None
        if cost is None and unit.logs.get(bb):
            cost = await _cost_from_logs(platform, bb_jobs, unit.logs[bb])
        out[bb] = cost
    return out


async def _cost_from_logs(platform: Platform, jobs: list[Job], locations: list[LogLocation]) -> Cost | None:
    """Estimate jobs' cost from their task logs at the platform's hourly rate; None if not possible."""
    rates = {platform.task_hourly_cost(j.target, **j.resource_overrides) for j in jobs}
    if len(rates) != 1 or (rate := rates.pop()) is None:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        # A folder each, as logs in different folders can share a name. A location that
        # cannot be read (e.g. a step-0 type this run did not build) is left out
        read_any = False
        for i, location in enumerate(locations):
            try:
                await platform.download(location.folder, location.pattern, Path(tmp) / str(i))
                read_any = True
            except Exception as e:  # noqa: BLE001 — storage helpers raise backend-specific errors
                logger.warning(f"Could not download task logs from {location.folder}: {e}")
        if not read_any:
            return None
        hours = [
            h
            for p in Path(tmp).glob("*/*")
            if (h := log_runtime_hours(p.read_text(errors="replace"))) is not None
        ]
    return Cost(
        amount=sum(hours) * rate,
        currency="USD",
        basis=f"estimated from {len(hours)} successful task logs at ${rate:.3f}/task-hour",
        estimated=True,
    )
