"""`--abort-all`: abort every running cloud job recorded in a run's state files.

Jobs are read from every phase and conditioning step state file, so the abort covers
the whole run. State files are not changed: a resumed run finds the aborted jobs
failed and records them as errors.
"""

from __future__ import annotations

import asyncio
from pathlib import Path

from loguru import logger

from fedcat.core.job import Job
from fedcat.core.structures import BioBank, JobState
from fedcat.core.tracking import connect_platforms, load_jobs
from fedcat.platforms.base import AbortError, Platform


def running_jobs(state_files: list[Path], conv) -> list[Job]:
    """Return the jobs recorded as RUNNING in state files, one per platform job ID.

    Args:
        state_files: Discovery and conditioning state files.
        conv: Converter from [`make_conv`][fedcat.core.structures.make_conv].

    Returns:
        The running jobs.
    """
    jobs: dict[tuple[str, str], Job] = {}
    for path in state_files:
        for job in load_jobs(path, conv):
            if job.state is JobState.RUNNING and job.id:
                jobs.setdefault((job.biobank.value, job.id), job)
    return list(jobs.values())


async def abort_jobs(jobs: list[Job], platforms: dict[BioBank, Platform]) -> int:
    """Abort jobs on their platforms concurrently.

    A job that cannot be aborted, typically because it has finished, is warned about.

    Args:
        jobs: Jobs to abort.
        platforms: Platform per biobank.

    Returns:
        Number of jobs aborted.
    """
    usable = await connect_platforms(platforms, {job.biobank for job in jobs})

    async def _abort(job: Job) -> bool:
        """Abort one job, returning whether the platform accepted it."""
        if job.biobank not in usable:
            return False
        try:
            await usable[job.biobank].abort(job.id)
        except AbortError as e:
            logger.warning(f"Could not abort {job.biobank.value} job {job.description!r} ({job.id}): {e}")
            return False
        logger.info(f"Aborted {job.biobank.value} job {job.description!r} ({job.id})")
        return True

    return sum(await asyncio.gather(*(_abort(job) for job in jobs)))
