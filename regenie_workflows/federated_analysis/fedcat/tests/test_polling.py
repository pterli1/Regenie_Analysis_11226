"""Tests for when polling of a running job starts.

A job submitted in this process cannot be finished yet, so its first poll waits a
full interval (``initial_delay=None``). A job that was already RUNNING when fedcat
picked it up (a resumed run) may have finished while nothing was watching, so it
waits only ``RESUME_POLL_DELAY``.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_polling.py -v
"""

from __future__ import annotations

import asyncio
import weakref

import pytest

from fedcat.core.job import RESUME_POLL_DELAY, Job
from fedcat.core.structures import BioBank, JobState
from fedcat.platforms.dnanexus import dx
from fedcat.platforms.gcs_cromwell import cromwell
from fedcat.tests.conftest import FakePlatform


class _FakeParent:
    """Minimal stand-in for DiscoveryRun/ConditioningStep as a job's parent."""

    def __init__(self, platform: FakePlatform, run_path) -> None:
        self.platforms = {BioBank.AOU: platform}
        self.run_path = run_path

    def persist_state(self, reason: str = "", max_delay: float | None = None) -> None:
        pass


def _run_job(state: JobState, tmp_path) -> FakePlatform:
    """Run a job from `state` on a fake platform, returning the platform."""
    platform = FakePlatform()
    job = Job(
        biobank=BioBank.AOU,
        description="job",
        id="id-1",
        state=state,
        target="wf",
        input_dict={"wf.x": 1},  # already built, so submission does not regenerate input
    )
    job._parent = weakref.ref(parent := _FakeParent(platform, tmp_path))
    assert parent  # keep the parent alive for the weakref
    asyncio.run(job.run())
    return platform


def test_resumed_running_job_is_polled_after_resume_delay(tmp_path):
    assert _run_job(JobState.RUNNING, tmp_path).waited == [("id-1", RESUME_POLL_DELAY)]


def test_job_submitted_in_this_process_waits_a_full_interval(tmp_path):
    platform = _run_job(JobState.READY, tmp_path)
    assert platform.submitted and platform.waited == [("fake-job-1", None)]


@pytest.mark.parametrize("initial_delay, expected_first_sleep", [(0.0, 0.0), (30, 30), (None, 300)])
def test_cromwell_wait_honours_initial_delay(monkeypatch, initial_delay, expected_first_sleep):
    sleeps: list[float] = []

    async def fake_sleep(seconds):
        sleeps.append(seconds)

    class _Client:
        async def status(self, workflow_id):
            return type("S", (), {"status": "Succeeded"})()

    monkeypatch.setattr(cromwell.asyncio, "sleep", fake_sleep)
    monkeypatch.setattr(cromwell, "_CROMWELL_CLIENT", _Client())
    status, ok = asyncio.run(
        cromwell.wait_for_cromwell_job("wf-1", poll_interval=300, initial_delay=initial_delay)
    )
    assert (status, ok) == ("Succeeded", True)
    assert sleeps == [expected_first_sleep]


@pytest.mark.parametrize("initial_delay, expected_first_sleep", [(0.0, 0.0), (30, 30), (None, 300)])
def test_dx_wait_honours_initial_delay(monkeypatch, initial_delay, expected_first_sleep):
    sleeps: list[float] = []

    async def fake_sleep(seconds):
        sleeps.append(seconds)

    monkeypatch.setattr(dx.asyncio, "sleep", fake_sleep)
    monkeypatch.setattr(
        dx.dxpy, "DXJob", lambda job_id: type("J", (), {"describe": lambda self: {"state": "done"}})()
    )
    status, ok = asyncio.run(dx.wait_for_dx_job("job-1", poll_interval=300, initial_delay=initial_delay))
    assert (status, ok) == ("done", True)
    assert sleeps == [expected_first_sleep]
