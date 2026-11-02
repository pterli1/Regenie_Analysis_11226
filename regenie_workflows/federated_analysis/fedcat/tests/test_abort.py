"""Tests for ``--abort-all``: finding running jobs in persisted state and aborting them."""

from __future__ import annotations

import asyncio
import json

from fedcat.core.abort import abort_jobs, running_jobs
from fedcat.core.conditioning import ConditioningStep
from fedcat.core.job import Job
from fedcat.core.structures import BioBank, CondType, JobState, make_conv
from fedcat.core.tracking import run_units
from fedcat.platforms.gcs_cromwell import cromwell_client
from fedcat.tests.conftest import FakePlatform
from fedcat.tests.test_pipeline_integration_guards import _make_aou_cfg


def _write_state(path, jobs: list[Job]) -> None:
    """Write a state file holding `jobs`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    conv = make_conv()
    path.write_text(json.dumps({"jobs": [conv.unstructure(j) for j in jobs]}))


def _job(desc: str, state: JobState, job_id: str = "", bb: BioBank = BioBank.AOU) -> Job:
    """Return a job with the given state and platform ID."""
    return Job(biobank=bb, description=desc, id=job_id, state=state)


def test_running_jobs_collects_across_discovery_and_conditioning_state(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    run_dir = tmp_path / "run"
    _write_state(
        run_dir / "discovery_state_step2_singles.json",
        [
            _job("done", JobState.COMPLETED, "id-done"),
            _job("running", JobState.RUNNING, "id-run"),
            _job("not_submitted", JobState.READY),
        ],
    )
    # A stale copy of the same job in another phase's state is only aborted once.
    _write_state(run_dir / "discovery_state_meta_analysis.json", [_job("running", JobState.RUNNING, "id-run")])
    cond_dir = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", CondType.SINGLEVAR).run_path
    _write_state(
        cond_dir / "singles_cond_step2" / "conditioning_state_singles_step2.json",
        [_job("cond", JobState.RUNNING, "id-cond")],
    )

    files = [u.state_file for u in run_units(cfg, run_dir, {BioBank.AOU: FakePlatform()})]
    assert len(files) == 3
    assert sorted(j.id for j in running_jobs(files, make_conv())) == ["id-cond", "id-run"]


def test_abort_jobs_reports_accepted_and_tolerates_failures(tmp_path):
    ok, rejecting = FakePlatform(), FakePlatform()
    rejecting.abort_error = "job already finished"
    jobs = [
        _job("a", JobState.RUNNING, "id-a", BioBank.AOU),
        _job("b", JobState.RUNNING, "id-b", BioBank.UKB),
    ]
    n = asyncio.run(abort_jobs(jobs, {BioBank.AOU: ok, BioBank.UKB: rejecting}))
    assert n == 1
    assert ok.aborted == ["id-a"]


def test_abort_jobs_skips_platforms_that_fail_to_connect():
    down = FakePlatform(connect_error="no auth")
    n = asyncio.run(abort_jobs([_job("a", JobState.RUNNING, "id-a")], {BioBank.AOU: down}))
    assert n == 0
    assert down.aborted == []


def test_wb_abort_runs_cancel(monkeypatch):
    calls = []

    async def fake_run(*args):
        calls.append(args)
        return 0, "", ""

    monkeypatch.setattr(cromwell_client, "_run_wb_command", fake_run)
    asyncio.run(cromwell_client.WbCliClient(output_bucket_id="b").abort("run-123"))
    assert calls == [("workflow", "job", "cancel", "--job-id=run-123")]
