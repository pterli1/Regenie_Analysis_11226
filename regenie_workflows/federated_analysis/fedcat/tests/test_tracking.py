"""Tests for task tracking and cost estimation (``--report`` and progress lines)."""

from __future__ import annotations

import asyncio
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pytest

from fedcat.core.conditioning import ConditioningStep
from fedcat.core.job import Job
from fedcat.core.structures import BioBank, CondType, JobState, WorkflowBackend, make_conv
from fedcat.core.tracking import LogLocation, RunUnit, log_runtime_hours, progress_line, run_units, unit_costs
from fedcat.platforms.base import Cost, TaskStatus
from fedcat.platforms.dnanexus import DxCfg, DxPlatform
from fedcat.platforms.gcs_cromwell import GcsCromwellCfg, GcsCromwellPlatform, cromwell, cromwell_client, pricing
from fedcat.platforms.gcs_cromwell.pricing import WORKFLOW_RESOURCES, DiskType, parse_disk
from fedcat.tests.conftest import FakePlatform
from fedcat.tests.test_pipeline_integration_guards import _make_aou_cfg

_WDL_DIR = Path(__file__).resolve().parents[3] / "aou_regenie_workflows" / "wdls"


def _job(state: JobState, job_id: str, bb: BioBank = BioBank.AOU, target: str = "wf") -> Job:
    """Return a job with the given state and platform ID."""
    return Job(biobank=bb, description=job_id, id=job_id, state=state, target=target)


def _log(start: str, end: str) -> str:
    """Return a task log running from `start` to `end`."""
    return f"[{start}] Starting task\n[{start}] chr: 1\nplain output line\n[{end}] Completed task\n"


# ── log timing ────────────────────────────────────────────────────────────────


def test_log_runtime_hours_spans_first_to_last_timestamp():
    assert log_runtime_hours(_log("2026-09-16 23:30:00", "2026-09-17 01:00:00")) == pytest.approx(1.5)


def test_log_runtime_hours_needs_two_timestamps():
    assert log_runtime_hours("[2026-09-16 23:30:00] Starting task\n") is None
    assert log_runtime_hours("no timestamps") is None


# ── run units ─────────────────────────────────────────────────────────────────


def test_run_units_cover_discovery_phases_and_conditioning_steps(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "discovery_state_step2_singles.json").write_text("{}")
    (run_dir / "discovery_state_meta_analysis.json").write_text("{}")
    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", CondType.AGGREGATE)
    for n in (10, 2):
        state = step.run_path / f"agg_cond_step{n}" / f"conditioning_state_agg_step{n}.json"
        state.parent.mkdir(parents=True)
        state.write_text("{}")

    units = run_units(cfg, run_dir, {BioBank.AOU: FakePlatform()})

    assert [u.label for u in units] == [
        "step2_singles",
        "meta_analysis",
        "TEST_PHENO agg step 2",
        "TEST_PHENO agg step 10",
    ]
    out = cfg.output_location(BioBank.AOU)
    assert units[0].logs[BioBank.AOU][0].folder == f"{out}/meta_analysis/single_variants"
    assert re.match(units[1].logs[BioBank.AOU][0].pattern, "testrun_meta_chr22.log")
    assert units[3].logs[BioBank.AOU][0].folder == f"{out}/meta_analysis/{step.dir_name}/agg_cond_step10"


# ── progress ──────────────────────────────────────────────────────────────────


def test_progress_line_adds_task_counts_of_running_jobs():
    platform = FakePlatform()
    platform.task_statuses = {
        "running": [TaskStatus.RUNNING, TaskStatus.RUNNING, TaskStatus.SUCCEEDED],
        "done": [TaskStatus.SUCCEEDED] * 50,
    }
    jobs = [_job(JobState.RUNNING, "running"), _job(JobState.COMPLETED, "done")]
    line = asyncio.run(progress_line(jobs, {BioBank.AOU: platform}))
    assert line == (
        "In progress: 1/2 jobs completed (1 running, 0 failed) | aou 1/3 tasks completed (2 running, 0 failed)"
    )


def test_progress_line_counts_tasks_that_have_not_started_yet():
    """A scattered job's total comes from its units, not from the tasks Cromwell knows about."""
    platform = FakePlatform()
    platform.task_statuses = {"wf-1": [TaskStatus.SUCCEEDED, TaskStatus.RUNNING]}
    job = _job(JobState.RUNNING, "wf-1")
    job.recipe = {"op": "step1", "units": [[0, "eur"], [0, "afr"], [1, "eur"], [1, "afr"]]}
    line = asyncio.run(progress_line([job], {BioBank.AOU: platform}))
    assert line.endswith("| aou 1/4 tasks completed (1 running, 0 failed)")


def test_cond_job_counts_only_the_tasks_its_strata_run():
    """A chrom that one stratum cannot condition is left out of that stratum's blocks."""
    job = _job(JobState.RUNNING, "wf-1")
    job.recipe = {"op": "cond", "strata": ["eur", "afr"], "input_list": [("1", "a"), ("2", "b")], "n_tasks": 3}
    assert job.expected_tasks == 3
    del job.recipe["n_tasks"]
    assert job.expected_tasks == 4


# ── cost ──────────────────────────────────────────────────────────────────────


def _unit(folder: str = "bucket/logs") -> RunUnit:
    """Return a run unit whose AoU task logs are in `folder`."""
    return RunUnit(
        "step2_singles", Path("state.json"), {BioBank.AOU: [LogLocation(folder, r".*\.workflow\.log$")]}
    )


def test_unit_costs_prefers_platform_reported_cost():
    platform = FakePlatform()
    platform.reported_cost = Cost(12.5, "GBP", "reported")
    platform.hourly_cost = 1.0
    costs = asyncio.run(unit_costs(_unit(), [_job(JobState.COMPLETED, "j1")], {BioBank.AOU: platform}))
    assert costs == {BioBank.AOU: Cost(12.5, "GBP", "reported")}


def test_unit_costs_estimates_from_task_logs():
    platform = FakePlatform()
    platform.hourly_cost = 0.5
    platform.files = {
        "bucket/logs": {
            "a.workflow.log": _log("2026-09-17 10:00:00", "2026-09-17 12:00:00"),
            "b.workflow.log": _log("2026-09-17 10:00:00", "2026-09-17 11:00:00"),
            "b.output_filelist": "not a log",
        }
    }
    cost = asyncio.run(unit_costs(_unit(), [_job(JobState.COMPLETED, "j1")], {BioBank.AOU: platform}))[
        BioBank.AOU
    ]
    assert cost.amount == pytest.approx(1.5)
    assert cost.currency == "USD" and cost.estimated
    assert "2 successful task logs" in cost.basis


def test_unit_costs_unavailable_without_rate_or_submitted_jobs():
    platform = FakePlatform()  # no reported cost, no hourly rate
    assert asyncio.run(unit_costs(_unit(), [_job(JobState.COMPLETED, "j1")], {BioBank.AOU: platform})) == {
        BioBank.AOU: None
    }
    unsubmitted = Job(biobank=BioBank.AOU, description="x", state=JobState.READY)
    assert asyncio.run(unit_costs(_unit(), [unsubmitted], {BioBank.AOU: platform})) == {}


# ── wb / Cromwell ─────────────────────────────────────────────────────────────


def test_wb_task_statuses_pages_until_short_page(monkeypatch):
    pages = {
        None: [{"id": "t1", "status": "RUNNING"}, {"id": "t2", "status": "COMPLETED"}],
        # --page is treated as inclusive here; the repeated t2 must not be double counted.
        "t2": [{"id": "t2", "status": "COMPLETED"}, {"id": "t3", "status": "PENDING"}],
        "t3": [{"id": "t4", "status": "FAILED"}],
    }
    calls = []

    async def fake_run(*args):
        calls.append(args)
        page = next((a.split("=", 1)[1] for a in args if a.startswith("--page=")), None)
        return 0, json.dumps(pages[page]), ""

    monkeypatch.setattr(cromwell_client, "_run_wb_command", fake_run)
    monkeypatch.setattr(cromwell_client, "_WB_TASK_PAGE_SIZE", 2)
    statuses = asyncio.run(cromwell_client.WbCliClient(output_bucket_id="b").task_statuses("job-1"))
    assert sorted(statuses) == ["COMPLETED", "FAILED", "PENDING", "RUNNING"]
    assert len(calls) == 3


def test_cromwell_metadata_task_statuses_use_latest_attempt(monkeypatch):
    client = cromwell_client.CromwellClient(cromwell_server="http://x")

    async def fake_metadata(workflow_id, **kwargs):
        return {
            "calls": {
                "wf.task": [
                    {"shardIndex": 0, "attempt": 1, "executionStatus": "RetryableFailure"},
                    {"shardIndex": 0, "attempt": 2, "executionStatus": "Running"},
                    {"shardIndex": 1, "attempt": 1, "executionStatus": "Done"},
                ]
            }
        }

    monkeypatch.setattr(client, "metadata", fake_metadata)
    assert sorted(asyncio.run(client.task_statuses("wf-1"))) == ["Done", "Running"]


def test_cromwell_task_statuses_report_preemptions_not_failures(monkeypatch):
    client = cromwell_client.CromwellClient(cromwell_server="http://x")
    preempted = [{"message": "Task wf.task:0:1 was preempted for the 1st time. The call will be restarted."}]
    real = [{"message": "Job exited with return code 1"}]

    async def fake_metadata(workflow_id, **kwargs):
        return {
            "calls": {
                "wf.task": [
                    {"shardIndex": 0, "attempt": 1, "executionStatus": "RetryableFailure", "failures": preempted},
                    {"shardIndex": 0, "attempt": 2, "executionStatus": "Running"},
                    # retry not yet in the metadata
                    {"shardIndex": 1, "attempt": 1, "executionStatus": "RetryableFailure", "failures": preempted},
                    {"shardIndex": 2, "attempt": 1, "executionStatus": "RetryableFailure", "failures": real},
                ]
            }
        }

    monkeypatch.setattr(client, "metadata", fake_metadata)
    statuses = asyncio.run(client.task_statuses("wf-1"))
    assert sorted(statuses) == ["Preempted", "Preempted", "RetryableFailure", "Running"]


def test_progress_line_shows_preempted_attempts_outside_the_task_total():
    platform = FakePlatform()
    platform.task_statuses = {"wf-1": [TaskStatus.RUNNING, TaskStatus.SUCCEEDED, TaskStatus.PREEMPTED]}
    line = asyncio.run(progress_line([_job(JobState.RUNNING, "wf-1")], {BioBank.AOU: platform}))
    assert line.endswith("| aou 1/2 tasks completed (1 running, 0 failed, 1 preempted)")


def test_cromwell_task_statuses_count_tasks_inside_sub_workflows(monkeypatch):
    """A nested scatter runs as a sub-workflow per outer shard; its shards are the tasks."""
    client = cromwell_client.CromwellClient(cromwell_server="http://x")
    asked: dict = {}

    def _sub(*statuses: str) -> dict:
        return {
            "calls": {
                "inner.run_regenie_step1": [
                    {"shardIndex": i, "attempt": 1, "executionStatus": s} for i, s in enumerate(statuses)
                ]
            }
        }

    async def fake_metadata(workflow_id, **kwargs):
        asked.update(kwargs)
        return {
            "calls": {
                "wf.ScatterAt_27": [
                    {
                        "shardIndex": 0,
                        "attempt": 1,
                        "executionStatus": "Done",
                        "subWorkflowMetadata": _sub("Done", "Done"),
                    },
                    {
                        "shardIndex": 1,
                        "attempt": 1,
                        "executionStatus": "Running",
                        "subWorkflowMetadata": _sub("Done", "Running"),
                    },
                    # started too recently to have called anything yet
                    {"shardIndex": 2, "attempt": 1, "executionStatus": "Running", "subWorkflowMetadata": {}},
                ]
            }
        }

    monkeypatch.setattr(client, "metadata", fake_metadata)
    statuses = asyncio.run(client.task_statuses("wf-1"))
    assert sorted(statuses) == ["Done", "Done", "Done", "Running", "Running"]
    # Cromwell only expands sub-workflows when subWorkflowId is one of the keys asked for
    assert asked["expand_subworkflows"] is True
    assert "subWorkflowId" in asked["include_keys"]


def test_cromwell_task_counts_map_statuses(monkeypatch):
    class FakeClient:
        async def task_statuses(self, workflow_id):
            return ["PENDING", "RUNNING", "COMPLETED", "COMPLETED", "CANCELLED", "SOMETHING_NEW"]

    monkeypatch.setattr(cromwell, "_CROMWELL_CLIENT", FakeClient())
    counts = asyncio.run(cromwell.cromwell_task_counts(["wf-1"]))
    assert counts == Counter(
        {TaskStatus.QUEUED: 1, TaskStatus.RUNNING: 1, TaskStatus.SUCCEEDED: 2, TaskStatus.FAILED: 1}
    )


# ── cost from Cromwell metadata ───────────────────────────────────────────────


def _attempt(**over) -> dict:
    """One Batch task attempt as Cromwell's metadata describes it: 8 cpu / 32 GB / 1000 HDD."""
    attempt = {
        "attempt": 1,
        "executionStatus": "Done",
        "preemptible": False,
        "runtimeAttributes": {"cpu": "8", "memory": "32 GB", "disks": "local-disk 1000 HDD"},
        "vmStartTime": "2026-09-23T10:00:00.000Z",
        "vmEndTime": "2026-09-23T12:00:00.000Z",
    }
    return {**attempt, **over}


def _on_demand_rate() -> float:
    """Return the hourly cost of a standard 8-CPU on-demand task."""
    return pricing.TaskResources(cpus=8, mem_gb=32, disk_gb=1000, preemptible=False).hourly_cost


def test_metadata_cost_prices_every_attempt_of_a_nested_scatter():
    """Each set's tasks live in a sub-workflow, and a preempted attempt was still paid for."""
    metadata = {
        "calls": {
            "wf.ScatterAt_27": [
                {
                    "shardIndex": 0,
                    "subWorkflowMetadata": {
                        "calls": {
                            "inner.step1": [
                                _attempt(),
                                _attempt(attempt=2, executionStatus="RetryableFailure"),
                            ]
                        }
                    },
                }
            ]
        }
    }
    cost = pricing.cost_from_metadata(metadata)
    assert cost.priced == 2 and cost.unpriced == 0
    assert cost.usd == pytest.approx(2 * 2 * _on_demand_rate())


def test_metadata_cost_uses_spot_rates_and_cromwells_own_rate():
    spot_rate = pricing.TaskResources(cpus=8, mem_gb=32, disk_gb=1000, preemptible=True).hourly_cost
    spot = pricing.cost_from_metadata({"calls": {"wf.t": [_attempt(preemptible=True)]}})
    assert spot.usd == pytest.approx(2 * spot_rate)
    assert spot.usd < 2 * _on_demand_rate()

    # without the per-attempt flag, the runtime attribute says which attempts got a spot VM
    no_flag = _attempt(attempt=2, runtimeAttributes={**_attempt()["runtimeAttributes"], "preemptible": "3"})
    del no_flag["preemptible"]
    assert pricing.cost_from_metadata({"calls": {"wf.t": [no_flag]}}).usd == pytest.approx(2 * spot_rate)

    # vmCostPerHour, when Cromwell's cost catalog service filled it in, wins over the table
    reported = pricing.cost_from_metadata({"calls": {"wf.t": [_attempt(vmCostPerHour=0.5)]}})
    assert reported.usd == pytest.approx(1.0)


def test_metadata_cost_bills_a_running_attempt_up_to_now():
    now = datetime(2026, 9, 23, 13, 0, tzinfo=timezone.utc)
    running = _attempt(executionStatus="Running")
    del running["vmEndTime"]
    assert pricing.cost_from_metadata({"calls": {"wf.t": [running]}}, now).usd == pytest.approx(
        3 * _on_demand_rate()
    )

    # one that ended without recording a time cannot be priced rather than billed to now
    stuck = _attempt(executionStatus="Failed")
    del stuck["vmEndTime"]
    unpriceable = pricing.cost_from_metadata({"calls": {"wf.t": [stuck]}}, now)
    assert (unpriceable.usd, unpriceable.priced, unpriceable.unpriced) == (0.0, 0, 1)


def test_metadata_cost_skips_attempts_it_cannot_size():
    no_sizes = _attempt()
    del no_sizes["runtimeAttributes"]
    cost = pricing.cost_from_metadata({"calls": {"wf.t": [no_sizes, _attempt()]}})
    assert (cost.priced, cost.unpriced) == (1, 1)


def test_job_costs_come_from_metadata_for_a_local_cromwell_cohort(monkeypatch):
    class FakeClient(cromwell_client.CromwellClient):
        def __init__(self):
            super().__init__(cromwell_server="http://x")

        async def metadata(self, workflow_id, **kwargs):
            assert kwargs["expand_subworkflows"] is True
            return {"calls": {"wf.t": [_attempt(), _attempt(attempt=2)]}}

    monkeypatch.setattr(cromwell, "_CROMWELL_CLIENT", FakeClient())
    platform = GcsCromwellPlatform(GcsCromwellCfg(bucket="gs://b", backend=WorkflowBackend.LOCAL_CROMWELL))
    cost = asyncio.run(platform.job_costs(["wf-1"]))
    assert cost.amount == pytest.approx(4 * _on_demand_rate())
    assert cost.currency == "USD" and cost.estimated
    assert "2 task attempts" in cost.basis


def test_job_costs_fall_back_when_the_cromwell_server_cannot_be_reached(monkeypatch):
    class DeadClient(cromwell_client.CromwellClient):
        def __init__(self):
            super().__init__(cromwell_server="http://x")

        async def metadata(self, workflow_id, **kwargs):
            raise cromwell_client.CromwellClientError("connection refused")

    monkeypatch.setattr(cromwell, "_CROMWELL_CLIENT", DeadClient())
    local = GcsCromwellPlatform(GcsCromwellCfg(bucket="gs://b", backend=WorkflowBackend.LOCAL_CROMWELL))
    assert asyncio.run(local.job_costs(["wf-1"])) is None
    # the wb CLI has no metadata to offer either, so that backend never tries
    wb = GcsCromwellPlatform(GcsCromwellCfg(bucket="gs://b", backend=WorkflowBackend.WB_CLI))
    assert asyncio.run(wb.job_costs(["wf-1"])) is None


def test_workflow_resources_match_wdl_runtime_defaults():
    if not _WDL_DIR.is_dir():
        pytest.skip("AoU WDLs not available")
    for stem, resources in WORKFLOW_RESOURCES.items():
        wdl = (_WDL_DIR / f"{stem}.wdl").read_text()
        assert f"Int n_cpu = {resources.cpus}" in wdl, stem
        assert f"Int mem_gb = {resources.mem_gb}" in wdl, stem
        assert f'String disk = "{resources.disk_gb} {resources.disk_type.value}"' in wdl, stem
        # The disk input is what reaches the runtime block, so sizes cannot be set twice.
        assert 'disks: "local-disk " + disk' in wdl, stem
        # every priced workflow takes the retry inputs fedcat may send, and retries
        # nothing unless the config asks for it
        assert "Int max_retries = 0" in wdl, stem
        tries = re.search(r"Int preemptible_tries = (\d+)", wdl)
        assert tries is not None and (int(tries.group(1)) > 0) is resources.preemptible, stem


def test_wdl_floats_have_no_exponent_literals():
    """Cromwell's WDL parser rejects exponent literals such as ``1e-6`` in WDL source."""
    if not _WDL_DIR.is_dir():
        pytest.skip("AoU WDLs not available")
    for wdl in _WDL_DIR.glob("*.wdl"):
        assert not re.search(r"Float \w+ = [-+]?[\d.]+[eE]", wdl.read_text()), wdl.name


# ── DNAnexus ──────────────────────────────────────────────────────────────────


def _patch_describe(monkeypatch, describes: dict[str, dict]):
    """Make dxpy describe jobs from `describes`."""
    import dxpy

    def fake_describe_executions(params, **kwargs):
        return {
            "results": [{"describe": {"id": i, **describes[i]}} for i in params["executions"] if i in describes]
        }

    monkeypatch.setattr(dxpy.api, "system_describe_executions", fake_describe_executions)


def test_dx_task_counts_and_reported_costs(monkeypatch):
    _patch_describe(
        monkeypatch,
        {
            "job-1": {"state": "runnable", "totalPrice": None},
            "job-2": {"state": "running"},
            "job-3": {"state": "done", "totalPrice": 1.25, "currency": {"code": "GBP"}},
            "job-4": {"state": "terminated", "totalPrice": 0.5, "currency": {"code": "GBP"}},
        },
    )
    platform = DxPlatform(DxCfg(project_id="project-000000000000000000000000"))
    ids = ["job-1", "job-2", "job-3", "job-4"]

    assert asyncio.run(platform.task_counts(ids)) == Counter(
        {TaskStatus.QUEUED: 1, TaskStatus.RUNNING: 1, TaskStatus.SUCCEEDED: 1, TaskStatus.FAILED: 1}
    )
    cost = asyncio.run(platform.job_costs(ids))
    assert (cost.amount, cost.currency, cost.estimated) == (1.75, "GBP", False)
    assert "2/4 jobs" in cost.basis


def test_dx_job_costs_none_when_prices_hidden(monkeypatch):
    _patch_describe(monkeypatch, {"job-1": {"state": "done"}})
    platform = DxPlatform(DxCfg(project_id="project-000000000000000000000000"))
    assert asyncio.run(platform.job_costs(["job-1"])) is None


def test_load_jobs_round_trip(tmp_path):
    from fedcat.core.tracking import load_jobs

    state = tmp_path / "discovery_state_step1.json"
    conv = make_conv()
    state.write_text(json.dumps({"jobs": [conv.unstructure(_job(JobState.RUNNING, "id-1"))]}))
    assert [(j.id, j.state) for j in load_jobs(state, conv)] == [("id-1", JobState.RUNNING)]


# ── --report ──────────────────────────────────────────────────────────────────


def test_report_summarises_units(tmp_path, monkeypatch):
    from loguru import logger

    import fedcat.__main__ as cli

    cfg = _make_aou_cfg(tmp_path)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    conv = make_conv()
    jobs = [_job(JobState.RUNNING, "wf-1", target="regenie_step2_single_variants_aou")]
    (run_dir / "discovery_state_step2_singles.json").write_text(
        json.dumps({"jobs": [conv.unstructure(j) for j in jobs]})
    )

    platform = FakePlatform()
    platform.task_statuses = {"wf-1": [TaskStatus.RUNNING, TaskStatus.SUCCEEDED, TaskStatus.SUCCEEDED]}
    platform.hourly_cost = 0.25
    folder = f"{cfg.output_location(BioBank.AOU)}/meta_analysis/single_variants"
    platform.files = {folder: {"a.workflow.log": _log("2026-09-17 10:00:00", "2026-09-17 14:00:00")}}
    monkeypatch.setattr("fedcat.platforms.build_platforms", lambda _cfg: {BioBank.AOU: platform})

    messages: list[str] = []
    sink = logger.add(lambda m: messages.append(m.record["message"]), level="INFO")
    try:
        asyncio.run(cli._report(cfg, run_dir, conv))
    finally:
        logger.remove(sink)

    report = "\n".join(messages)
    assert "step2_singles: 0/1 jobs completed (1 running, 0 failed)" in report
    assert "aou  2/3 tasks completed (1 running, 0 failed); cost: 1.00 USD (estimated from 1 successful" in report
    assert "total cost: 1.00 USD" in report
    assert "price compute and disk at list rates" in report


def test_report_falls_back_to_storage_when_the_job_service_is_unreachable(tmp_path, monkeypatch):
    """A stopped Cromwell server costs the report its task counts, not the whole run."""
    from loguru import logger

    import fedcat.__main__ as cli

    cfg = _make_aou_cfg(tmp_path)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    conv = make_conv()
    jobs = [_job(JobState.COMPLETED, "wf-1", target="regenie_step2_single_variants_aou")]
    (run_dir / "discovery_state_step2_singles.json").write_text(
        json.dumps({"jobs": [conv.unstructure(j) for j in jobs]})
    )

    platform = FakePlatform(submit_connect_error="no Cromwell server on http://localhost:8000")
    platform.hourly_cost = 0.25
    folder = f"{cfg.output_location(BioBank.AOU)}/meta_analysis/single_variants"
    platform.files = {folder: {"a.workflow.log": _log("2026-09-17 10:00:00", "2026-09-17 14:00:00")}}
    monkeypatch.setattr("fedcat.platforms.build_platforms", lambda _cfg: {BioBank.AOU: platform})

    messages: list[str] = []
    sink = logger.add(lambda m: messages.append(m.record["message"]), level="INFO")
    try:
        asyncio.run(cli._report(cfg, run_dir, conv))
    finally:
        logger.remove(sink)

    report = "\n".join(messages)
    assert "aou  unavailable; cost: 1.00 USD (estimated from 1 successful task logs" in report


# ── disk specs ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "spec, expected",
    [
        ("200 SSD", (200, DiskType.SSD)),
        ("1000 HDD", (1000, DiskType.HDD)),
        ("local-disk 500 SSD", (500, DiskType.SSD)),
        ("350 ssd", (350, DiskType.SSD)),
        ("350GB SSD", None),
        ("350 NVME", None),
        ("350", None),
        ("", None),
    ],
)
def test_parse_disk(spec, expected):
    assert parse_disk(spec) == expected


def test_unparseable_disk_falls_back_to_the_workflow_default():
    platform = GcsCromwellPlatform(GcsCromwellCfg(bucket="gs://b"))
    target = "regenie_step2_single_variants_aou"
    assert platform.task_hourly_cost(target, disk="nonsense") == platform.task_hourly_cost(target)


def test_disk_is_priced_by_type_and_size():
    platform = GcsCromwellPlatform(GcsCromwellCfg(bucket="gs://b"))
    target = "regenie_step2_single_variants_aou"  # 1000 HDD by default
    small_ssd = platform.task_hourly_cost(target, disk="200 SSD")
    same_size_hdd = platform.task_hourly_cost(target, disk="200 HDD")
    assert small_ssd > same_size_hdd
    # 200 GB SSD costs the same as 850 GB HDD: 0.17 vs 0.04 per GB-month.
    assert small_ssd == pytest.approx(platform.task_hourly_cost(target, disk="850 HDD"), rel=1e-3)


# ── log line context ──────────────────────────────────────────────────────────


def test_log_label_orders_stratum_stage_pheno_and_current_step():
    import contextvars

    from fedcat.core import logging as fl

    class _Step:
        step_number = 3

    step = _Step()

    def _inside() -> list[str]:
        fl.COVAR_STRATUM_VAR.set("female")
        fl.STAGE_VAR.set("cond_singles")
        fl.PHENO_VAR.set("T2D")
        fl.COND_STEP_VAR.set(step)
        before = fl._log_label()
        step.step_number = 4  # read when logged, so an advancing step is always current
        return [before, fl._log_label()]

    assert contextvars.Context().run(_inside) == [
        "[female cond_singles T2D step3] ",
        "[female cond_singles T2D step4] ",
    ]
    assert contextvars.Context().run(fl._log_label) == ""


def test_log_context_is_scoped_to_its_block():
    import contextvars

    from fedcat.core import logging as fl

    def _inside() -> list[str]:
        fl.STAGE_VAR.set("step0_subset")
        with fl.log_context(fl.PHENO_VAR, "BMI"):
            inside = fl._log_label()
        return [inside, fl._log_label()]

    assert contextvars.Context().run(_inside) == ["[step0_subset BMI] ", "[step0_subset] "]


def test_step0_subset_costs_cover_both_type_folders(tmp_path):
    """Each type logs to its own threshold folder under the same names; both count, markers do not."""
    cfg = _make_aou_cfg(tmp_path)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "discovery_state_step0_subset.json").write_text("{}")
    platform = FakePlatform()
    platform.hourly_cost = 1.0
    (unit,) = run_units(cfg, run_dir, {BioBank.AOU: platform})
    folders = [loc.folder for loc in unit.logs[BioBank.AOU]]
    meta = f"{cfg.output_location(BioBank.AOU)}/meta_analysis/all_biobank_meta"
    assert folders == [f"{meta}/{cfg.step0_subset_dirname(t)}" for t in CondType]
    log = f"{cfg.run_prefix}_step0_subset_chr1.log"
    platform.files = {
        folders[0]: {log: _log("2026-09-17 10:00:00", "2026-09-17 11:00:00")},
        folders[1]: {log: _log("2026-09-17 10:00:00", "2026-09-17 12:00:00")},
    }
    cost = asyncio.run(unit_costs(unit, [_job(JobState.COMPLETED, "j1")], {BioBank.AOU: platform}))[BioBank.AOU]
    assert cost.amount == pytest.approx(3.0)
