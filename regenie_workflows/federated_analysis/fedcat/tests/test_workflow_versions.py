"""Versioned Workbench workflows: registered as <wdl>-v<version>, found by the WDL path they run."""

from __future__ import annotations

import asyncio
import json
import subprocess
import weakref

import pytest
from loguru import logger

from fedcat.core.structures import BioBank, WorkflowBackend, versioned_workflow_name, wdl_workflow_name
from fedcat.platforms.base import JobSpec, SubmitError
from fedcat.platforms.gcs_cromwell import platform as platform_mod
from fedcat.platforms.gcs_cromwell import workflow_setup
from fedcat.platforms.gcs_cromwell.platform import GcsCromwellCfg, GcsCromwellPlatform
from fedcat.platforms.gcs_cromwell.workflow_setup import WorkflowSetupError, resolve_registered_workflow


@pytest.mark.parametrize(
    "workflow, wdl_name, versioned",
    [
        # wb cli doesn't seem to parse dots in workflow ids, so the version's dots become dashes
        ("regenie_step1_scatter_aou", "regenie_step1_scatter_aou", "regenie_step1_scatter_aou-v1-6-2"),
        ("wdls/regenie_step1_scatter_aou.wdl", "regenie_step1_scatter_aou", "regenie_step1_scatter_aou-v1-6-2"),
        # a name already carrying a version is kept, written either way; Path.stem would have cut a ".1"
        ("regenie_step1_scatter_aou-v1-7-1", "regenie_step1_scatter_aou", "regenie_step1_scatter_aou-v1-7-1"),
        ("regenie_step1_scatter_aou-v1.7.1", "regenie_step1_scatter_aou", "regenie_step1_scatter_aou-v1.7.1"),
    ],
)
def test_workflow_names(workflow, wdl_name, versioned):
    assert wdl_workflow_name(workflow) == wdl_name
    assert versioned_workflow_name(workflow, "1.6.2") == versioned


def _wf(wf_id: str, version: str, bucket: str = "uuid-mine", wdl: str = "step1", by: str = "me") -> dict:
    """Return a registered workflow running a WDL staged at a version."""
    path = f"workflows/regenie_wgs_workflow/v{version}/wdls/{wdl}.wdl"
    return {"id": wf_id, "createdBy": by, "bucketSource": {"bucketResourceId": bucket, "path": path}}


REGISTERED = [
    _wf("step1", "1.6.2"),  # registered before versioned names
    _wf("my_own_name", "1.7.0"),  # any name: matched on the WDL it runs
    _wf("step1-v1.8.0", "1.8.0"),
    _wf("step1-v1.8.0-theirs", "1.8.0", bucket="uuid-theirs", by="them"),
    _wf("step1-v1.9.0-theirs", "1.9.0", bucket="uuid-theirs", by="them"),
]


@pytest.mark.parametrize(
    "version, bucket_uuid, expected",
    [("1.6.2", None, "step1"), ("1.7.0", None, "my_own_name"), ("1.8.0", "uuid-mine", "step1-v1.8.0")],
)
def test_registered_workflow_is_found_by_the_wdl_it_runs(version, bucket_uuid, expected):
    assert resolve_registered_workflow("step1", version, REGISTERED, bucket_uuid) == expected


def test_unregistered_or_ambiguous_workflows_are_reported():
    with pytest.raises(WorkflowSetupError, match="setup-workflows --version 2.0.0"):
        resolve_registered_workflow("step1", "2.0.0", REGISTERED)
    # only someone else's bucket has 1.9.0 staged
    with pytest.raises(WorkflowSetupError, match="setup-workflows --version 1.9.0"):
        resolve_registered_workflow("step1", "1.9.0", REGISTERED, "uuid-mine")
    # two users staged 1.8.0 into their own buckets, and this one could not be told apart
    with pytest.raises(WorkflowSetupError, match=r"step1-v1.8.0 \(by me\), step1-v1.8.0-theirs \(by them\)"):
        resolve_registered_workflow("step1", "1.8.0", REGISTERED, None)


def test_workspace_bucket_uuid_from_wb_resource_list(monkeypatch):
    # fields as `wb resource list --type=GCS_BUCKET --format=JSON` gives them
    resources = [
        {"uuid": "uuid-theirs", "id": "Their-Bucket", "bucketName": "their-bucket-wb-x"},
        {"uuid": "uuid-mine", "id": "My-Bucket", "bucketName": "my-bucket-wb-x"},
    ]
    stdout = json.dumps(resources)
    monkeypatch.setattr(workflow_setup, "_run_wb", lambda *_: subprocess.CompletedProcess([], 0, stdout, ""))
    assert workflow_setup.workspace_bucket_uuid("gs://my-bucket-wb-x/") == "uuid-mine"
    assert workflow_setup.workspace_bucket_uuid("other", "My-Bucket") == "uuid-mine"
    assert workflow_setup.workspace_bucket_uuid("unknown") is None


def _platform(monkeypatch, backend=WorkflowBackend.WB_CLI, version="1.7.0") -> GcsCromwellPlatform:
    """Return a GCS platform whose workspace has the REGISTERED workflows."""
    monkeypatch.setattr(platform_mod, "registered_workflows", lambda: REGISTERED)
    monkeypatch.setattr(platform_mod, "workspace_bucket_uuid", lambda *_: "uuid-mine")
    return GcsCromwellPlatform(GcsCromwellCfg(bucket="b", backend=backend, dependencies_version=version))


def test_target_is_the_wdl_name_on_both_backends(monkeypatch):
    for backend in WorkflowBackend:
        assert _platform(monkeypatch, backend).canonical_workflow_ref("wdls/step1-v1.7.0.wdl") == "step1"


def test_wb_check_targets_and_submit_use_the_registered_workflow(monkeypatch, tmp_path):
    platform = _platform(monkeypatch)
    assert asyncio.run(platform.check_targets(["step1"])) == []
    problems = asyncio.run(platform.check_targets(["step2"]))
    assert (
        len(problems) == 1
        and "no workflow registered for workflows/regenie_wgs_workflow/v1.7.0/wdls/step2.wdl" in problems[0]
    )

    submitted = {}

    async def _fake_submit(**kwargs):
        submitted.update(kwargs)
        return "job-1"

    monkeypatch.setattr(platform_mod, "submit_cromwell_job", _fake_submit)
    spec = JobSpec(description="d", target="step1", inputs={}, scratch_dir=tmp_path)
    assert asyncio.run(platform.submit(spec)) == "job-1"
    assert submitted["workflow_name"] == "my_own_name"
    with pytest.raises(SubmitError, match="step2"):
        asyncio.run(platform.submit(JobSpec(description="d", target="step2", inputs={}, scratch_dir=tmp_path)))


def test_job_inputs_and_target_keep_the_wdl_name(tmp_path):
    """WDL inputs are keyed by the workflow name inside the WDL, never by the registered name."""
    import attrs

    from fedcat.core.engine import DiscoveryRun
    from fedcat.core.job import Job
    from fedcat.tests.conftest import make_config

    cfg = make_config(local_run_path=str(tmp_path))
    cfg.aou_bucket = "gs://bkt"
    workflows = cfg.cohorts[BioBank.AOU].workflows
    cfg.cohorts[BioBank.AOU].workflows = attrs.evolve(workflows, regenie_step1="regenie_step1_scatter_aou-v1.6.2")
    disc = DiscoveryRun(config=cfg, run_path=tmp_path, variant_density_blockfile=tmp_path / "b.tsv")
    job = Job(biobank=BioBank.AOU, description="x", recipe={"op": "step1", "units": [[0, "eur"]]})
    job._parent = weakref.ref(disc)
    asyncio.run(job.generate_input())
    assert job.target == "regenie_step1_scatter_aou"
    assert all(k.startswith("regenie_step1_scatter_aou.") for k in job.input_dict)


def test_resource_inputs_and_pricing_use_the_wdl_name(monkeypatch):
    from fedcat.core.structures import JobResources

    platform = _platform(monkeypatch)
    target = platform.canonical_workflow_ref("regenie_step2_aggregates_aou-v1.7.0")
    assert platform.resource_inputs(target, JobResources(n_cpu=4)) == {"regenie_step2_aggregates_aou.n_cpu": 4}
    assert platform.task_hourly_cost(target) is not None


def test_setup_workflows_registers_versioned_names_at_versioned_paths():
    messages: list[str] = []
    sink = logger.add(lambda m: messages.append(m.record["message"]), level="INFO")
    try:
        workflow_setup.setup_workflows(bucket="gs://bkt", bucket_id="bkt-id", version="9.9.9", dry_run=True)
    finally:
        logger.remove(sink)
    assert "- regenie_step1_scatter_aou-v9-9-9" in messages
    assert any(
        "wb workflow create --workflow=regenie_step1_scatter_aou-v9-9-9 --bucket-id=bkt-id "
        "--path=workflows/regenie_wgs_workflow/v9.9.9/wdls/regenie_step1_scatter_aou.wdl" in m
        for m in messages
    )


def test_wb_failures_carry_wbs_own_message(monkeypatch):
    """wb reports why a call failed on its output, which a bare CalledProcessError would hide."""

    def _fail(cmd, **_kwargs):
        raise subprocess.CalledProcessError(
            2, cmd, output="", stderr="[ERROR] Request could not be parsed or was invalid"
        )

    monkeypatch.setattr(workflow_setup.subprocess, "run", _fail)
    with pytest.raises(WorkflowSetupError, match="workflow create --workflow=x.*exit 2.*could not be parsed"):
        workflow_setup._run_wb("workflow", "create", "--workflow=x")
    # callers that tolerate failure still do
    assert workflow_setup.workspace_bucket_uuid("b") is None
    with pytest.raises(WorkflowSetupError, match="Failed to list existing workflows"):
        workflow_setup.registered_workflows()
