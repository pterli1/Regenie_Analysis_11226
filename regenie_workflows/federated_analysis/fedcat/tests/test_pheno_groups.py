"""Tests for running many phenotypes: step1 sets, step2 groups, workflow splitting,
skipping completed tasks, the concurrency cap and step2 block-length ceilings.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_pheno_groups.py -v
"""

from __future__ import annotations

import asyncio
import contextlib
import json
import shutil
import subprocess
import weakref
from pathlib import Path

import pytest

from fedcat.core.engine import DiscoveryRun, _step2_block_length
from fedcat.core.job import Job
from fedcat.core.structures import BioBank, JobState, Phase
from fedcat.platforms.gcs_cromwell.cromwell_client import (
    WB_MAX_INPUTS_BYTES,
    CromwellClientError,
    WbCliClient,
)
from fedcat.tests.conftest import FakePlatform
from fedcat.tests.test_autodetect_and_overwrite import _make_cfg

_REPO = Path(__file__).resolve().parents[3]
_WDLS = _REPO / "aou_regenie_workflows" / "wdls"
_SHARED = _REPO / "shared_resources" / "shared_functions.sh"

PHENOS = ["P1", "P2", "P3", "P4", "P5"]
_RUNS: list[DiscoveryRun] = []
# six 100 kb blocks on chr1: a 200 kb ceiling gives three step2 blocks
VAR_BLOCKS = [("1", i * 100_000, (i + 1) * 100_000) for i in range(6)]


def _cfg(tmp_path, **overrides):
    """Return a config with the test phenotypes and any overrides."""
    return _make_cfg(local_run_path=str(tmp_path), pheno_names=PHENOS, **overrides)


def _disc(cfg, tmp_path, phase: Phase = Phase.STEP1, platform: FakePlatform | None = None) -> DiscoveryRun:
    """Return a discovery run at `phase`, on `platform` if given."""
    disc = DiscoveryRun(
        config=cfg,
        run_path=tmp_path / "run",
        variant_density_blockfile=tmp_path / "blocks.tsv",
        current_phase=phase,
    )
    if platform is not None:
        disc.platforms = {bb: platform for bb in cfg.cohorts}
    _RUNS.append(disc)  # jobs hold their run by weakref
    return disc


def _adopt(disc: DiscoveryRun, jobs: list[Job]) -> list[Job]:
    """Make `jobs` the run's jobs."""
    for job in jobs:
        job._parent = weakref.ref(disc)
    disc.jobs = jobs
    return jobs


def _inputs(job: Job) -> dict:
    """Build a job's inputs, returning them without their workflow prefix."""
    asyncio.run(job.generate_input())
    return {k.rsplit(".", 1)[-1]: v for k, v in job.input_dict.items()}


def _step2_jobs(disc: DiscoveryRun, bb: BioBank) -> list[Job]:
    """Return a cohort's step 2 singles jobs, as the run's jobs."""
    return _adopt(disc, disc._build_step2_jobs(bb, disc.config.cohorts[bb], True, VAR_BLOCKS, None))


# ── phenotype grouping ────────────────────────────────────────────────────────


def test_ungrouped_run_keeps_its_names(tmp_path):
    cfg = _cfg(tmp_path)
    assert cfg.step1_pheno_sets == [PHENOS]
    assert cfg.step1_run_prefixes == ["test"]
    assert cfg.step2_pheno_groups == [PHENOS]
    assert cfg.step2_group_label(0) == ""


def test_sets_and_groups_split_in_pheno_order(tmp_path):
    cfg = _cfg(tmp_path, step1_pheno_group_size=2, step2_pheno_group_size=3)
    assert cfg.step1_pheno_sets == [["P1", "P2"], ["P3", "P4"], ["P5"]]
    assert cfg.step1_run_prefixes == ["test_set1", "test_set2", "test_set3"]
    assert cfg.step2_pheno_groups == [["P1", "P2", "P3"], ["P4", "P5"]]
    assert [cfg.step2_group_label(g) for g in range(2)] == ["g1", "g2"]


def test_step2_group_reads_every_step1_set_holding_its_phenotypes(tmp_path):
    """Groups need not align with sets; a group spanning two sets merges both pred lists."""
    cfg = _cfg(tmp_path, step1_pheno_group_size=2, step2_pheno_group_size=3)
    assert cfg.step1_prefixes_for(["P1", "P2", "P3"]) == ["test_set1", "test_set2"]
    assert cfg.step1_prefixes_for(["P4", "P5"]) == ["test_set2", "test_set3"]


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"pheno_names": ["P1", "bad name"]}, "pheno_names may only contain"),
        ({"pheno_names": ["P1", "P1;id"]}, "pheno_names may only contain"),
        ({"pheno_names": ["P1", "P1"]}, "duplicates"),
        ({"step1_pheno_group_size": -1}, "step1_pheno_group_size"),
        ({"max_step2_singles_block_length": 50_000}, "max_step2_singles_block_length"),
    ],
)
def test_validate_rejects(overrides, message, tmp_path):
    with pytest.raises(ValueError, match=message):
        _make_cfg(local_run_path=str(tmp_path), **overrides).validate()


def test_validate_rejects_negative_cohort_limits(tmp_path):
    cfg = _cfg(tmp_path)
    cfg.cohorts[BioBank.AOU].max_concurrent_jobs = -1
    with pytest.raises(ValueError, match="max_concurrent_jobs"):
        cfg.validate()


# ── block length ceilings ─────────────────────────────────────────────────────


def test_singles_ceiling_caps_the_size_derived_length():
    """A small cohort would otherwise get whole-chromosome blocks."""
    assert _step2_block_length("wdl_scatter_v1", True, False, 1000) > 100_000_000
    assert _step2_block_length("wdl_scatter_v1", True, False, 1000, 10_000_000) == 10_000_000


def test_aggregates_ceiling_replaces_the_built_in_one():
    size = 1000  # 1.5e11 / size is far above the built-in 20 Mb ceiling
    assert _step2_block_length("wdl_scatter_v1", False, False, size) == 20_000_000
    assert _step2_block_length("wdl_scatter_v1", False, False, size, 50_000_000) == 50_000_000
    assert _step2_block_length("dx_applet_v1", False, True, size, 50_000_000) == 50_000_000
    assert _step2_block_length("wdl_scatter_v1", False, False, size, 5_000_000) == 5_000_000


def test_block_ceiling_reaches_the_built_blocks(tmp_path):
    cfg = _cfg(tmp_path, max_step2_singles_block_length=200_000)
    jobs = _step2_jobs(_disc(cfg, tmp_path, Phase.STEP2_SINGLES), BioBank.AOU)
    assert [u[2:] for u in jobs[0].units] == [
        ("1", "1", "200000"),
        ("1", "200001", "400000"),
        ("1", "400001", "300000000"),
    ]


# ── step1 ─────────────────────────────────────────────────────────────────────


def test_step1_wdl_splits_over_workflows(tmp_path):
    cfg = _cfg(tmp_path, step1_pheno_group_size=2)
    cfg.cohorts[BioBank.AOU].max_tasks_per_workflow = 2
    disc = _disc(cfg, tmp_path)
    jobs = _adopt(disc, disc._build_step1_jobs(BioBank.AOU, cfg.cohorts[BioBank.AOU]))
    assert [j.description for j in jobs] == [
        "test_aou_regenie_step1_all_strata_part1of2",
        "test_aou_regenie_step1_all_strata_part2of2",
    ]
    first, second = _inputs(jobs[0]), _inputs(jobs[1])
    assert first["pheno_sets"] == ["P1,P2", "P3,P4"]
    assert first["run_prefixes"] == ["test_set1", "test_set2"]
    assert first["strata_by_set"] == [["eur"], ["eur"]]
    assert second["pheno_sets"] == ["P5"] and second["run_prefixes"] == ["test_set3"]


def test_step1_wdl_is_one_workflow_without_a_limit(tmp_path):
    cfg = _cfg(tmp_path, step1_pheno_group_size=2)
    disc = _disc(cfg, tmp_path)
    jobs = disc._build_step1_jobs(BioBank.AOU, cfg.cohorts[BioBank.AOU])
    assert [j.description for j in jobs] == ["test_aou_regenie_step1_all_strata"]


def test_step1_dx_runs_one_job_per_set(tmp_path):
    cfg = _cfg(tmp_path, step1_pheno_group_size=2)
    disc = _disc(cfg, tmp_path)
    jobs = _adopt(disc, disc._build_step1_jobs(BioBank.UKB, cfg.cohorts[BioBank.UKB]))
    assert [j.description for j in jobs] == [f"test_ukb_eur_set{i}_regenie_step1" for i in (1, 2, 3)]
    assert [(_inputs(j)["phenos_to_run"], _inputs(j)["run_prefix"]) for j in jobs] == [
        ("P1,P2", "test_set1"),
        ("P3,P4", "test_set2"),
        ("P5", "test_set3"),
    ]


# ── step2 ─────────────────────────────────────────────────────────────────────


def test_step2_wdl_inputs_carry_each_group(tmp_path):
    cfg = _cfg(tmp_path, step1_pheno_group_size=2, step2_pheno_group_size=3)
    inputs = _inputs(_step2_jobs(_disc(cfg, tmp_path, Phase.STEP2_SINGLES), BioBank.AOU)[0])
    assert inputs["pheno_groups"] == ["P1,P2,P3", "P4,P5"]
    assert inputs["step1_run_prefixes"] == ["test_set1,test_set2", "test_set2,test_set3"]
    assert inputs["group_labels"] == ["g1", "g2"]
    assert inputs["units"][1] == [["eur", "1", "1", "300000000", "None"]]
    for legacy in ("ancestries_array", "singles_blocks", "phenos_to_run", "step1_run_prefix"):
        assert legacy not in inputs


def test_step2_wdl_splits_by_group_first(tmp_path):
    """Units are ordered by group, so a workflow carries few groups' phenotype lists."""
    cfg = _cfg(tmp_path, step2_pheno_group_size=3, max_step2_singles_block_length=200_000)
    cfg.cohorts[BioBank.AOU].max_tasks_per_workflow = 3
    jobs = _step2_jobs(_disc(cfg, tmp_path, Phase.STEP2_SINGLES), BioBank.AOU)
    assert [_inputs(j)["pheno_groups"] for j in jobs] == [["P1,P2,P3"], ["P4,P5"]]


def test_step2_dx_tags_grouped_jobs(tmp_path):
    cfg = _cfg(tmp_path, step1_pheno_group_size=2, step2_pheno_group_size=3)
    jobs = _step2_jobs(_disc(cfg, tmp_path, Phase.STEP2_SINGLES), BioBank.UKB)
    assert [j.description for j in jobs] == [
        "test_ukb_eur_g1_singles_chr1_1_300000000",
        "test_ukb_eur_g2_singles_chr1_1_300000000",
    ]
    inputs = _inputs(jobs[1])
    assert inputs["pheno_group"] == "g2"
    assert inputs["phenos_to_run"] == "P4,P5"
    assert inputs["step1_run_prefix"] == "test_set2,test_set3"


def test_ungrouped_step2_dx_payload_is_unchanged(tmp_path):
    jobs = _step2_jobs(_disc(_cfg(tmp_path), tmp_path, Phase.STEP2_SINGLES), BioBank.UKB)
    inputs = _inputs(jobs[0])
    assert "pheno_group" not in inputs
    assert inputs["step1_run_prefix"] == "test"
    assert jobs[0].description == "test_ukb_eur_singles_chr1_1_300000000"


def test_expected_output_lists_carry_the_group_tag(tmp_path):
    cfg = _cfg(tmp_path, step2_pheno_group_size=3)
    disc = _disc(cfg, tmp_path, Phase.STEP2_SINGLES)
    _step2_jobs(disc, BioBank.AOU)
    assert disc._expected_step2_output_lists(is_singles=True)[BioBank.AOU] == {
        "test_eur_chr1_1_300000000_g1_step2_singlevar.output_filelist",
        "test_eur_chr1_1_300000000_g2_step2_singlevar.output_filelist",
    }


def test_conditioning_reads_the_step1_set_holding_its_phenotype(tmp_path):
    cfg = _cfg(tmp_path, step1_pheno_group_size=2)
    assert cfg.step1_prefixes_for(["P3"]) == ["test_set2"]


def test_changed_grouping_fails_the_job_rather_than_misassigning(tmp_path):
    cfg = _cfg(tmp_path, step2_pheno_group_size=3)
    disc = _disc(cfg, tmp_path, Phase.STEP2_SINGLES)
    job = _step2_jobs(disc, BioBank.UKB)[1]
    cfg.step2_pheno_group_size = 0
    asyncio.run(job.generate_input())
    assert job.state == JobState.ERROR


# ── skipping completed tasks ──────────────────────────────────────────────────


def _no_persist(*_args, **_kwargs) -> None:
    """Do nothing, in place of saving state."""


def test_completed_units_are_left_out_of_a_resubmission(tmp_path, monkeypatch):
    monkeypatch.setattr(DiscoveryRun, "persist_state", _no_persist)
    platform = FakePlatform()
    cfg = _cfg(tmp_path, step2_pheno_group_size=3, max_step2_singles_block_length=200_000)
    disc = _disc(cfg, tmp_path, Phase.STEP2_SINGLES, platform)
    job = _step2_jobs(disc, BioBank.AOU)[0]
    asyncio.run(job.generate_input())
    done = job.units[:4]
    folder = f"{cfg.cohorts[BioBank.AOU].output_location}/meta_analysis/single_variants"
    platform.listing[folder] = [disc._unit_marker(job, u) for u in done] + ["unrelated.regenie.gz"]

    asyncio.run(disc._skip_completed_units())

    assert job.pending_units == job.units[4:]
    assert job.state == JobState.READY
    inputs = _inputs(job)
    assert inputs["group_labels"] == ["g2"]
    assert [tuple(u[:4]) for u in inputs["units"][0]] == [u[1:] for u in job.units[4:]]


def test_a_fully_completed_job_is_not_submitted(tmp_path, monkeypatch):
    monkeypatch.setattr(DiscoveryRun, "persist_state", _no_persist)
    platform = FakePlatform()
    cfg = _cfg(tmp_path, step1_pheno_group_size=2)
    disc = _disc(cfg, tmp_path, Phase.STEP1, platform)
    jobs = _adopt(disc, disc._build_step1_jobs(BioBank.UKB, cfg.cohorts[BioBank.UKB]))
    for job in jobs:
        asyncio.run(job.generate_input())
    platform.listing["/s1"] = [disc._unit_marker(jobs[0], jobs[0].units[0])]

    asyncio.run(disc._skip_completed_units())

    assert [j.state for j in jobs] == [JobState.COMPLETED, JobState.READY, JobState.READY]


def test_listing_failure_submits_everything(tmp_path, monkeypatch):
    monkeypatch.setattr(DiscoveryRun, "persist_state", _no_persist)

    class Unlistable(FakePlatform):
        async def list_files(self, folder, pattern=None):
            raise RuntimeError("denied")

    disc = _disc(_cfg(tmp_path), tmp_path, Phase.STEP1, Unlistable())
    jobs = _adopt(disc, disc._build_step1_jobs(BioBank.UKB, disc.config.cohorts[BioBank.UKB]))
    asyncio.run(jobs[0].generate_input())
    asyncio.run(disc._skip_completed_units())
    assert jobs[0].state == JobState.READY and "pending" not in jobs[0].recipe


# ── concurrency cap ───────────────────────────────────────────────────────────


def test_max_concurrent_jobs_caps_running_jobs(tmp_path, monkeypatch):
    monkeypatch.setattr(DiscoveryRun, "persist_state", _no_persist)
    real_sleep = asyncio.sleep

    async def no_stagger(_delay, *args, **kwargs):
        await real_sleep(0)

    class Tracking(FakePlatform):
        running = peak = 0

        async def wait(self, job_id, **kwargs):
            Tracking.running += 1
            Tracking.peak = max(Tracking.peak, Tracking.running)
            await real_sleep(0.01)
            Tracking.running -= 1
            return "Succeeded", True

    monkeypatch.setattr(asyncio, "sleep", no_stagger)
    monkeypatch.setattr("fedcat.core.engine.periodic_progress", lambda _msg: contextlib.nullcontext())
    platform = Tracking()
    cfg = _cfg(tmp_path, step1_pheno_group_size=1)
    cfg.cohorts[BioBank.UKB].max_concurrent_jobs = 2
    del cfg.cohorts[BioBank.AOU]
    disc = _disc(cfg, tmp_path, Phase.STEP1, platform)
    disc.jobs = _adopt(disc, disc._build_step1_jobs(BioBank.UKB, cfg.cohorts[BioBank.UKB]))

    asyncio.run(disc._run_phase(reset_jobs=None, reset_jobs_list=[]))

    assert len(platform.submitted) == 5
    assert Tracking.peak == 2
    assert all(j.state == JobState.COMPLETED for j in disc.jobs)


# ── workflow and script plumbing ──────────────────────────────────────────────


@pytest.mark.parametrize("wdl", ["regenie_step2_single_variants_aou.wdl", "regenie_step2_aggregates_aou.wdl"])
def test_step2_wdls_scatter_groups_and_pass_the_tag(wdl):
    path = _WDLS / wdl
    if not path.is_file():
        pytest.skip("AoU WDLs not available")
    text = path.read_text()
    for needle in (
        "Array[String] pheno_groups",
        "Array[String] step1_run_prefixes",
        "Array[String] group_labels",
        "Array[Array[Array[String]]] units",
        "scatter(unit in units[idx])",
        "pheno_group = group_labels[idx],",
        '"~{pheno_group}"',
        "maxRetries: max_retries",
    ):
        assert needle in text, needle


def test_step1_wdl_scatters_sets_and_uploads_its_log_last():
    path = _WDLS / "regenie_step1_scatter_aou.wdl"
    if not path.is_file():
        pytest.skip("AoU WDLs not available")
    text = path.read_text()
    assert "scatter(stratum in strata_by_set[set_idx])" in text
    assert text.index("gcloud storage cp step1_out/*") < text.index("gcloud storage cp ${log_file}")
    assert 'gcloud storage objects describe "${bucket}/${step1_directory}/${log_file}" >/dev/null 2>&1' in text


def _bash(script: str, cwd: Path) -> subprocess.CompletedProcess:
    """Run a bash script with the shared workflow functions sourced."""
    return subprocess.run(
        ["bash", "-c", f"set -euo pipefail; source {_SHARED}; {script}"],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,  # tests read the return code
    )


@pytest.fixture
def pred_lists(tmp_path):
    """A directory of two step 1 sets' prediction lists."""
    if shutil.which("bash") is None or not _SHARED.is_file():
        pytest.skip("bash or shared scripts not available")
    (tmp_path / "s1_Continuous_RINT_Step1_eur_pred.list").write_text("A /w/s1_1.loco.gz\nB /w/s1_2.loco.gz\n")
    (tmp_path / "s2_Continuous_RINT_Step1_eur_pred.list").write_text("C /w/s2_1.loco.gz\n")
    return tmp_path


def test_pred_lists_of_several_sets_are_merged(pred_lists):
    result = _bash("merge_step1_pred_lists s1,s2 Continuous RINT eur A,C merged", pred_lists)
    assert result.returncode == 0, result.stderr
    assert (pred_lists / "merged").read_text().split("\n")[:3] == [
        "A /w/s1_1.loco.gz",
        "B /w/s1_2.loco.gz",
        "C /w/s2_1.loco.gz",
    ]


def test_a_phenotype_without_predictions_fails_the_task(pred_lists):
    """regenie itself would only warn and drop it."""
    result = _bash("merge_step1_pred_lists s1 Continuous RINT eur A,C merged", pred_lists)
    assert result.returncode != 0
    assert "no step1 predictions for: C" in result.stderr


def test_step1_sources_lists_one_glob_per_prefix(pred_lists):
    result = _bash("step1_sources gs://b/dir s1,s2 Continuous RINT eur", pred_lists)
    assert result.stdout.split() == [
        "gs://b/dir/s1_Continuous_RINT_Step1_eur_*",
        "gs://b/dir/s2_Continuous_RINT_Step1_eur_*",
    ]


def test_oversized_wb_inputs_fail_with_advice(tmp_path):
    inputs = tmp_path / "inputs.json"
    inputs.write_text(json.dumps({"wf.units": ["x" * WB_MAX_INPUTS_BYTES]}))
    client = WbCliClient(output_bucket_id="bucket-id")
    with pytest.raises(CromwellClientError, match="max_tasks_per_workflow"):
        asyncio.run(client.submit_workflow(workflow_name="wf", inputs_json_path=inputs))
