"""Tests for resuming a run: config comparison, conditioning step autodetection, --overwrite-config,
and the workflows and resources jobs take from the config."""

from __future__ import annotations

import pytest

from fedcat.core.engine import (
    DiscoveryRun,
    adopt_run_controls,
    compare_configs,
    config_differences,
    format_config_differences,
)
from fedcat.core.structures import BioBank, CondType, DiscoveryConfig, MetaAnalysisConfig, RunEnv
from fedcat.tests.conftest import cohort_raw, make_config


def _make_cfg(**overrides) -> DiscoveryConfig:
    """AoU and UKB, eur only, with short workflow names so job inputs read e.g. ``s1.mem_gb``."""
    short = ["s1", "s2s", "s2a", "cs", "ca"]
    workflows = dict(zip(cohort_raw("aou")["workflows"], short))
    cohorts = {
        "aou": cohort_raw(
            "aou",
            workflows={k: f"{v}.wdl" for k, v in workflows.items()},
            step1_directory="step1",
            root_results_directory="results",
        ),
        "ukb": cohort_raw(
            "ukb",
            workflows={k: f"app-{i}" for i, k in enumerate(workflows, 1)},
            step1_directory="s1",
            root_results_directory="ukb_res",
        ),
    }
    fields = {
        "meta_analysis_config": MetaAnalysisConfig(
            location=RunEnv.AOU, meta_workflow="meta.wdl", subset_workflow="subset.wdl"
        ),
        "pheno_names": ["PHENO_A"],
        "local_run_path": "/tmp/cond",
        "run_prefix": "test",
        "cond_suffix": "sfx",
        "pgs_prefix": None,
    }
    return make_config(cohorts, **(fields | overrides))


def test_identical_configs_have_no_differences():
    assert compare_configs(_make_cfg(), _make_cfg()) == []


def test_differences_are_reported_by_path_old_to_new():
    persisted = _make_cfg(cond_p_agg_retest=1e-6, cond_suffix="x")
    current = _make_cfg(cond_p_agg_retest=1e-5, cond_suffix="y")
    current.cohorts[BioBank.AOU].workflows.regenie_step1 = "new_step1.wdl"
    current.cohorts[BioBank.AOU].pheno_file = "different.tsv"
    diffs = config_differences(current, persisted)
    assert {path for path, _, _ in diffs} == {
        "cond_p_agg_retest",
        "cond_suffix",
        "cohorts.aou.workflows.regenie_step1",
        "cohorts.aou.pheno_file",
    }
    assert ("cond_p_agg_retest", 1e-6, 1e-5) in diffs
    assert "  cond_p_agg_retest: 1e-06 → 1e-05" in format_config_differences(diffs).splitlines()


def test_derived_values_are_not_differences():
    current = _make_cfg()
    current.cohorts[BioBank.AOU].output_location = "something/else"
    assert compare_configs(current, _make_cfg()) == []


def test_run_controls_are_not_differences_and_follow_the_current_config():
    persisted = _make_cfg(cond_max_steps=5)
    current = _make_cfg(cond_max_steps=10)
    assert compare_configs(current, persisted) == []
    adopt_run_controls(persisted, current)
    assert persisted.cond_max_steps == 10


@pytest.mark.parametrize(
    "entries, cond_type, expected",
    [
        ([], CondType.SINGLEVAR, 0),
        (["singles_cond_step0/"], CondType.SINGLEVAR, 0),
        (["singles_cond_step0/", "singles_cond_step1/", "singles_cond_step3/"], CondType.SINGLEVAR, 3),
        (["agg_cond_step0/", "agg_cond_step2/"], CondType.AGGREGATE, 2),
        # other types, files and unrelated folders are ignored
        (
            [
                "agg_cond_step5/",
                "singles_cond_step99",
                "random/",
                "singles_cond_step_bad/",
                "singles_cond_step2/",
            ],
            CondType.SINGLEVAR,
            2,
        ),
    ],
)
def test_autodetect_cond_step(tmp_path, entries, cond_type, expected):
    for entry in entries:
        (tmp_path / entry.rstrip("/")).mkdir() if entry.endswith("/") else (tmp_path / entry).touch()
    assert DiscoveryRun._autodetect_cond_step(tmp_path, cond_type) == expected
    assert DiscoveryRun._autodetect_cond_step(tmp_path / "missing", cond_type) == 0


# ===========================================================================
# adopt_config (--overwrite-config)
# ===========================================================================


class TestAdoptConfig:
    """A run persists local paths fixed when it was created; ``--overwrite-config``
    must re-derive them, or a moved ``local_run_path`` keeps reading the old one."""

    def _run(self, cfg: DiscoveryConfig, tmp_path) -> DiscoveryRun:
        """Create a discovery run for `cfg`."""
        return DiscoveryRun(
            config=cfg,
            run_path=cfg.discovery_run_path,
            variant_density_blockfile=tmp_path / "blocks.tsv",
        )

    def test_run_path_follows_local_run_path(self, tmp_path):
        old, new = tmp_path / "old", tmp_path / "new"
        disc = self._run(_make_cfg(local_run_path=str(old)), tmp_path)
        assert disc.run_path == old / "test_discovery_run"

        disc.adopt_config(_make_cfg(local_run_path=str(new)))
        assert disc.run_path == new / "test_discovery_run"

    def test_run_path_follows_run_prefix(self, tmp_path):
        disc = self._run(_make_cfg(local_run_path=str(tmp_path)), tmp_path)
        disc.adopt_config(_make_cfg(local_run_path=str(tmp_path), run_prefix="test_female"))
        assert disc.run_path == tmp_path / "test_female_discovery_run"

    def test_blockfiles_under_the_old_path_are_rebased(self, tmp_path):
        old, new = tmp_path / "old", tmp_path / "new"
        disc = self._run(_make_cfg(local_run_path=str(old)), tmp_path)
        disc.run_blockfile = old / "blocks" / "run.tsv"
        disc.variant_density_blockfile = tmp_path / "shared" / "density.tsv"

        disc.adopt_config(_make_cfg(local_run_path=str(new)))
        assert disc.run_blockfile == new / "blocks" / "run.tsv"
        # Outside local_run_path: an explicitly supplied file is left where it is.
        assert disc.variant_density_blockfile == tmp_path / "shared" / "density.tsv"

    def test_jobs_are_marked_stale(self, tmp_path):
        from fedcat.core.job import Job

        disc = self._run(_make_cfg(local_run_path=str(tmp_path / "old")), tmp_path)
        disc.jobs = [Job(biobank=BioBank.AOU, description="j", input_dict={"a": 1})]
        disc.adopt_config(_make_cfg(local_run_path=str(tmp_path / "new")))
        assert disc.jobs[0]._input_stale


class TestCondStepAdoptConfig:
    """A conditioning step adopting a new config with `--overwrite-config`."""

    def _step(self, cfg: DiscoveryConfig):
        """Create a singles conditioning step for `cfg`."""
        from fedcat.core.conditioning import ConditioningStep

        return ConditioningStep.from_discovery_config(cfg, "PHENO_A", CondType.SINGLEVAR)

    def test_paths_and_snapshots_follow_the_new_config(self, tmp_path):
        old, new = tmp_path / "old", tmp_path / "new"
        step = self._step(_make_cfg(local_run_path=str(old)))
        assert step.run_path == old / "test_PHENO_A_singlevar_conditioning_sfx"

        step.adopt_config(_make_cfg(local_run_path=str(new), run_prefix="test_male"))
        assert step.base_path == new
        assert step.run_prefix == "test_male"
        assert step.run_path == new / "test_male_PHENO_A_singlevar_conditioning_sfx"
        assert step.input_tarball_path == new / "PHENO_A_singles_cond_step0_subset.tar.gz"

    def test_tarball_outside_the_run_path_is_left_alone(self, tmp_path):
        old, new = tmp_path / "old", tmp_path / "new"
        elsewhere = tmp_path / "elsewhere" / "subset.tar.gz"
        step = self._step(_make_cfg(local_run_path=str(old)))
        step.input_tarball_path = elsewhere

        step.adopt_config(_make_cfg(local_run_path=str(new)))
        assert step.input_tarball_path == elsewhere


# ===========================================================================
# Workflow + resources follow the config
# ===========================================================================


class TestWorkflowFollowsConfig:
    """``target`` and the job's resources are derived from the current config on every
    input generation, so ``--overwrite-config`` moves jobs onto the new workflow/size."""

    def _job(self, parent, op: str, bb: BioBank):
        """Return a job of `op` on `bb`, linked to `parent`."""
        import weakref

        from fedcat.core.job import Job

        job = Job(biobank=bb, description=op, recipe={"op": op})
        job._parent = weakref.ref(parent)
        return job

    def _run(self, cfg: DiscoveryConfig, tmp_path) -> DiscoveryRun:
        """Create a discovery run for `cfg`."""
        return DiscoveryRun(
            config=cfg, run_path=cfg.discovery_run_path, variant_density_blockfile=tmp_path / "blocks.tsv"
        )

    def _applied(self, parent, op: str, bb: BioBank):
        """Return a job of `op` with the workflow and resources from its parent's config applied."""
        job = self._job(parent, op, bb)
        job._apply_config_workflow()
        return job

    def test_wdl_targets_are_stems_and_dx_targets_are_applet_ids(self, tmp_path):
        disc = self._run(_make_cfg(local_run_path=str(tmp_path)), tmp_path)
        cases = [
            ("step1", BioBank.AOU, "s1"),
            ("step2_singles", BioBank.AOU, "s2s"),
            ("step2_agg", BioBank.AOU, "s2a"),
            ("meta", BioBank.AOU, "meta"),
            ("step1", BioBank.UKB, "app-1"),
            ("step2_singles", BioBank.UKB, "app-2"),
            ("step2_agg", BioBank.UKB, "app-3"),
        ]
        for op, bb, expected in cases:
            assert self._applied(disc, op, bb).target == expected, op

    def test_conditioning_target_follows_the_step_type(self, tmp_path):
        from fedcat.core.conditioning import ConditioningStep

        cfg = _make_cfg(local_run_path=str(tmp_path))
        for cond_type, aou_target, ukb_target in [
            (CondType.SINGLEVAR, "cs", "app-4"),
            (CondType.AGGREGATE, "ca", "app-5"),
        ]:
            step = ConditioningStep.from_discovery_config(cfg, "PHENO_A", cond_type)
            assert self._applied(step, "cond", BioBank.AOU).target == aou_target
            assert self._applied(step, "cond", BioBank.UKB).target == ukb_target

    def test_unknown_op_is_left_alone(self, tmp_path):
        disc = self._run(_make_cfg(local_run_path=str(tmp_path)), tmp_path)
        job = self._job(disc, "not_an_op", BioBank.UKB)
        job.target, job.instance_type = "kept", "mem3_ssd1_v2_x2"
        job._apply_config_workflow()
        assert (job.target, job.instance_type) == ("kept", "mem3_ssd1_v2_x2")

    def test_dx_jobs_get_the_default_instance_type_when_unconfigured(self, tmp_path):
        from fedcat.core.structures import JobKind
        from fedcat.platforms.dnanexus.platform import DEFAULT_INSTANCE_TYPES

        disc = self._run(_make_cfg(local_run_path=str(tmp_path)), tmp_path)
        for op, kind in [("step1", JobKind.STEP1), ("step2_singles", JobKind.STEP2_SINGLES)]:
            assert self._applied(disc, op, BioBank.UKB).instance_type == DEFAULT_INSTANCE_TYPES[kind]

    def test_configured_resources_reach_each_platform(self, tmp_path):
        from fedcat.core.conditioning import ConditioningStep
        from fedcat.core.structures import JobKind, JobResources

        cfg = _make_cfg(local_run_path=str(tmp_path))
        cfg.cohorts[BioBank.UKB].workflows.resources = {
            JobKind.STEP1: JobResources(instance_type="mem3_ssd1_v2_x32"),
            JobKind.COND_SINGLES: JobResources(instance_type="mem3_ssd1_v2_x64"),
        }
        cfg.cohorts[BioBank.AOU].workflows.resources = {
            JobKind.STEP1: JobResources(n_cpu=16, mem_gb=128),
            JobKind.COND_SINGLES: JobResources(mem_gb=96),
        }
        disc = self._run(cfg, tmp_path)
        step = ConditioningStep.from_discovery_config(cfg, "PHENO_A", CondType.SINGLEVAR)

        assert self._applied(disc, "step1", BioBank.UKB).instance_type == "mem3_ssd1_v2_x32"
        assert self._applied(step, "cond", BioBank.UKB).instance_type == "mem3_ssd1_v2_x64"
        # AoU sizes are workflow inputs, prefixed with the workflow name.
        assert self._applied(disc, "step1", BioBank.AOU).input_dict == {"s1.n_cpu": 16, "s1.mem_gb": 128}
        # Only what the config sets is passed; the rest keeps the WDL's own defaults.
        assert self._applied(step, "cond", BioBank.AOU).input_dict == {"cs.mem_gb": 96}

    def test_configured_retries_and_priority_reach_each_platform(self, tmp_path):
        from fedcat.core.structures import JobKind, JobResources

        cfg = _make_cfg(local_run_path=str(tmp_path))
        cfg.cohorts[BioBank.AOU].workflows.resources = {
            JobKind.STEP1: JobResources(max_retries=2, preemptible=0),
        }
        cfg.cohorts[BioBank.UKB].workflows.resources = {JobKind.STEP1: JobResources(priority="low")}
        disc = self._run(cfg, tmp_path)

        assert self._applied(disc, "step1", BioBank.AOU).input_dict == {
            "s1.max_retries": 2,
            "s1.preemptible_tries": 0,
        }
        assert self._applied(disc, "step1", BioBank.UKB).priority == "low"
        # unconfigured kinds keep the priority the job was created with
        assert self._applied(disc, "step2_singles", BioBank.UKB).priority == "normal"

    def test_raised_size_reaches_a_retried_job(self, tmp_path):
        """The OOM path: raise the size in the config, rerun with --overwrite-config."""
        import asyncio

        from fedcat.core.structures import JobKind, JobResources

        disc = self._run(_make_cfg(local_run_path=str(tmp_path)), tmp_path)
        aou_job = self._job(disc, "step1", BioBank.AOU)
        aou_job.recipe["units"] = [[0, "eur"]]
        ukb_job = self._job(disc, "step1", BioBank.UKB)
        disc.jobs = [aou_job, ukb_job]
        asyncio.run(aou_job.generate_input())
        assert "s1.mem_gb" not in aou_job.input_dict  # WDL default until configured

        bigger = _make_cfg(local_run_path=str(tmp_path))
        bigger.cohorts[BioBank.AOU].workflows.resources = {JobKind.STEP1: JobResources(mem_gb=128)}
        bigger.cohorts[BioBank.UKB].workflows.resources = {
            JobKind.STEP1: JobResources(instance_type="mem3_ssd1_v2_x32")
        }
        disc.adopt_config(bigger)
        assert aou_job._input_stale

        asyncio.run(aou_job.generate_input())
        assert aou_job.input_dict["s1.mem_gb"] == 128
        assert self._applied(disc, "step1", BioBank.UKB).instance_type == "mem3_ssd1_v2_x32"

    def test_cost_estimate_uses_the_configured_size(self, tmp_path):
        """A job sized up in the config is priced at that size, not the WDL default."""
        from fedcat.core.job import Job
        from fedcat.platforms.gcs_cromwell import GcsCromwellCfg, GcsCromwellPlatform

        platform = GcsCromwellPlatform(GcsCromwellCfg(bucket="gs://b"))
        target = "regenie_step2_aggregates_aou"  # 8 cpu / 64 GB in the WDL
        default = platform.task_hourly_cost(target)
        job = Job(biobank=BioBank.AOU, description="j", target=target, input_dict={f"{target}.mem_gb": 128})

        assert job.resource_overrides == {"mem_gb": 128}
        assert platform.task_hourly_cost(target, **job.resource_overrides) > default

    def test_configured_disk_reaches_the_workflow_and_the_estimate(self, tmp_path):
        """The pwas case: a smaller disk for a subset run, priced as configured."""
        from fedcat.core.job import Job
        from fedcat.core.structures import JobKind, JobResources
        from fedcat.platforms.gcs_cromwell import GcsCromwellCfg, GcsCromwellPlatform

        cfg = _make_cfg(local_run_path=str(tmp_path))
        cfg.cohorts[BioBank.AOU].workflows.resources = {
            JobKind.STEP2_SINGLES: JobResources(disk="200 SSD"),
        }
        disc = self._run(cfg, tmp_path)
        job = self._job(disc, "step2_singles", BioBank.AOU)
        job.recipe["units"] = []
        job._apply_config_workflow()
        assert job.input_dict == {"s2s.disk": "200 SSD"}

        # Priced from the configured disk, not the WDL's 1000 HDD default.
        platform = GcsCromwellPlatform(GcsCromwellCfg(bucket="gs://b"))
        target = "regenie_step2_single_variants_aou"
        priced = Job(
            biobank=BioBank.AOU, description="j", target=target, input_dict={f"{target}.disk": "200 SSD"}
        )
        assert priced.resource_overrides == {"disk": "200 SSD"}
        assert platform.task_hourly_cost(target, **priced.resource_overrides) < platform.task_hourly_cost(target)

    def test_every_job_kind_names_a_config_field(self):
        """JobKind values are the config field holding that job's workflow — the
        invariant workflow_field() (and the `resources` config keys) rely on."""
        import attrs

        from fedcat.core.structures import CohortWorkflows, JobKind, MetaAnalysisConfig

        cohort_fields = {f.name for f in attrs.fields(CohortWorkflows)}
        meta_fields = {f.name for f in attrs.fields(MetaAnalysisConfig)}
        for kind in JobKind:
            expected = meta_fields if kind in (JobKind.META, JobKind.STEP0_SUBSET) else cohort_fields
            assert kind.value in expected, kind

    def test_workflow_field_reads_each_kind(self, tmp_path):
        from fedcat.core.structures import JobKind, workflow_field

        cfg = _make_cfg(local_run_path=str(tmp_path))
        found = {kind: workflow_field(cfg, BioBank.UKB, kind) for kind in JobKind}
        assert found == {
            JobKind.STEP1: "app-1",
            JobKind.STEP2_SINGLES: "app-2",
            JobKind.STEP2_AGG: "app-3",
            JobKind.META: "meta.wdl",
            JobKind.STEP0_SUBSET: "subset.wdl",
            JobKind.COND_SINGLES: "app-4",
            JobKind.COND_AGG: "app-5",
        }


def test_restarting_from_an_earlier_phase_sets_later_state_aside(tmp_path):
    """Later phases' state is moved to .old at once, so a resume before they are reached starts from the earlier phase."""
    from fedcat.__main__ import _set_aside_later_phases
    from fedcat.core.structures import Phase

    for phase in (Phase.STEP2_SINGLES, Phase.META, Phase.COND_SINGLES):
        (tmp_path / f"discovery_state_{phase.value}.json").write_text("{}")
    _set_aside_later_phases(tmp_path, Phase.STEP2_SINGLES)
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(
        [
            f"discovery_state_{Phase.STEP2_SINGLES.value}.json",
            f"discovery_state_{Phase.META.value}.json.old",
            f"discovery_state_{Phase.COND_SINGLES.value}.json.old",
        ]
    )
