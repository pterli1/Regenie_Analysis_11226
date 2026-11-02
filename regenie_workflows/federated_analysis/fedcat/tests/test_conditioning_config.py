"""Tests for the config structures, their validation and serialisation, and what conditioning steps and
jobs take from them."""

from __future__ import annotations

import weakref
from enum import Enum

import attrs
import pytest

from fedcat.core.conditioning import ConditioningStep
from fedcat.core.engine import DiscoveryRun
from fedcat.core.job import Job
from fedcat.core.structures import (
    PHASE_STRINGS_ORDERED,
    BioBank,
    CohortCfg,
    CohortWorkflows,
    CondType,
    CovarStratification,
    DiscoveryConfig,
    JobKind,
    JobResources,
    MetaMode,
    Phase,
    TraitType,
    Transform,
    make_conv,
    stratum_descrip,
)
from fedcat.tests.conftest import cohort_raw, make_config


def test_phase_order():
    assert PHASE_STRINGS_ORDERED == [
        "step1",
        "step2_singles",
        "step2_aggregates",
        "meta_analysis",
        "step0_subset",
        "conditioning_singles",
        "conditioning_aggregates",
    ]


def test_stratum_descrip():
    assert (stratum_descrip(BioBank.AOU, "eur"), stratum_descrip(BioBank.UKB, "afr")) == ("eur_aou", "afr_ukb")


# ── Output locations ──────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "transform, trait_type, folder",
    [
        ("RINT", "Continuous", "rint"),
        ("MCC", "Continuous", "mcc"),
        ("RAW", "Continuous", "raw"),
        ("RAW", "Binary", "binary"),
    ],
)
def test_output_location(tmp_path, transform, trait_type, folder):
    cfg = make_config(trait_transform=Transform(transform), trait_type=TraitType(trait_type))
    expected = {
        BioBank.AOU: f"dummy-bucket/aou_results/testrun_{folder}_wgs",
        BioBank.UKB: f"/ukb_results/testrun_{folder}_wgs",
    }
    assert {bb: cfg.output_location(bb) for bb in cfg.cohorts} == expected
    # a run stores them on its cohorts
    run = DiscoveryRun(config=cfg, run_path=tmp_path, variant_density_blockfile=tmp_path / "blocks.tsv")
    assert {bb: c.output_location for bb, c in run.config.cohorts.items()} == expected


# ── Validation ────────────────────────────────────────────────────────────────


def test_validate_fills_only_missing_cohort_sizes():
    cfg = make_config()
    cfg.cohorts[BioBank.AOU].cohort_sizes = {"eur": 9999}
    cfg.validate()
    assert cfg.cohorts[BioBank.AOU].cohort_sizes == {"eur": 9999, "afr": 100_000}


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda c: setattr(c, "pheno_names", []), "pheno_names"),
        (lambda c: setattr(c.cohorts[BioBank.AOU], "strata", []), "no strata"),
        (lambda c: setattr(c.cohorts[BioBank.AOU], "strata", ["eur_hq"]), "strata"),
        (lambda c: setattr(c.cohorts[BioBank.AOU], "strata", ["eur", "eur"]), "strata"),
        (lambda c: setattr(c.cohorts[BioBank.AOU], "strata", ["eur.1"]), "strata"),
        (lambda c: setattr(c.cohorts[BioBank.AOU], "cat_vars", ""), "cat_vars"),
        (lambda c: setattr(c.cohorts[BioBank.AOU], "cat_vars", "a, b"), "cat_vars"),
        (
            lambda c: setattr(
                c.cohorts[BioBank.AOU].workflows, "resources", {JobKind.STEP1: JobResources(preemptible=-1)}
            ),
            "preemptible cannot be negative",
        ),
        (
            lambda c: setattr(
                c.cohorts[BioBank.AOU].workflows,
                "resources",
                {JobKind.STEP2_AGG: JobResources(priority="urgent")},
            ),
            "priority must be one of",
        ),
        *[
            (lambda c, f=f: setattr(c, "aggregate_freqs", [f]), "aggregate_freqs entry")
            for f in ("0.0010", "1e-3", "1", "0.00001", ".001")
        ],
    ],
)
def test_validate_rejects(change, message):
    cfg = make_config()
    change(cfg)
    with pytest.raises(ValueError, match=message):
        cfg.validate()


def test_validate_accepts():
    cfg = make_config(cohorts={"aou": cohort_raw("aou", ["all-samples", "EUR2"])})
    cfg.aggregate_freqs = ["0.0001", "0.001", "0.01", "0.5"]
    cfg.cohorts[BioBank.AOU].workflows.resources = {
        JobKind.STEP1: JobResources(max_retries=0, preemptible=0, priority="low")
    }
    cfg.validate()
    assert cfg.cohorts[BioBank.AOU].cohort_sizes == {"all-samples": 100_000, "EUR2": 100_000}


def test_validate_reports_every_excluded_covar_name_with_a_separator():
    cfg = make_config(exclude_covars=["good", "bad,name"], exclude_covars_dict={"1": ["also bad"]})
    with pytest.raises(ValueError) as exc:
        cfg.validate()
    assert "bad,name" in str(exc.value) and "also bad" in str(exc.value)


@pytest.mark.parametrize(
    "fields, message",
    [
        ({"cond_engine": "manccojo"}, "cond_engine"),
        ({"cond_mancojo_effect_mode": "gcta"}, "cond_mancojo_effect_mode"),
        ({"cond_p_agg_select": 1e-5}, "cond_p_agg_select"),  # looser than the retest threshold
        ({"cond_p_singles_select": 1e-5}, "cond_p_singles_select"),
    ],
)
def test_config_rejects_on_construction(fields, message):
    with pytest.raises(ValueError, match=message):
        attrs.evolve(make_config(), **fields)


def test_conditioning_defaults():
    cfg = make_config()
    assert (cfg.cond_p_agg_retest, cfg.cond_p_agg_select) == (1e-6, 1e-7)
    assert (cfg.cond_p_singles_retest, cfg.cond_p_singles_select) == (1e-6, 5e-8)
    assert (cfg.cond_rsq_threshold, cfg.cond_check_coll, cfg.cond_max_steps, cfg.cond_min_arm_ac) == (
        0.9,
        True,
        -1,
        0,
    )
    assert (cfg.cond_engine, cfg.cond_mancojo_effect_mode) == ("regenie", "imputeNA")


def test_exclude_covars_stack_pgs_global_and_per_chromosome():
    cfg = make_config(
        chromosomes_to_run=["1", "X"],
        exclude_covars=["batch"],
        exclude_covars_dict={"X": ["sex"]},
        pgs_prefix="PGS_chr",
    )
    cfg.chromosomes_no_pgs = ["X"]
    assert (cfg.exclude_covars_for("1"), cfg.exclude_covars_for("X")) == ("PGS_chr1,batch", "batch,sex")
    # step1 fits every chromosome at once, so only the global exclusions apply
    assert cfg.step1_exclude_covars == "batch"
    cfg.exclude_covars, cfg.exclude_covars_dict = [], {}
    assert (cfg.exclude_covars_for("X"), cfg.step1_exclude_covars) == ("None", "None")


# ── Serialisation ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize("enum_cls", [BioBank, TraitType, Transform, CondType, Phase, MetaMode, JobKind])
def test_enums_serialise_by_value(enum_cls: type[Enum]):
    conv = make_conv()
    for member in enum_cls:
        assert conv.unstructure(member) == member.value
        assert conv.structure(member.value, enum_cls) is member


def test_config_round_trips(tmp_path):
    """A config with non-default values everywhere, and a minimal one, survive a JSON round trip."""
    conv = make_conv()
    full = make_config(
        cond_p_agg_retest=1e-5,
        cond_p_singles_retest=1e-7,
        cond_rsq_threshold=0.8,
        cond_check_coll=False,
        cond_max_steps=5,
        cond_min_arm_ac=5,
        cond_engine="mancojo",
        cond_mancojo_effect_mode="GCTA",
        exclude_covars=["batch"],
        exclude_covars_dict={"X": ["sex"]},
        stratify_by_covar=CovarStratification(covar_name="sex", covar_values_names={"F": "female", "M": "male"}),
    )
    full.cohorts[BioBank.UKB].workflows.resources = {
        JobKind.COND_SINGLES: JobResources(instance_type="mem2_ssd1_v2_x4"),
        JobKind.STEP2_AGG: JobResources(n_cpu=16, mem_gb=128),
    }
    minimal = make_config(meta_analysis_config=None, pgs_prefix=None, cond_suffix=None, chromosomes_no_pgs=[])
    for cfg in (full, minimal):
        assert conv.loads(conv.dumps(cfg), DiscoveryConfig) == cfg
    # state files hold the per-stratum chromosome sets
    sets = {BioBank.AOU: {"eur": {"1", "X"}}, BioBank.UKB: {"afr": set()}}
    assert conv.structure(conv.unstructure(sets), dict[BioBank, dict[str, set[str]]]) == sets


def test_cohort_optional_fields_default():
    raw = cohort_raw("aou")
    for key in ("cond_file", "cohort_sizes"):
        del raw[key]
    cohort = make_conv().structure(raw, CohortCfg)
    assert (cohort.cond_file, cohort.output_location, cohort.cohort_sizes) == ("", "", {})
    assert (
        make_conv()
        .structure(
            {k: f"applet-{k}" for k in CohortWorkflows.__annotations__ if k != "resources"}, CohortWorkflows
        )
        .resources
        == {}
    )


# ── ConditioningStep.from_discovery_config ────────────────────────────────────


def test_conditioning_step_from_config(tmp_path):
    cfg = make_config(local_run_path=str(tmp_path), chromosomes_to_run=["1", "2", "X"], cond_suffix="sfx")
    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", CondType.SINGLEVAR)
    assert step.config is cfg
    assert step.dir_name == "TEST_PHENO_singlevar_conditioning_sfx"
    # local runs share local_run_path between run prefixes; the cloud name, which the scripts rebuild, does not
    assert step.run_path == tmp_path / "testrun_TEST_PHENO_singlevar_conditioning_sfx"
    assert step.input_tarball_path == tmp_path / "TEST_PHENO_singles_cond_step0_subset.tar.gz"
    assert (step.pheno.name, step.pheno.trait_type, step.pheno.transform) == (
        "TEST_PHENO",
        TraitType.CONTINUOUS,
        Transform.RINT,
    )
    assert step.chromosomes_to_run == ["1", "2", "X"]
    assert step.chroms_cannot_condition == {
        BioBank.AOU: {"eur": set(), "afr": set()},
        BioBank.UKB: {"eur": set()},
    }

    agg = ConditioningStep.from_discovery_config(
        cfg, "TEST_PHENO", CondType.AGGREGATE, tarball_path=tmp_path / "t.tar.gz"
    )
    assert agg.dir_name == "TEST_PHENO_aggregate_conditioning_sfx"
    assert agg.input_tarball_path == tmp_path / "t.tar.gz"


def test_conditioning_run_paths_differ_between_covar_strata(tmp_path):
    cfg = make_config(local_run_path=str(tmp_path))
    cfg.stratify_by_covar = CovarStratification(covar_name="sex", covar_values_names={"F": "female", "M": "male"})
    paths = {
        ConditioningStep.from_discovery_config(c, "TEST_PHENO", CondType.SINGLEVAR).run_path
        for _, c in cfg.covar_stratum_configs()
    }
    assert len(paths) == 3


def test_mark_jobs_stale(tmp_path):
    disc = DiscoveryRun(
        config=make_config(), run_path=tmp_path, variant_density_blockfile=tmp_path / "blocks.tsv"
    )
    disc.jobs = [Job(biobank=BioBank.AOU, description="a"), Job(biobank=BioBank.UKB, description="b")]
    disc.mark_jobs_stale()
    assert all(j._input_stale for j in disc.jobs)


# ── Job inputs ────────────────────────────────────────────────────────────────


def _cond_job_inputs(cfg: DiscoveryConfig, cond_type: CondType, step_number: int) -> dict[BioBank, dict]:
    """Input names (without their workflow prefix) of a conditioning job on each platform."""
    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", cond_type)
    step.step_number = step_number
    recipes = {
        BioBank.AOU: {
            "op": "cond",
            "strata": ["eur", "afr"],
            "testing_ids_file_desc": "gs://b/ids",
            "input_list": [("1", "ID")],
        },
        BioBank.UKB: {
            "op": "cond",
            "stratum": "eur",
            "testing_ids_file_desc": "file-ids",
            "input_list": [("1", "ID")],
        },
    }
    out = {}
    for bb, recipe in recipes.items():
        job = Job(biobank=bb, description=f"{bb.value}_cond", recipe=recipe)
        job._parent = weakref.ref(step)
        job._gen_cond_wdl() if bb is BioBank.AOU else job._gen_cond_dx()
        out[bb] = {k.rsplit(".", 1)[-1]: v for k, v in job.input_dict.items()}
    return out


@pytest.mark.parametrize("cond_type", list(CondType))
@pytest.mark.parametrize("mode, rule", [("GCTA", True), ("imputeNA", False), ("removeNA", True)])
def test_cond_jobs_pass_engine_settings_and_step(tmp_path, cond_type, mode, rule):
    """Jobs name their outputs by fedcat's step, and hand the engine, its effect mode and the ACAT rule to the scripts."""
    cfg = make_config(
        local_run_path=str(tmp_path),
        cond_engine="mancojo",
        cond_mancojo_effect_mode=mode,
        acat_singles_stratum_rare=rule,
    )
    for bb, inputs in _cond_job_inputs(cfg, cond_type, step_number=3).items():
        assert (inputs["cond_step"], inputs["cond_engine"], inputs["cond_mancojo_effect_mode"]) == (
            3,
            "mancojo",
            mode,
        ), bb
        # only the aggregate scripts run ACAT
        assert inputs.get("acat_singles_stratum_rare") == (rule if cond_type is CondType.AGGREGATE else None), bb


def test_step1_jobs_exclude_only_covariates_excluded_on_every_chromosome(tmp_path):
    cfg = make_config(exclude_covars=["batch"], exclude_covars_dict={"1": ["local_pc"]}, pgs_prefix="PGS_chr")
    disc = DiscoveryRun(config=cfg, run_path=tmp_path, variant_density_blockfile=tmp_path / "blocks.tsv")
    recipes = {
        BioBank.AOU: {"op": "step1", "units": [[0, "eur"]]},
        BioBank.UKB: {"op": "step1", "stratum": "eur", "set": 0},
    }
    for bb, recipe in recipes.items():
        job = Job(biobank=bb, description=f"{bb.value}_step1", recipe=recipe)
        job._parent = weakref.ref(disc)
        job._gen_step1_wdl() if bb is BioBank.AOU else job._gen_step1_dx()
        assert [v for k, v in job.input_dict.items() if k.rsplit(".", 1)[-1] == "exclcovars"] == ["batch"], bb
