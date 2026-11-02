"""Tests for the preflight checks.

Covers the egress check, the order of the check scopes and stopping after an error,
`--skip-validation`, and remote file errors (malformed reference or not found).
"""

from __future__ import annotations

import asyncio

import pytest

from fedcat.core.checks import (
    AUTH,
    CONFIG,
    REMOTE_FILES,
    WORKFLOWS,
    CheckContext,
    RunFlags,
    Severity,
    check_cohort_annotations,
    check_cond_suffix,
    check_egress_policy,
    check_meta_annotations,
    check_meta_config_present,
    check_meta_mode,
    check_remote_files,
    check_step0_between_meta_and_cond,
    run_checks,
)
from fedcat.core.structures import BioBank, MetaMode, Phase, RunEnv
from fedcat.platforms.base import PathSyntaxError
from fedcat.tests.conftest import FakePlatform
from fedcat.tests.test_pipeline_integration_guards import _make_aou_cfg


def _ctx(cfg, platform, **kw):
    """Return a check context for an AoU-only config, with optional phases, meta mode and flags."""
    return CheckContext(
        config=cfg,
        platforms={BioBank.AOU: platform},
        phases_to_run=kw.pop("phases", [Phase.STEP1]),
        meta_mode=kw.pop("meta_mode", MetaMode.ALL),
        flags=kw.pop("flags", RunFlags()),
    )


# ── egress ────────────────────────────────────────────────────────────────────


def test_egress_check_passes_for_aou_env_aou(tmp_path):
    cfg = _make_aou_cfg(tmp_path)  # local_env=AOU
    findings = asyncio.run(check_egress_policy(_ctx(cfg, FakePlatform())))
    assert findings == []


def test_egress_check_flags_forbidden_cohort_with_wrong_env(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    cfg.local_env = RunEnv.LOCAL  # AoU forbids egress => fedcat must run with local_env AOU
    findings = asyncio.run(check_egress_policy(_ctx(cfg, FakePlatform())))
    assert len(findings) == 1
    assert findings[0].severity is Severity.ERROR
    assert findings[0].scope == CONFIG
    assert "egress" in findings[0].label


# ── runner: scope ordering & flags ─────────────────────────────────────────────


def test_run_checks_short_circuits_on_auth_error(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    platform = FakePlatform(kind="gcs_cromwell", connect_error="dxpy not set up", exists_fn=lambda _p: True)
    findings = asyncio.run(run_checks(_ctx(cfg, platform)))

    assert any(f.scope == AUTH and f.severity is Severity.ERROR for f in findings)
    # Workflows + remote_files must be skipped after an auth error.
    assert not any(f.scope in (WORKFLOWS, REMOTE_FILES) for f in findings)
    assert platform.exists_seen == []


def test_run_checks_skips_workflows_and_remote_when_skip_validation(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    platform = FakePlatform(kind="gcs_cromwell", exists_fn=lambda _p: True)
    findings = asyncio.run(run_checks(_ctx(cfg, platform, flags=RunFlags(skip_validation=True))))

    assert not any(f.scope in (WORKFLOWS, REMOTE_FILES) for f in findings)
    assert platform.exists_seen == []  # remote existence never probed


def test_run_checks_full_pass_is_clean(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    platform = FakePlatform(kind="gcs_cromwell", exists_fn=lambda _p: True)
    findings = asyncio.run(run_checks(_ctx(cfg, platform)))
    assert findings == []
    assert platform.exists_seen  # remote existence was probed


# ── remote error classification ────────────────────────────────────────────────


def test_remote_check_reports_not_found(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    platform = FakePlatform(kind="gcs_cromwell", exists_fn=lambda p: "shared_input" not in p)
    findings = asyncio.run(check_remote_files(_ctx(cfg, platform)))
    msgs = [f.message for f in findings]
    assert any("not found" in m for m in msgs)
    assert not any("malformed" in m for m in msgs)


def test_remote_check_distinguishes_malformed_reference(tmp_path):
    cfg = _make_aou_cfg(tmp_path)

    def _exists(p: str) -> bool:
        if "shared_input" in p:
            raise PathSyntaxError("not a valid reference")
        return True

    platform = FakePlatform(kind="gcs_cromwell", exists_fn=_exists)
    findings = asyncio.run(check_remote_files(_ctx(cfg, platform)))
    assert any("malformed reference" in f.message for f in findings)


# ── phase-aware config requirements ───────────────────────────────────────────


def _labels(findings) -> list[str]:
    """Return the labels of findings."""
    return [f.label for f in findings]


def test_singles_only_run_needs_no_meta_config_or_annotations(tmp_path):
    """A step1 + step2_singles run should not require any aggregate/meta config."""
    cfg = _make_aou_cfg(tmp_path)
    cfg.meta_analysis_config = None
    cfg.cohorts[BioBank.AOU].regenie_annotations = []
    ctx = _ctx(cfg, FakePlatform(), phases=[Phase.STEP1, Phase.STEP2_SINGLES])

    assert asyncio.run(check_meta_config_present(ctx)) == []
    assert asyncio.run(check_cohort_annotations(ctx)) == []
    assert asyncio.run(check_meta_annotations(ctx)) == []


def test_cohort_annotations_required_for_aggregate_phases(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    cfg.cohorts[BioBank.AOU].regenie_annotations = []

    for phase in (Phase.STEP2_AGG, Phase.COND_AGG):
        findings = asyncio.run(check_cohort_annotations(_ctx(cfg, FakePlatform(), phases=[phase])))
        assert _labels(findings) == ["aou.regenie_annotations"], phase
        assert findings[0].severity is Severity.ERROR


def test_cohort_annotations_satisfied_when_present(tmp_path):
    cfg = _make_aou_cfg(tmp_path)  # fixture supplies one annotation triplet
    ctx = _ctx(cfg, FakePlatform(), phases=[Phase.STEP2_AGG])
    assert asyncio.run(check_cohort_annotations(ctx)) == []


def test_meta_config_required_only_for_meta_phase(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    cfg.meta_analysis_config = None

    assert asyncio.run(check_meta_config_present(_ctx(cfg, FakePlatform(), phases=[Phase.STEP1]))) == []
    findings = asyncio.run(check_meta_config_present(_ctx(cfg, FakePlatform(), phases=[Phase.META])))
    assert _labels(findings) == ["meta_analysis_config"]


def test_meta_workflow_still_required_when_meta_config_present(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    cfg.meta_analysis_config.meta_workflow = ""
    findings = asyncio.run(check_meta_config_present(_ctx(cfg, FakePlatform(), phases=[Phase.META])))
    assert _labels(findings) == ["meta_analysis_config.meta_workflow"]


def test_subset_workflow_required_only_for_step0_phase(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    cfg.meta_analysis_config.subset_workflow = ""
    assert asyncio.run(check_meta_config_present(_ctx(cfg, FakePlatform(), phases=[Phase.META]))) == []
    findings = asyncio.run(check_meta_config_present(_ctx(cfg, FakePlatform(), phases=[Phase.STEP0_SUBSET])))
    assert _labels(findings) == ["meta_analysis_config.subset_workflow"]


@pytest.mark.parametrize(
    "phases, ok",
    [
        ([Phase.META, Phase.COND_SINGLES], False),
        ([Phase.META, Phase.STEP0_SUBSET, Phase.COND_SINGLES], True),
        ([Phase.STEP0_SUBSET, Phase.COND_AGG], True),
        ([Phase.COND_AGG], True),
        ([Phase.META], True),
    ],
)
def test_step0_phase_required_between_meta_and_conditioning(tmp_path, phases, ok):
    findings = asyncio.run(
        check_step0_between_meta_and_cond(_ctx(_make_aou_cfg(tmp_path), FakePlatform(), phases=phases))
    )
    assert (findings == []) is ok


@pytest.mark.parametrize(
    "phases, expected",
    [
        ([Phase.STEP1, Phase.STEP2_SINGLES, Phase.META], MetaMode.SINGLES_ONLY),
        ([Phase.META, Phase.COND_SINGLES], MetaMode.SINGLES_ONLY),
        ([Phase.STEP2_AGG, Phase.META, Phase.COND_AGG], MetaMode.AGGREGATES_ONLY),
        ([Phase.STEP2_SINGLES, Phase.META, Phase.COND_AGG], MetaMode.ALL),
        ([Phase.META], MetaMode.ALL),
        (list(Phase), MetaMode.ALL),
    ],
)
def test_meta_mode_from_phases(phases, expected):
    assert MetaMode.from_phases(phases) is expected


@pytest.mark.parametrize("have_cond_file", [True, False])
@pytest.mark.parametrize("mode", list(MetaMode))
def test_meta_stages_match_legacy_job_flags(mode, have_cond_file):
    """``MetaMode.stages`` reproduces the flag logic previously inlined in the meta job builders."""
    stages = mode.stages(have_cond_file)
    assert stages.aggregates == (mode in (MetaMode.ALL, MetaMode.AGGREGATES_ONLY))
    assert stages.uncond_singles == (
        mode in (MetaMode.UNCOND_SINGLES_ONLY, MetaMode.ALL, MetaMode.SINGLES_ONLY)
        or (not have_cond_file and mode == MetaMode.AGGREGATES_ONLY)
    )
    assert stages.cond_singles == (
        have_cond_file
        and mode in (MetaMode.COND_SINGLES_ONLY, MetaMode.ALL, MetaMode.SINGLES_ONLY, MetaMode.AGGREGATES_ONLY)
    )
    assert stages.step0_subset == (mode != MetaMode.COND_SINGLES_ONLY)


def test_meta_mode_auto_always_passes_check(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    for phases in ([Phase.META, Phase.COND_SINGLES], [Phase.META, Phase.COND_AGG], list(Phase)):
        ctx = _ctx(cfg, FakePlatform(), phases=phases, meta_mode=MetaMode.from_phases(phases))
        assert asyncio.run(check_meta_mode(ctx)) == [], phases


def test_meta_mode_check_flags_modes_that_starve_conditioning(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    ctx = _ctx(cfg, FakePlatform(), phases=[Phase.META, Phase.COND_AGG], meta_mode=MetaMode.SINGLES_ONLY)
    assert _labels(asyncio.run(check_meta_mode(ctx))) == ["meta_mode"]

    # Not a meta run: the mode is irrelevant.
    ctx = _ctx(cfg, FakePlatform(), phases=[Phase.COND_AGG], meta_mode=MetaMode.SINGLES_ONLY)
    assert asyncio.run(check_meta_mode(ctx)) == []


def test_meta_mode_check_cond_singles_only_needs_cond_file(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    for cohort in cfg.cohorts.values():
        cohort.cond_file = None
    ctx = _ctx(cfg, FakePlatform(), phases=[Phase.META], meta_mode=MetaMode.COND_SINGLES_ONLY)
    assert _labels(asyncio.run(check_meta_mode(ctx))) == ["meta_mode"]


def test_meta_annotations_required_only_for_aggregate_meta_modes(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    cfg.meta_analysis_config.schema_directory = None
    cfg.meta_analysis_config.regenie_annotations = None

    for mode in (MetaMode.ALL, MetaMode.AGGREGATES_ONLY):
        ctx = _ctx(cfg, FakePlatform(), phases=[Phase.META], meta_mode=mode)
        assert _labels(asyncio.run(check_meta_annotations(ctx))) == ["meta_analysis_config.annotations"], mode

    for mode in (MetaMode.SINGLES_ONLY, MetaMode.COND_SINGLES_ONLY, MetaMode.UNCOND_SINGLES_ONLY):
        ctx = _ctx(cfg, FakePlatform(), phases=[Phase.META], meta_mode=mode)
        assert asyncio.run(check_meta_annotations(ctx)) == [], mode


def test_meta_annotations_reject_both_sources_regardless_of_phase(tmp_path):
    cfg = _make_aou_cfg(tmp_path)  # fixture sets schema_directory
    cfg.meta_analysis_config.regenie_annotations = cfg.cohorts[BioBank.AOU].regenie_annotations
    ctx = _ctx(cfg, FakePlatform(), phases=[Phase.STEP1])
    findings = asyncio.run(check_meta_annotations(ctx))
    assert _labels(findings) == ["meta_analysis_config.annotations"]
    assert "cannot set both" in findings[0].message


def test_cond_suffix_required_only_for_conditioning_phases(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    cfg.cond_suffix = None

    assert asyncio.run(check_cond_suffix(_ctx(cfg, FakePlatform(), phases=[Phase.STEP1]))) == []
    for phase in (Phase.COND_SINGLES, Phase.COND_AGG):
        findings = asyncio.run(check_cond_suffix(_ctx(cfg, FakePlatform(), phases=[phase])))
        assert _labels(findings) == ["cond_suffix"], phase


def test_remote_checks_skip_annotations_for_singles_only_run(tmp_path):
    """Annotation and meta-schema files are not stat-ed when no aggregate phase runs."""
    cfg = _make_aou_cfg(tmp_path)
    platform = FakePlatform(kind="gcs_cromwell", exists_fn=lambda _p: True)
    ctx = _ctx(cfg, platform, phases=[Phase.STEP1, Phase.STEP2_SINGLES])

    assert asyncio.run(check_remote_files(ctx)) == []
    # aou_annotations/* are the cohort triplets, wgs_schema/* the meta parquet schema.
    assert not any("aou_annotations" in seen or "wgs_schema" in seen for seen in platform.exists_seen)


def test_remote_checks_include_annotations_for_aggregate_phases(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    platform = FakePlatform(kind="gcs_cromwell", exists_fn=lambda _p: True)
    ctx = _ctx(cfg, platform, phases=[Phase.STEP2_AGG, Phase.META])

    assert asyncio.run(check_remote_files(ctx)) == []
    assert any("aou_annotations" in seen for seen in platform.exists_seen)
    assert any("wgs_schema" in seen for seen in platform.exists_seen)
