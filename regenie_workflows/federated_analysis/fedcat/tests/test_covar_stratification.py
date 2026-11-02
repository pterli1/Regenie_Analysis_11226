"""Covariate-stratified discovery runs.

Covers the config derivation (:meth:`DiscoveryConfig.covar_stratum_configs`), the
``--covar-stratum`` selection, the preflight invariants, and an end-to-end CLI dry-run
asserting that each covariate stratum gets its own isolated run with the
stratification inputs on its jobs.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_covar_stratification.py -v
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path
from unittest.mock import patch

import pytest

import fedcat.__main__ as cli
from fedcat.core.checks import CheckContext, RunFlags, Severity, check_covar_stratification
from fedcat.core.structures import (
    BioBank,
    DiscoveryConfig,
    MetaMode,
    Phase,
    make_conv,
)
from fedcat.tests.test_dry_run import _load_disc_obj, _run_cli_dry_run, make_min_config

STRAT = {"covar_name": "sex", "covar_values_names": {"2": "female", "1": "male"}, "also_run_combined": True}


def _config(tmp_path: Path, **overrides) -> DiscoveryConfig:
    """Structure the shared minimal test config, with top-level keys overridden."""
    raw = json.loads(make_min_config(tmp_path).read_text())
    raw.update(overrides)
    return make_conv().structure(raw, DiscoveryConfig)


# ---------------------------------------------------------------------------
# Derivation
# ---------------------------------------------------------------------------


def test_no_stratification_is_a_single_untouched_run(tmp_path):
    """Without the config block, the run is exactly what it always was."""
    cfg = _config(tmp_path)
    strata = cfg.covar_stratum_configs()

    assert [label for label, _ in strata] == [""]
    assert strata[0][1] is cfg, "unstratified config should be passed through, not copied"
    assert cfg.covar_filter == ("", "")


def test_covar_stratum_configs_derives_three_isolated_runs(tmp_path):
    cfg = _config(tmp_path, stratify_by_covar=STRAT)
    strata = cfg.covar_stratum_configs()

    assert [label for label, _ in strata] == ["", "female", "male"]
    by_label = dict(strata)

    # The combined run keeps the bare prefix so it lines up with an existing
    # unstratified run; the strata are suffixed, which isolates every path.
    assert by_label[""].run_prefix == "testrun"
    assert by_label["female"].run_prefix == "testrun_female"
    assert by_label["male"].run_prefix == "testrun_male"
    # step1 must be refit per stratum, so its prefix is suffixed too.
    assert by_label["female"].effective_step1_run_prefix == "testrun_female"
    assert by_label["male"].effective_step1_run_prefix == "testrun_male"

    assert by_label["female"].covar_filter == ("sex", "2")
    assert by_label["male"].covar_filter == ("sex", "1")
    assert by_label[""].covar_filter == ("", "")


def test_explicit_step1_run_prefix_is_suffixed_per_stratum(tmp_path):
    """A shared step1 fit cannot be reused across strata, so it is suffixed too."""
    cfg = _config(tmp_path, stratify_by_covar=STRAT, step1_run_prefix="shared_step1")
    by_label = dict(cfg.covar_stratum_configs())

    assert by_label[""].effective_step1_run_prefix == "shared_step1"
    assert by_label["female"].effective_step1_run_prefix == "shared_step1_female"
    assert by_label["male"].effective_step1_run_prefix == "shared_step1_male"


def test_stratifying_covar_is_excluded_from_covariates(tmp_path):
    """Sex is constant within a stratum, so regenie must not be given it."""
    cfg = _config(tmp_path, stratify_by_covar=STRAT, exclude_covars=["other_covar"])
    by_label = dict(cfg.covar_stratum_configs())

    assert by_label["female"].exclude_covars == ["other_covar", "sex"]
    assert "sex" in by_label["female"].exclude_covars_for("1")
    # ...but not for the combined run, where it still varies.
    assert by_label[""].exclude_covars == ["other_covar"]
    assert "sex" not in by_label[""].exclude_covars_for("1")


def test_exclude_covars_not_duplicated_when_already_listed(tmp_path):
    cfg = _config(tmp_path, stratify_by_covar=STRAT, exclude_covars=["sex"])
    assert dict(cfg.covar_stratum_configs())["female"].exclude_covars == ["sex"]


@pytest.mark.parametrize(
    ("cat_vars", "expected"),
    [
        ("seq_center,sex", "seq_center"),
        ("sex,seq_center", "seq_center"),
        ("sex", "NONE"),
        ("seq_center", "seq_center"),
        ("NONE", "NONE"),
        ("", "NONE"),
        ("seq_center, sex", "seq_center"),
    ],
)
def test_stratifying_covar_stripped_from_cat_vars(tmp_path, cat_vars, expected):
    """A name cannot be both a categorical covariate and excluded."""
    cfg = _config(tmp_path, stratify_by_covar=STRAT)
    for cohort in cfg.cohorts.values():
        cohort.cat_vars = cat_vars

    by_label = dict(cfg.covar_stratum_configs())
    for bb in cfg.cohorts:
        assert by_label["female"].cohorts[bb].cat_vars == expected
        assert by_label["male"].cohorts[bb].cat_vars == expected
        # the combined run is untouched
        assert by_label[""].cohorts[bb].cat_vars == cat_vars


def test_also_run_combined_false_drops_the_combined_run(tmp_path):
    cfg = _config(tmp_path, stratify_by_covar={**STRAT, "also_run_combined": False})
    assert [label for label, _ in cfg.covar_stratum_configs()] == ["female", "male"]


def test_derivation_does_not_mutate_the_source_config(tmp_path):
    """Derivation must be repeatable — it runs again on every resume."""
    cfg = _config(tmp_path, stratify_by_covar=STRAT)
    first = {label: c.run_prefix for label, c in cfg.covar_stratum_configs()}

    assert cfg.run_prefix == "testrun"
    assert cfg.cohorts[BioBank.AOU].cat_vars == "seq_center"
    assert {label: c.run_prefix for label, c in cfg.covar_stratum_configs()} == first


# ---------------------------------------------------------------------------
# Validation and preflight
# ---------------------------------------------------------------------------


def test_covar_with_many_levels_derives_one_run_per_level(tmp_path):
    """Any number of levels can be run, e.g. a disease state with several categories."""
    strat = {"covar_name": "t2d_status", "covar_values_names": {"0": "control", "1": "prediabetic", "2": "case"}}
    cfg = _config(tmp_path, stratify_by_covar=strat)
    by_label = dict(cfg.covar_stratum_configs())

    assert list(by_label) == ["", "control", "prediabetic", "case"]
    assert by_label["prediabetic"].run_prefix == "testrun_prediabetic"
    assert by_label["prediabetic"].covar_filter == ("t2d_status", "1")
    assert by_label["case"].covar_filter == ("t2d_status", "2")


@pytest.mark.parametrize(
    ("strat", "expected"),
    [
        ({**STRAT, "covar_name": ""}, "covar_name"),
        ({**STRAT, "covar_name": "a covar"}, "covar_name"),
        ({**STRAT, "covar_values_names": {}}, "cannot be empty"),
        ({**STRAT, "covar_values_names": {"": "female", "1": "male"}}, "empty value"),
        ({**STRAT, "covar_values_names": {"2": "female", "1": "female"}}, "duplicate names"),
        ({**STRAT, "covar_values_names": {"2": "combined", "1": "male"}}, "cannot be 'combined'"),
        ({**STRAT, "covar_values_names": {"2": "fe male", "1": "male"}}, "may only contain"),
    ],
)
def test_validate_rejects_bad_stratification(tmp_path, strat, expected):
    cfg = _config(tmp_path, stratify_by_covar=strat)
    with pytest.raises(ValueError, match=expected):
        cfg.validate()


@pytest.mark.parametrize("interaction", ["sex", "sex[1]"])
def test_validate_rejects_interaction_with_the_stratifying_covar(tmp_path, interaction):
    cfg = _config(tmp_path, stratify_by_covar=STRAT, interaction=interaction)
    with pytest.raises(ValueError, match="stratifying covariate"):
        cfg.validate()


def _findings(cfg: DiscoveryConfig):
    """Return the config-scope findings of a config."""
    ctx = CheckContext(
        config=cfg,
        platforms={},
        phases_to_run=[Phase.STEP1],
        variant_density_blockfile=None,
        meta_mode=MetaMode.ALL,
        flags=RunFlags(dry_run=True, skip_validation=True),
    )
    return asyncio.run(check_covar_stratification(ctx))


def test_preflight_passes_on_derived_configs(tmp_path):
    cfg = _config(tmp_path, stratify_by_covar=STRAT)
    for _label, derived in cfg.covar_stratum_configs():
        assert _findings(derived) == []


def test_preflight_catches_covar_left_in_the_covariates(tmp_path):
    """Guards a hand-edited config that sets covar_stratum without the derivation."""
    cfg = _config(tmp_path, stratify_by_covar=STRAT)
    female = dict(cfg.covar_stratum_configs())["female"]

    female.exclude_covars = []
    female.cohorts[BioBank.AOU].cat_vars = "seq_center,sex"
    findings = _findings(female)

    assert len(findings) == 2
    assert all(f.severity is Severity.ERROR for f in findings)
    assert {f.label for f in findings} == {"stratify_by_covar.covar_name", "aou.cat_vars"}


# ---------------------------------------------------------------------------
# --covar-stratum selection
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("requested", "expected"),
    [
        (None, ["", "female", "male"]),
        (["female"], ["female"]),
        (["combined"], [""]),
        (["male", "combined"], ["", "male"]),
    ],
)
def test_select_covar_strata(tmp_path, requested, expected):
    cfg = _config(tmp_path, stratify_by_covar=STRAT)
    assert [label for label, _ in cli._select_covar_strata(cfg, requested)] == expected


def test_select_covar_strata_requires_stratified_config(tmp_path):
    cfg = _config(tmp_path)
    with pytest.raises(SystemExit):
        cli._select_covar_strata(cfg, ["female"])


def test_select_covar_strata_rejects_a_covar_stratum_the_config_does_not_produce(tmp_path):
    cfg = _config(tmp_path, stratify_by_covar={**STRAT, "also_run_combined": False})
    with pytest.raises(SystemExit):
        cli._select_covar_strata(cfg, ["combined"])


# ---------------------------------------------------------------------------
# End-to-end dry run
# ---------------------------------------------------------------------------


def _stratified_config(tmp_path: Path, **strat_overrides) -> Path:
    """Write a minimal config with `stratify_by_covar` set, returning its path."""
    cfg_path = make_min_config(tmp_path)
    raw = json.loads(cfg_path.read_text())
    raw["stratify_by_covar"] = {**STRAT, **strat_overrides}
    cfg_path.write_text(json.dumps(raw, indent=2))
    return cfg_path


def _run_dir(cfg_path: Path, run_prefix: str) -> Path:
    """Return the discovery run directory of a run prefix."""
    return Path(json.loads(cfg_path.read_text())["local_run_path"]) / f"{run_prefix}_discovery_run"


def test_dry_run_builds_every_stratum_independently(tmp_path):
    """One invocation produces one complete, separate run per stratum."""
    cfg_path = _stratified_config(tmp_path)

    _run_cli_dry_run(cfg_path, "step1")

    conv = make_conv()
    for run_prefix, expected in [
        ("testrun", ("", "")),
        ("testrun_female", ("sex", "2")),
        ("testrun_male", ("sex", "1")),
    ]:
        run_dir = _run_dir(cfg_path, run_prefix)
        assert run_dir.is_dir(), f"no run directory for {run_prefix}"
        disc = _load_disc_obj(run_dir, "step1", conv)

        assert disc.config.run_prefix == run_prefix
        assert disc.config.covar_filter == expected
        assert disc.jobs, f"no jobs built for {run_prefix}"

        for job in disc.jobs:
            keys = {k.rsplit(".", 1)[-1]: v for k, v in job.input_dict.items()}
            if expected == ("", ""):
                # An unstratified run's payload is unchanged, so it stays
                # compatible with workflow versions that predate these inputs.
                assert "stratify_covar" not in keys and "stratify_value" not in keys
            else:
                assert (keys["stratify_covar"], keys["stratify_value"]) == expected


def test_dry_run_stratum_flag_runs_only_that_stratum(tmp_path):
    cfg_path = _stratified_config(tmp_path)

    _run_cli_dry_run(cfg_path, "step1", extra_args=["--covar-stratum", "female"])

    assert _run_dir(cfg_path, "testrun_female").is_dir()
    assert not _run_dir(cfg_path, "testrun_male").exists()
    assert not _run_dir(cfg_path, "testrun").exists()


# ---------------------------------------------------------------------------
# Parallel launch
# ---------------------------------------------------------------------------


class _FakeRun:
    """Stands in for a prepared DiscoveryRun, recording when it ran."""

    def __init__(self, config: DiscoveryConfig, order: list[str], fail: bool = False) -> None:
        self.config = config
        self._order = order
        self._fail = fail
        self.had_error = False

    async def run(self, **_kwargs) -> None:
        self._order.append(f"start:{self.config.run_prefix}")
        await asyncio.sleep(0)  # yield, so overlapping runs interleave
        if self._fail:
            raise RuntimeError("boom")
        self._order.append(f"end:{self.config.run_prefix}")


def _discover_args(cfg_path: Path, **overrides):
    """Return parsed `discover` arguments for a config, with attributes overridden."""
    parser = cli._build_parser()
    args = parser.parse_args(["discover", "--config", str(cfg_path), "--logging-level", "NOFILE"])
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


class _FakeWriteQueue:
    """Stub for the global STATE_WRITE_QUEUE singleton.

    The real one binds an ``asyncio.Event`` to the loop it is started in, so the
    repeated ``asyncio.run`` calls in these tests would strand its worker.
    """

    def __init__(self, order: list[str]) -> None:
        self._order = order

    def start(self) -> None:
        self._order.append("queue:start")

    async def stop(self) -> None:
        self._order.append("queue:stop")


def _run_discover_with_fakes(cfg_path: Path, order: list[str], failing: frozenset[str] = frozenset()):
    """Drive _discover_main with preparation stubbed out, exercising the launch path."""
    prepared: list[_FakeRun] = []

    async def _fake_prepare(_args, config_obj, _conv, _bucket_cache):
        run = _FakeRun(config_obj, order, fail=config_obj.run_prefix in failing)
        prepared.append(run)
        return run

    with (
        patch.object(cli, "_prepare_discovery_run", _fake_prepare),
        patch.object(cli, "STATE_WRITE_QUEUE", _FakeWriteQueue(order)),
    ):
        asyncio.run(cli._discover_main(_discover_args(cfg_path)))
    return prepared


def test_strata_run_concurrently(tmp_path):
    """All strata are launched together rather than one after another."""
    cfg_path = _stratified_config(tmp_path)
    order: list[str] = []

    prepared = _run_discover_with_fakes(cfg_path, order)

    assert [r.config.run_prefix for r in prepared] == ["testrun", "testrun_female", "testrun_male"]
    # The write queue wraps the whole launch, and every run starts before any
    # of them finishes, i.e. they are concurrent rather than sequential.
    assert order[0] == "queue:start"
    assert order[-1] == "queue:stop"
    assert order[1:4] == ["start:testrun", "start:testrun_female", "start:testrun_male"]
    assert set(order[4:-1]) == {"end:testrun", "end:testrun_female", "end:testrun_male"}


def test_one_failing_stratum_does_not_stop_the_others(tmp_path):
    """A stratum is an independent run; its failure must not cancel its siblings."""
    cfg_path = _stratified_config(tmp_path)
    order: list[str] = []

    with pytest.raises(SystemExit) as exit_info:
        _run_discover_with_fakes(cfg_path, order, failing={"testrun_female"})

    assert exit_info.value.code == 1, "a failed stratum makes the invocation exit non-zero"
    assert "end:testrun" in order
    assert "end:testrun_male" in order
    assert "end:testrun_female" not in order
    assert order[-1] == "queue:stop", "the write queue must still be shut down cleanly"


def test_unstratified_run_still_runs(tmp_path):
    """The single-run path is preserved for configs without stratification."""
    cfg_path = make_min_config(tmp_path)
    order: list[str] = []

    prepared = _run_discover_with_fakes(cfg_path, order)

    assert [r.config.run_prefix for r in prepared] == ["testrun"]
    assert order == ["queue:start", "start:testrun", "end:testrun", "queue:stop"]


def test_dry_run_step2_stratified_jobs_carry_the_filter(tmp_path):
    """Step2 fans out over strata/blocks; every job must carry the covariate stratum."""
    cfg_path = _stratified_config(tmp_path, also_run_combined=False)

    _run_cli_dry_run(cfg_path, "step2_singles", mock_step2=True, tmp_dir=tmp_path)

    conv = make_conv()
    for run_prefix, value in [("testrun_female", "2"), ("testrun_male", "1")]:
        disc = _load_disc_obj(_run_dir(cfg_path, run_prefix), "step2_singles", conv)
        assert disc.jobs
        for job in disc.jobs:
            keys = {k.rsplit(".", 1)[-1]: v for k, v in job.input_dict.items()}
            assert keys["stratify_covar"] == "sex"
            assert keys["stratify_value"] == value
            # the excluded sex covariate reaches regenie's --covarExcludeList
            excl = keys.get("exclcovars") or keys.get("units")
            assert "sex" in str(excl)
