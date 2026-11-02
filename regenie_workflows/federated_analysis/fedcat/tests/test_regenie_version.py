"""Tests for the ``regenie_version`` config option.

The version reaches every regenie job as a workflow/applet input; the AoU tasks
select ``/usr/local/bin/regenie-<version>`` from the fedcat_meta image and the
UKB applets their bundled ``regenie_v<version>...`` binary.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_regenie_version.py -v
"""

from __future__ import annotations

import asyncio
import json
import re
import weakref
from pathlib import Path

import pytest

from fedcat.core.engine import DiscoveryRun
from fedcat.core.job import Job
from fedcat.core.structures import DEFAULT_REGENIE_VERSION, REGENIE_VERSIONS, BioBank
from fedcat.tests.test_autodetect_and_overwrite import _make_cfg

_REPO = Path(__file__).resolve().parents[3]
_DOCKERFILE = _REPO / "federated_analysis" / "docker" / "Dockerfile"
_WDL_DIR = _REPO / "aou_regenie_workflows" / "wdls"
_APPLET_DIR = _REPO / "ukb_regenie_apps"
_REGENIE_WDLS = [
    "regenie_step1_scatter_aou.wdl",
    "regenie_step2_single_variants_aou.wdl",
    "regenie_step2_aggregates_aou.wdl",
    "regenie_conditioning_singles_aou.wdl",
    "regenie_conditioning_agg_aou.wdl",
]
_REGENIE_APPLETS = {
    "regenie_step1": "regenie_step1.sh",
    "regenie_step2_single_variants": "regenie_step2_single_variants.sh",
    "regenie_step2_aggregates": "regenie_step2_aggregates.sh",
    "regenie_conditioning_singles": "singles_cond_run_step.sh",
    "regenie_conditioning_agg": "agg_cond_run_step.sh",
}

_RECIPES = {
    ("step1", BioBank.AOU): {"op": "step1", "units": [[0, "eur"]]},
    ("step1", BioBank.UKB): {"op": "step1", "stratum": "eur"},
    ("step2_singles", BioBank.AOU): {"op": "step2_singles", "units": [[0, "eur", "1", "1", "1000"]]},
    ("step2_singles", BioBank.UKB): {
        "op": "step2_singles",
        "stratum": "eur",
        "block": ["1", "1", "1000", "None"],
    },
    ("step2_agg", BioBank.AOU): {"op": "step2_agg", "units": [[0, "eur", "1", "1", "1000"]]},
    ("step2_agg", BioBank.UKB): {"op": "step2_agg", "stratum": "eur", "block": ["1", "1", "1000", "None"]},
}


def _job_inputs(cfg, biobank: BioBank, op: str, tmp_path) -> dict[str, object]:
    """Build one job's submission payload, keyed by its bare input name."""
    disc = DiscoveryRun(
        config=cfg, run_path=cfg.discovery_run_path, variant_density_blockfile=tmp_path / "blocks.tsv"
    )
    job = Job(biobank=biobank, description=op, recipe=_RECIPES[(op, biobank)])
    job._parent = weakref.ref(disc)
    disc.jobs = [job]
    asyncio.run(job.generate_input())
    return {k.rsplit(".", 1)[-1]: v for k, v in job.input_dict.items()}


# ── config ────────────────────────────────────────────────────────────────────


def test_default_is_4_1_3(tmp_path):
    assert DEFAULT_REGENIE_VERSION == "4.1.3" and DEFAULT_REGENIE_VERSION in REGENIE_VERSIONS
    assert _make_cfg(local_run_path=str(tmp_path)).regenie_version == "4.1.3"


@pytest.mark.parametrize("version", REGENIE_VERSIONS)
def test_validate_accepts_supported_versions(version, tmp_path):
    _make_cfg(local_run_path=str(tmp_path), regenie_version=version).validate()


@pytest.mark.parametrize("version", ["4.0", "4.1.4", "v4.1", "4.1; id", ""])
def test_validate_rejects_other_versions(version, tmp_path):
    with pytest.raises(ValueError, match="regenie_version"):
        _make_cfg(local_run_path=str(tmp_path), regenie_version=version).validate()


# ── job payloads ──────────────────────────────────────────────────────────────


@pytest.mark.parametrize("op", ["step1", "step2_singles", "step2_agg"])
@pytest.mark.parametrize("biobank", [BioBank.AOU, BioBank.UKB])
def test_version_reaches_regenie_jobs(op, biobank, tmp_path):
    cfg = _make_cfg(local_run_path=str(tmp_path), regenie_version="4.1")
    assert _job_inputs(cfg, biobank, op, tmp_path)["regenie_version"] == "4.1"


# ── image and workflow plumbing ───────────────────────────────────────────────


def test_image_installs_every_supported_version():
    if not _DOCKERFILE.is_file():
        pytest.skip("Dockerfile not available")
    dockerfile = _DOCKERFILE.read_text()
    installed = set(
        re.findall(
            r"^COPY --from=\S+ /usr/local/bin/regenie /usr/local/bin/regenie-(\S+)$", dockerfile, re.MULTILINE
        )
    )
    assert installed == set(REGENIE_VERSIONS)
    assert f"ln -s regenie-{DEFAULT_REGENIE_VERSION} /usr/local/bin/regenie" in dockerfile


@pytest.mark.parametrize("name", _REGENIE_WDLS)
def test_wdl_declares_forwards_and_selects_the_version(name):
    wdl_path = _WDL_DIR / name
    if not wdl_path.is_file():
        pytest.skip("AoU WDLs not available")
    wdl = wdl_path.read_text()
    assert f'String regenie_version = "{DEFAULT_REGENIE_VERSION}"' in wdl  # workflow input
    assert "regenie_version = regenie_version," in wdl  # passed to the task
    assert re.search(r"^\s*String regenie_version$", wdl, re.MULTILINE)  # task input
    assert "ln -sf /usr/local/bin/regenie-~{regenie_version} /usr/local/bin/regenie" in wdl


@pytest.mark.parametrize(("applet", "script"), _REGENIE_APPLETS.items())
def test_applet_declares_and_selects_the_version(applet, script):
    if not (_APPLET_DIR / applet).is_dir():
        pytest.skip("UKB applets not available")
    spec = json.loads((_APPLET_DIR / applet / "dxapp.json").read_text())
    entry = next(i for i in spec["inputSpec"] if i["name"] == "regenie_version")
    assert entry["optional"] and entry["default"] == DEFAULT_REGENIE_VERSION and entry["class"] == "string"
    src = (_APPLET_DIR / applet / "src" / script).read_text()
    assert f"regenie_version=${{regenie_version:-{DEFAULT_REGENIE_VERSION}}}" in src
    assert "mv regenie_v${regenie_version}.gz_x86_64_ubuntu20_mkl /usr/local/bin/regenie" in src
