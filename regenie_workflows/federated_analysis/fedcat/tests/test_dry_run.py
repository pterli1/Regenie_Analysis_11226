"""Dry runs of ``fedcat discover``: the jobs each phase builds and the state it saves."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

import fedcat.__main__ as cli
from fedcat.__main__ import _load_discovery_state
from fedcat.core.engine import DiscoveryRun
from fedcat.core.structures import BioBank, JobState, MetaMode, Phase, make_conv
from fedcat.platforms.gcs_cromwell.buckets import WorkspaceBucket
from fedcat.tests.conftest import cohort_raw


def make_min_config(tmp_dir: Path, *, cohorts: dict | None = None, **overrides) -> Path:
    """Write a minimal valid config file (AoU eur+afr and UKB eur unless *cohorts* is given) and return its path."""
    config = {
        "local_env": "aou",
        "meta_analysis_config": {
            "location": "aou",
            "meta_workflow": "aou_regenie_workflows/wdls/meta_analysis_discovery_aou.wdl",
            "subset_workflow": "aou_regenie_workflows/wdls/step0_subsetting_only_aou.wdl",
            "schema_directory": "wgs_annotations/v2026.2",
        },
        "cohorts": cohorts or {"aou": cohort_raw("aou", ["eur", "afr"]), "ukb": cohort_raw("ukb")},
        "pheno_names": ["TEST_PHENO"],
        "trait_type": "Continuous",
        "trait_transform": "RINT",
        "chromosomes_to_run": ["1"],
        "chromosomes_no_pgs": ["1"],
        "local_run_path": str(tmp_dir / "cond_base"),
        "run_prefix": "testrun",
        "aou_bucket": "dummy-bucket",
        "cond_suffix": "testsuffix",
        "pgs_prefix": "XX_SCORE_CHR",
    } | overrides
    cfg_path = tmp_dir / "config.json"
    cfg_path.write_text(json.dumps(config, indent=2))
    return cfg_path


_FAKE_BLOCKS_HEADER = "CHROM\tSTART\tEND\n"
_FAKE_BLOCKS_ROWS = "1\t1\t100000\n1\t100001\t200000\n"


def _write_fake_blocks_file(dest_dir: Path) -> Path:
    """Create a minimal variant-blocks TSV in *dest_dir* so step2 job-building works offline."""
    p = dest_dir / "variant_density_blocks.tsv"
    p.write_text(_FAKE_BLOCKS_HEADER + _FAKE_BLOCKS_ROWS)
    return p


async def _noop_async(*_args, **_kwargs) -> None:
    """No-op async coroutine for mocking async setup functions in tests."""


def _run_cli_dry_run(
    cfg_path: Path,
    phase: str,
    *,
    mock_step2: bool = False,
    tmp_dir: Path | None = None,
    extra_args: list[str] | None = None,
) -> None:
    """Invoke the unified CLI in dry-run mode, simulating ``sys.argv``.

    Pass ``mock_step2=True`` (and ``tmp_dir``) for step2 phases that need
    variant-block file download to be bypassed, and ``extra_args`` for any
    further CLI flags under test.
    """
    argv_backup = sys.argv.copy()
    sys.argv = [
        "python -m fedcat",
        "discover",
        "--config",
        str(cfg_path),
        "--start-phase",
        phase,
        "--dry-run",
        "--skip-validation",
        "--logging-level",
        "NOFILE",
        *(extra_args or []),
    ]
    try:
        with (
            patch(
                "fedcat.platforms.gcs_cromwell.buckets.list_controlled_buckets",
                return_value=[WorkspaceBucket(resource_id="dummy-bucket-id", bucket_name="dummy-bucket")],
            ),
            patch("fedcat.platforms.dnanexus.platform.setup_dxpy", _noop_async),
            patch("fedcat.platforms.dnanexus.platform.get_file_id", side_effect=lambda x, _proj: x),
            patch("fedcat.platforms.dnanexus.platform.get_file_path", side_effect=lambda x, _proj: x),
        ):
            if mock_step2:
                assert tmp_dir is not None, "tmp_dir required when mock_step2=True"
                _write_fake_blocks_file(tmp_dir)
                orig_dir = os.getcwd()
                os.chdir(tmp_dir)
                try:
                    cli.cli()
                finally:
                    os.chdir(orig_dir)
            else:
                cli.cli()
    finally:
        sys.argv = argv_backup


def _load_disc_obj(run_dir: Path, phase: str, conv) -> DiscoveryRun:
    """Load a DiscoveryRun from the canonical state JSON or the latest snapshot."""
    state_file = run_dir / f"discovery_state_{phase}.json"
    if state_file.exists():
        return conv.loads(state_file.read_text(), DiscoveryRun)

    # Dry-run path: only a snapshot file is written.
    snap_dir = run_dir / "snapshots"
    snaps = sorted(snap_dir.glob(f"snapshot_{phase}_dry_run_no_overwrite_*.json"))
    assert snaps, f"No dry-run snapshot found for phase '{phase}' in {snap_dir}"
    return conv.loads(snaps[-1].read_text(), DiscoveryRun)


def _dry_run(tmp_path: Path, phase: str, **config) -> DiscoveryRun:
    """Dry-run *phase* on a minimal config with *config* overrides and load the state it saved."""
    cfg_path = make_min_config(tmp_path, **config)
    _run_cli_dry_run(cfg_path, phase, mock_step2=phase.startswith("step2"), tmp_dir=tmp_path)
    run_prefix = json.loads(cfg_path.read_text())["run_prefix"]
    return _load_disc_obj(tmp_path / "cond_base" / f"{run_prefix}_discovery_run", phase, make_conv())


def _input_values(job, name: str) -> list:
    """Return a job's input values named `name`, whatever their workflow prefix."""
    return [v for k, v in job.input_dict.items() if k.rsplit(".", 1)[-1] == name]


def _run_prefix_vals(job) -> list[str]:
    """Output prefixes on a job; a step1 WDL lists one per phenotype set."""
    return _input_values(job, "run_prefix") + [p for v in _input_values(job, "run_prefixes") for p in v]


def _step1_run_prefix_vals(job) -> list[str]:
    """Step1 prefixes a job reads; a step2 WDL lists comma-joined prefixes per phenotype group."""
    joined = _input_values(job, "step1_run_prefix") + [
        p for v in _input_values(job, "step1_run_prefixes") for p in v
    ]
    return [p for v in joined for p in v.split(",")]


def test_step1(tmp_path):
    disc = _dry_run(tmp_path, "step1")
    assert disc.current_phase is Phase.STEP1 and disc.meta_mode is MetaMode.ALL
    assert disc.config.aou_bucket_resource_id == "dummy-bucket-id"  # resolved from the workspace's buckets
    assert all(j.state in (JobState.READY, JobState.EMPTY) and "testrun" in j.description for j in disc.jobs)
    aou = [j for j in disc.jobs if j.biobank is BioBank.AOU]
    ukb = [j for j in disc.jobs if j.biobank is BioBank.UKB]
    assert len(ukb) == 1  # one per UKB stratum
    assert aou and all("bucket/path/aou_pheno.tsv.gz" in " ".join(map(str, j.input_dict.values())) for j in aou)
    assert all(_run_prefix_vals(j) and set(_run_prefix_vals(j)) == {"testrun"} for j in disc.jobs)
    # the saved state loads back unchanged
    conv = make_conv()
    again = conv.loads(conv.dumps(disc), DiscoveryRun)
    assert [j.description for j in again.jobs] == [
        j.description for j in disc.jobs
    ] and again.config == disc.config


def test_step2(tmp_path):
    """Step2 jobs read step1 output under run_prefix by default; only aggregate jobs take the AAF bins."""
    singles = _dry_run(tmp_path, "step2_singles", aggregate_freqs=["0.001", "0.01"])
    assert singles.jobs and singles.current_phase is Phase.STEP2_SINGLES
    assert all(set(_step1_run_prefix_vals(j)) == {"testrun"} for j in singles.jobs)
    assert not any(_input_values(j, "aaf_bins") for j in singles.jobs)
    expected = singles._expected_step2_output_lists(is_singles=True)
    names = [n for ns in expected.values() for n in ns]
    assert names and all(n.endswith(".output_filelist") for n in names)

    aggs = _dry_run(tmp_path, "step2_aggregates", aggregate_freqs=["0.001", "0.01"])
    assert aggs.jobs and aggs.current_phase is Phase.STEP2_AGG
    assert all(_input_values(j, "aaf_bins") == ["0.001,0.01"] for j in aggs.jobs)
    assert {j.description for j in aggs.jobs}.isdisjoint(j.description for j in singles.jobs)


def test_step1_run_prefix_names_step1_output_and_is_read_by_step2(tmp_path):
    step1 = _dry_run(tmp_path, "step1", step1_run_prefix="prevstep1")
    assert step1.jobs and all(set(_run_prefix_vals(j)) == {"prevstep1"} for j in step1.jobs)
    step2 = _dry_run(tmp_path, "step2_singles", step1_run_prefix="prevstep1")
    for j in step2.jobs:
        assert set(_step1_run_prefix_vals(j)) == {"prevstep1"} and set(_run_prefix_vals(j)) == {"testrun"}


@pytest.mark.parametrize("rule", [True, False])
def test_meta_and_step0_subset_pass_the_acat_rule_and_bins(tmp_path, rule):
    meta = _dry_run(tmp_path, "meta_analysis", aggregate_freqs=["0.001", "0.01"], acat_singles_stratum_rare=rule)
    assert meta.jobs and all(_input_values(j, "aaf_bins") == ["0.001,0.01"] for j in meta.jobs)
    assert all(_input_values(j, "acat_singles_stratum_rare") == [rule] for j in meta.jobs)
    # subset jobs check the meta ACAT was built with the run's rule
    subset = _dry_run(tmp_path, "step0_subset", acat_singles_stratum_rare=rule)
    assert subset.jobs and all(_input_values(j, "acat_singles_stratum_rare") == [rule] for j in subset.jobs)


@pytest.mark.parametrize("phase", ["conditioning_singles", "conditioning_aggregates"])
def test_conditioning_phases_build_no_discovery_jobs(tmp_path, phase):
    """Conditioning steps build their own jobs as they run."""
    disc = _dry_run(tmp_path, phase, cond_suffix="my_special_suffix")
    assert disc.jobs == [] and disc.current_phase is Phase(phase)
    assert disc.config.cond_suffix == "my_special_suffix"


def test_aou_only_config(tmp_path):
    disc = _dry_run(tmp_path, "step1", cohorts={"aou": cohort_raw("aou", ["eur", "afr"])}, run_prefix="aouonly")
    assert disc.jobs and all(j.biobank is BioBank.AOU for j in disc.jobs)


def test_state_predating_acat_singles_rule_keeps_legacy_behaviour(tmp_path):
    """A run saved before acat_singles_stratum_rare resumes on the rule its step-0 ACAT used."""
    _dry_run(tmp_path, "step1")
    snap = max((tmp_path / "cond_base" / "testrun_discovery_run" / "snapshots").glob("snapshot_step1_*.json"))
    raw = json.loads(snap.read_text())
    assert raw["config"]["acat_singles_stratum_rare"] is True
    del raw["config"]["acat_singles_stratum_rare"]
    legacy = tmp_path / "legacy_state.json"
    legacy.write_text(json.dumps(raw))
    assert _load_discovery_state(legacy, make_conv()).config.acat_singles_stratum_rare is False
    # a state that records the field is loaded as saved
    assert _load_discovery_state(snap, make_conv()).config.acat_singles_stratum_rare is True
