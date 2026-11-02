"""Tests for the "no conditioning variants in this stratum" marker.

Selection happens on the meta-analysis, so a stratum can legitimately lack the selected
signal outright. The job reports that with a ``_no_conditioning_variants`` marker and
returns unconditioned results, which are the conditional results in that case. What a
stratum must never do is condition once and then stop finding its conditioning variants.

No cloud credentials or downloads are used.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_no_conditioning_variants.py -v
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from fedcat.core.conditioning import ConditioningStep
from fedcat.core.structures import (
    BioBank,
    CondType,
    StepState,
)
from fedcat.tests.conftest import cohort_raw, make_config

CHROM = "1"
STRATUM = "eur"
BB = BioBank.AOU

MARKER_SUFFIX = "_no_conditioning_variants"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_step(tmp_path: Path, cond_type: CondType, step_number: int) -> ConditioningStep:
    """Build a ConditioningStep rooted at *tmp_path* with its results dir created."""
    cfg = make_config(cohorts={"aou": cohort_raw("aou", ["eur", "afr"])}, local_run_path=str(tmp_path))
    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", cond_type)
    step.step_number = step_number
    (step.step_path / f"regenie_output_{BB.value}").mkdir(parents=True, exist_ok=True)
    return step


def _results_dir(step: ConditioningStep) -> Path:
    """Return the step's download directory for the test cohort."""
    return step.step_path / f"regenie_output_{BB.value}"


def _job_prefix(step: ConditioningStep) -> str:
    """Return the name prefix of the test stratum's job outputs."""
    kind = "aggregate" if step.cond_type is CondType.AGGREGATE else "singlevar"
    return f"TEST_PHENO_{STRATUM}_chr{CHROM}_annot1_{kind}_conditioning_step{step.step_number}"


def _write_regenie(step: ConditioningStep) -> Path:
    """Write a non-empty regenie result.

    The aggregate checker treats a 2-line file as empty and the singles checker a
    1-line file, so three data lines is unambiguously non-empty for both.
    """
    path = _results_dir(step) / f"{_job_prefix(step)}.regenie.gz"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("CHROM GENPOS ID BETA SE LOG10P\n")
        fh.write("1 100 var1 0.1 0.01 3.0\n")
        fh.write("1 200 var2 0.2 0.02 4.0\n")
    return path


def _write_marker(step: ConditioningStep, cause: str = "NO_CONDITIONING_VARIANTS_PRESENT_IN_STRATUM") -> Path:
    """Write a `_no_conditioning_variants` marker with a cause."""
    path = _results_dir(step) / f"{_job_prefix(step)}{MARKER_SUFFIX}"
    path.write_text(f"{cause} {step.step_number} 7\n")
    return path


def _write_collinearity(step: ConditioningStep) -> Path:
    """Singles require a collinearity file unless the chrom is excused."""
    path = _results_dir(step) / f"{_job_prefix(step)}_variant_collinearity"
    path.write_text("ID MAX_RSQ\nvar1 0.1\n")
    return path


def _check(step: ConditioningStep) -> None:
    """Check the test stratum's results for the step's type."""
    if step.cond_type is CondType.AGGREGATE:
        step._check_chrom_results_aggs(BB, STRATUM, CHROM)
    else:
        step._check_chrom_results_singles(BB, STRATUM, CHROM)


def _setup_results(step: ConditioningStep, *, marker: bool, cause: str | None = None) -> None:
    """Write REGENIE results (and collinearity for singles), with a marker if `marker`."""
    _write_regenie(step)
    if step.cond_type is CondType.SINGLEVAR:
        _write_collinearity(step)
    if marker:
        _write_marker(step, cause) if cause else _write_marker(step)


# ---------------------------------------------------------------------------
# The marker is tolerated when the stratum has never conditioned
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cond_type", [CondType.AGGREGATE, CondType.SINGLEVAR])
def test_marker_tolerated_when_never_conditioned(tmp_path: Path, cond_type: CondType) -> None:
    """A stratum that never carried the signal runs unconditioned without erroring."""
    step = _make_step(tmp_path, cond_type, step_number=1)
    _setup_results(step, marker=True)

    _check(step)

    assert step.state != StepState.ERROR
    # Not excused from future steps: it still contributes its (unconditioned) results.
    assert CHROM not in step.chroms_cannot_condition[BB][STRATUM]
    # And it has still never conditioned, so a later marker stays legitimate.
    assert CHROM not in step.chroms_ever_conditioned[BB][STRATUM]


@pytest.mark.parametrize("cond_type", [CondType.AGGREGATE, CondType.SINGLEVAR])
def test_marker_is_not_an_unknown_file_type(tmp_path: Path, cond_type: CondType) -> None:
    """The suffix must be claimed before the else-branch that errors on unknown files."""
    step = _make_step(tmp_path, cond_type, step_number=1)
    _setup_results(step, marker=True)
    # A genuinely unrecognised file still errors, which proves the check above is
    # doing something rather than the whole branch being unreachable.
    _check(step)
    assert step.state != StepState.ERROR

    (_results_dir(step) / f"{_job_prefix(step)}_something_unexpected").write_text("x")
    _check(step)
    assert step.state == StepState.ERROR


def test_arm_ac_filter_cause_is_also_tolerated(tmp_path: Path) -> None:
    """Both causes take the same path; only the recorded token differs."""
    step = _make_step(tmp_path, CondType.AGGREGATE, step_number=2)
    _setup_results(step, marker=True, cause="ALL_CONDITIONING_VARIANTS_DROPPED_BY_ARM_AC_FILTER")

    _check(step)

    assert step.state != StepState.ERROR
    assert CHROM not in step.chroms_ever_conditioned[BB][STRATUM]


# ---------------------------------------------------------------------------
# Losing conditioning variants after having them is still an error
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cond_type", [CondType.AGGREGATE, CondType.SINGLEVAR])
def test_marker_errors_when_stratum_previously_conditioned(tmp_path: Path, cond_type: CondType) -> None:
    """Conditioning once and then finding nothing means inputs and genotypes diverged."""
    step = _make_step(tmp_path, cond_type, step_number=2)
    step.chroms_ever_conditioned[BB][STRATUM].add(CHROM)
    _setup_results(step, marker=True)

    _check(step)

    assert step.state == StepState.ERROR


# ---------------------------------------------------------------------------
# Recording that a stratum has conditioned
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cond_type", [CondType.AGGREGATE, CondType.SINGLEVAR])
def test_real_results_record_ever_conditioned(tmp_path: Path, cond_type: CondType) -> None:
    step = _make_step(tmp_path, cond_type, step_number=1)
    _setup_results(step, marker=False)

    _check(step)

    assert step.state != StepState.ERROR
    assert CHROM in step.chroms_ever_conditioned[BB][STRATUM]


@pytest.mark.parametrize("cond_type", [CondType.AGGREGATE, CondType.SINGLEVAR])
def test_step0_does_not_record_ever_conditioned(tmp_path: Path, cond_type: CondType) -> None:
    """Step 0 is the unconditioned baseline, so it never counts as having conditioned."""
    step = _make_step(tmp_path, cond_type, step_number=0)
    _setup_results(step, marker=False)

    _check(step)

    assert step.state != StepState.ERROR
    assert CHROM not in step.chroms_ever_conditioned[BB][STRATUM]


def test_ever_conditioned_survives_state_round_trip(tmp_path: Path) -> None:
    """The invariant is only enforceable if it persists across steps."""
    from fedcat.core.structures import make_conv

    step = _make_step(tmp_path, CondType.AGGREGATE, step_number=1)
    step.chroms_ever_conditioned[BB][STRATUM].add(CHROM)

    conv = make_conv()
    restored = conv.loads(conv.dumps(step), ConditioningStep)

    assert CHROM in restored.chroms_ever_conditioned[BB][STRATUM]


# ---------------------------------------------------------------------------
# The marker has to survive the download filter to be seen at all
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cond_type", [CondType.AGGREGATE, CondType.SINGLEVAR])
def test_marker_matched_by_download_pattern(tmp_path: Path, cond_type: CondType) -> None:
    """A suffix absent from the download regex would never reach the results dir."""
    import asyncio
    import re

    step = _make_step(tmp_path, cond_type, step_number=1)
    captured: dict[str, str] = {}

    class _FakePlatform:
        async def download(self, *, folder: str, pattern: str, dest: Path, max_workers: int) -> None:
            captured["pattern"] = pattern

    step.platforms = {BB: _FakePlatform()}
    asyncio.run(step._download_cohort(BB, STRATUM))

    marker_name = f"{_job_prefix(step)}{MARKER_SUFFIX}"
    assert re.match(captured["pattern"], marker_name), f"{marker_name!r} not matched by {captured['pattern']!r}"


# ---------------------------------------------------------------------------
# An empty aggregate build only ends conditioning when nothing else came back
# ---------------------------------------------------------------------------


def _write_empty_regenie(step: ConditioningStep, suffix: str) -> Path:
    """Header-only aggregate output: regenie ran but every mask was below --minMAC."""
    path = _results_dir(step) / f"{_job_prefix(step)}_{suffix}.regenie.gz"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write('##MASKS=<Mask1="a,b">\n')
        fh.write("CHROM GENPOS ID BETA SE LOG10P\n")
    return path


def _write_no_aggs_marker(step: ConditioningStep) -> Path:
    """Write a `_no_aggregates_remaining` marker."""
    path = _results_dir(step) / f"{_job_prefix(step)}_no_aggregates_remaining"
    path.write_text("NO_MASKS_ABOVE_MINMAC_AFTER_CONDITIONING\n")
    return path


def _write_singles_cma(step: ConditioningStep) -> Path:
    """Write a single-variant manc-cojo result."""
    path = _results_dir(step) / f"{_job_prefix(step)}_singlevars_mancojo.cma.cojo"
    path.write_text("Chr\tSNP\tbp\nvar1\t1\t100\n")
    return path


def test_no_aggs_marker_alone_stops_conditioning(tmp_path: Path) -> None:
    step = _make_step(tmp_path, CondType.AGGREGATE, step_number=1)
    _write_no_aggs_marker(step)

    _check(step)

    assert step.state != StepState.ERROR
    assert CHROM in step.chroms_cannot_condition[BB][STRATUM]


def test_no_aggs_marker_with_singles_keeps_stratum(tmp_path: Path) -> None:
    """manc-cojo with every mask below --minMAC: singles still feed ACAT."""
    step = _make_step(tmp_path, CondType.AGGREGATE, step_number=1)
    _write_no_aggs_marker(step)
    _write_singles_cma(step)

    _check(step)

    assert step.state != StepState.ERROR
    assert CHROM not in step.chroms_cannot_condition[BB][STRATUM]
    assert CHROM in step.chroms_ever_conditioned[BB][STRATUM]


def test_empty_regenie_alone_stops_conditioning(tmp_path: Path) -> None:
    step = _make_step(tmp_path, CondType.AGGREGATE, step_number=1)
    _write_empty_regenie(step, "annot1_mac10")

    _check(step)

    assert step.state != StepState.ERROR
    assert CHROM in step.chroms_cannot_condition[BB][STRATUM]


def test_empty_build_beside_results_keeps_stratum(tmp_path: Path) -> None:
    """An empty mac10 build next to a burden build with rows is not 'no aggregates'."""
    step = _make_step(tmp_path, CondType.AGGREGATE, step_number=1)
    _write_empty_regenie(step, "annot1_mac10")
    _write_regenie(step)

    _check(step)

    assert step.state != StepState.ERROR
    assert CHROM not in step.chroms_cannot_condition[BB][STRATUM]


# ---------------------------------------------------------------------------
# Conditioning jobs carry the ACAT singles rule to the scripts
# ---------------------------------------------------------------------------
