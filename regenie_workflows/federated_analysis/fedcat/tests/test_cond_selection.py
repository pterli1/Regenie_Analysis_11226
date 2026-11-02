"""Conditioning selection: a unit is only selected per chromosome when it passes the select threshold."""

from __future__ import annotations

import asyncio

import polars as pl

from fedcat.core.conditioning import ConditioningStep
from fedcat.core.io import TOPHITS_DTYPES
from fedcat.core.logging import STATE_WRITE_QUEUE
from fedcat.core.structures import CondType, StepState
from fedcat.tests.test_pipeline_integration_guards import _make_aou_cfg

# Defaults: agg retest 1e-6 (LOG10P 6), select 1e-7 (LOG10P 7); singles select 5e-8 (LOG10P ~7.3)


def _completed_step0(
    tmp_path, hits: list[tuple[str, str, float]], cond_type: CondType = CondType.AGGREGATE
) -> ConditioningStep:
    """Step 0, done, with (ID, CHROM, LOG10P_OUT) results."""
    cfg = _make_aou_cfg(tmp_path)
    cfg.chromosomes_to_run = sorted({c for _, c, _ in hits})
    step = ConditioningStep.from_discovery_config(cfg, cfg.pheno_names[0], cond_type)
    step.step_path.mkdir(parents=True)
    n = len(hits)
    pl.DataFrame(
        {
            "ID": [h[0] for h in hits],
            "LOG10P_IN": [20.0] * n,
            "A1FREQ_IN": [0.001] * n,
            "BETA_IN": [1.0] * n,
            "SE_IN": [0.1] * n,
            "TEST": ["BURDEN" if cond_type is CondType.AGGREGATE else "ADD"] * n,
            "CHROM": [h[1] for h in hits],
            "START": [1] * n,
            "END": [2] * n,
            "A1FREQ_OUT": [0.001] * n,
            "BETA_OUT": [1.0] * n,
            "SE_OUT": [0.1] * n,
            "LOG10P_OUT": [h[2] for h in hits],
            "STEP": [0] * n,
            "JUMPER": [False] * n,
        }
    ).write_csv(step.step_path / "summary_output", separator="\t")
    step.state = StepState.COMPLETED
    return step


def _next_input_chroms(step: ConditioningStep) -> dict[str, list[str]]:
    """Return the IDs of the next step's `tophits_input` per chromosome."""
    df = pl.read_csv(step.step_path / "tophits_input", separator="\t", has_header=False, schema=TOPHITS_DTYPES)
    return {c: sorted(g["ID"].to_list()) for (c,), g in df.group_by("CHROM")}


def test_chrom_with_only_retest_hits_is_not_carried_on(tmp_path):
    step = _completed_step0(
        tmp_path,
        [
            ("A.m.singleton", "1", 6.5),
            ("B.m.singleton", "1", 6.2),
            ("C.m.singleton", "2", 9.0),
            ("D.m.singleton", "2", 8.0),
        ],
    )
    step.increment()
    assert step.step_number == 1
    assert _next_input_chroms(step) == {"2": ["C.m.singleton", "D.m.singleton"]}
    assert step.passing_ids == {"1": [], "2": []}


def test_lone_selectable_unit_with_retest_hits_carries_on(tmp_path):
    step = _completed_step0(tmp_path, [("A.m.singleton", "1", 9.0), ("B.m.singleton", "1", 6.5)])
    step.increment()
    assert step.step_number == 1
    assert _next_input_chroms(step) == {"1": ["A.m.singleton", "B.m.singleton"]}
    assert step.passing_ids == {"1": []}


def test_lone_selectable_unit_with_nothing_to_retest_passes(tmp_path):
    step = _completed_step0(tmp_path, [("A.m.singleton", "1", 9.0), ("A.m.0.001", "1", 6.5)])
    step.increment()
    assert not step.required and step.step_number == 0
    assert step.passing_ids == {"1": ["A.m"]}


def test_nothing_selectable_ends_conditioning(tmp_path):
    step = _completed_step0(tmp_path, [("A.m.singleton", "1", 6.5), ("B.m.singleton", "2", 6.9)])
    step.increment()
    assert not step.required and step.step_number == 0
    assert step.passing_ids == {"1": [], "2": []}


def _collinear_singles_step(tmp_path, hits, max_rsq: dict[str, float]) -> ConditioningStep:
    """Return a completed singles step 0 with the given collinearity per variant."""
    step = _completed_step0(tmp_path, hits, CondType.SINGLEVAR)
    (step.step_path / "collinearity").mkdir()
    rows = "".join(f"{vid}\t0.1\t{r}\n" for vid, r in max_rsq.items())
    (step.step_path / "collinearity" / "all_cohorts_max_collinearity").write_text(
        f"ID\tcand_rsq\tmax_rsq\n{rows}"
    )
    return step


def test_collinear_singles_are_neither_selected_nor_retested(tmp_path):
    """A variant at or above cond_rsq_threshold (0.9) with a conditioned one leaves the next step."""
    hits = [("1:1:A:G", "1", 9.0), ("1:2:A:G", "1", 8.5), ("1:3:A:G", "1", 6.5)]
    step = _collinear_singles_step(tmp_path, hits, {"1:1:A:G": 0.1, "1:2:A:G": 0.95})  # 1:3 has none
    step.increment()
    assert _next_input_chroms(step) == {"1": ["1:1:A:G", "1:3:A:G"]}


def test_lone_selectable_single_left_by_collinearity_passes(tmp_path):
    step = _collinear_singles_step(
        tmp_path / "a", [("1:1:A:G", "1", 9.0), ("1:2:A:G", "1", 8.5)], {"1:2:A:G": 0.9}
    )
    step.increment()
    assert not step.required and step.passing_ids == {"1": ["1:1:A:G"]}


def _singles_step1(tmp_path, hits: list[tuple[str, float]]) -> ConditioningStep:
    """Singles step 1 input with (ID, LOG10P_IN) on chr1."""
    cfg = _make_aou_cfg(tmp_path)
    cfg.chromosomes_to_run = ["1"]
    step = ConditioningStep.from_discovery_config(cfg, cfg.pheno_names[0], CondType.SINGLEVAR, step_number=1)
    step.step_path.mkdir(parents=True)
    n = len(hits)
    pl.DataFrame(
        {
            "ID": [h[0] for h in hits],
            "LOG10P_IN": [h[1] for h in hits],
            "A1FREQ_IN": [0.01] * n,
            "BETA_IN": [1.0] * n,
            "SE_IN": [0.1] * n,
            "TEST": ["ADD"] * n,
            "CHROM": ["1"] * n,
            "START": [1] * n,
            "END": [1] * n,
        }
    ).write_csv(step.step_path / "tophits_input", separator="\t", include_header=False)
    step.state = StepState.SORT_INPUT
    return step


def _sort(step: ConditioningStep) -> None:
    """Run a step's SORT_INPUT handler."""

    async def _drive() -> None:
        STATE_WRITE_QUEUE.start()
        try:
            await step._sort_input()
        finally:
            await STATE_WRITE_QUEUE.stop()

    asyncio.run(_drive())


def test_sort_input_selects_the_top_log10p_unit(tmp_path):
    step = _singles_step1(tmp_path, [("1:1:A:G", 6.5), ("1:2:A:G", 8.0), ("1:3:A:G", 6.2)])
    _sort(step)
    assert step.state == StepState.GENERATE_JOBS
    assert step.passing_ids["1"] == ["1:2:A:G"]
    assert step.testing_ids["1"] == ["1:1:A:G", "1:3:A:G"]


def test_sort_input_errors_when_a_chrom_has_nothing_selectable(tmp_path):
    step = _singles_step1(tmp_path, [("1:1:A:G", 6.5), ("1:2:A:G", 6.2)])
    _sort(step)
    assert step.state == StepState.ERROR
    assert step.passing_ids["1"] == []


def test_rust_panic_in_a_step_sets_error(tmp_path, monkeypatch):
    """pyo3 raises a fedcat_rs panic as a BaseException; the step must still end in ERROR."""

    class _Panic(BaseException):
        pass

    step = _singles_step1(tmp_path, [("1:1:A:G", 8.0)])

    async def _panic(_self) -> None:
        raise _Panic("rust panicked")

    monkeypatch.setattr(ConditioningStep, "run", _panic)

    async def _drive() -> bool:
        STATE_WRITE_QUEUE.start()
        try:
            return await step.run_to_completion()
        finally:
            await STATE_WRITE_QUEUE.stop()

    assert asyncio.run(_drive()) is False
    assert step.state == StepState.ERROR
