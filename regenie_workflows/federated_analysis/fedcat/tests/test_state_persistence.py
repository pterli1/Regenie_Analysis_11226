"""State persistence: the write queue's ordering and shutdown, and state files matching their names."""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

import polars as pl
import pytest
from attrs import define

from fedcat.core.conditioning import ConditioningStep
from fedcat.core.engine import CondOutcome
from fedcat.core.logging import STATE_WRITE_QUEUE, StateWriteQueue, write_text_atomic
from fedcat.core.persistence import PersistentState
from fedcat.core.structures import CondType, StepState, make_conv
from fedcat.tests.test_pipeline_integration_guards import _make_aou_cfg, _make_disc_run


def test_same_key_writes_run_in_order_behind_a_slow_one():
    order: list[str] = []

    async def _drive() -> None:
        queue = StateWriteQueue(min_interval=0)
        queue.start()
        release = asyncio.Event()

        async def _slow() -> None:
            await release.wait()
            order.append("old")

        async def _fast() -> None:
            order.append("new")

        queue.enqueue("k", _slow, max_delay=0)
        await asyncio.sleep(0.05)  # the worker starts the slow write
        queue.enqueue("k", _fast, max_delay=0)
        flushed = asyncio.create_task(queue.flush("k"))
        await asyncio.sleep(0.05)
        assert order == []  # the newer write waits for the one already running
        release.set()
        await flushed
        await queue.stop()

    asyncio.run(_drive())
    assert order == ["old", "new"]


def test_stop_waits_for_a_write_already_started():
    done: list[bool] = []

    async def _drive() -> None:
        queue = StateWriteQueue(min_interval=0)
        queue.start()

        async def _write() -> None:
            await asyncio.sleep(0.1)
            done.append(True)

        queue.enqueue("k", _write, max_delay=0)
        await asyncio.sleep(0.05)  # started, not finished
        await queue.stop()

    asyncio.run(_drive())
    assert done == [True]


def test_pending_writes_coalesce_to_the_latest():
    ran: list[int] = []

    async def _drive() -> None:
        queue = StateWriteQueue(min_interval=3600)
        queue.start()

        def _write(i: int):
            async def _fn() -> None:
                ran.append(i)

            return _fn

        queue.enqueue("k", _write(0), max_delay=0)
        await asyncio.sleep(0.05)  # first write runs; the next ones are rate-limited
        queue.enqueue("k", _write(1))
        queue.enqueue("k", _write(2))
        await queue.stop()

    asyncio.run(_drive())
    assert ran == [0, 2]


def test_write_text_atomic_creates_parents_and_leaves_no_temp_file(tmp_path):
    path = tmp_path / "missing" / "state.json"
    write_text_atomic(path, "{}")
    assert path.read_text() == "{}"
    assert [p.name for p in path.parent.iterdir()] == ["state.json"]


@define
class _Stepped(PersistentState):
    """Persistent state whose file is named by its step, like a conditioning step."""

    root: Path
    step: int = 0

    def _state_file(self) -> Path:
        return self.root / f"step{self.step}" / "state.json"

    def _state_queue_key(self) -> str:
        return "stepped"


def _saved_step(path: Path) -> int:
    """Return the step saved in a state file."""
    return json.loads(path.read_text())["step"]


def test_a_pending_write_lands_in_the_file_for_the_state_it_holds(tmp_path):
    async def _drive() -> None:
        STATE_WRITE_QUEUE.start()
        try:
            obj = _Stepped(root=tmp_path)
            obj.persist_state()
            obj.step = 1  # moves on before the write runs
            await obj.flush_state()
        finally:
            await STATE_WRITE_QUEUE.stop()

    asyncio.run(_drive())
    assert not (tmp_path / "step0" / "state.json").exists()
    assert _saved_step(tmp_path / "step1" / "state.json") == 1


def test_flushing_before_moving_on_keeps_each_step_file_final(tmp_path):
    async def _drive() -> None:
        STATE_WRITE_QUEUE.start()
        try:
            obj = _Stepped(root=tmp_path)
            obj.persist_state()
            await obj.flush_state()
            obj.step = 1
            obj.persist_state()
            await obj.flush_state()
        finally:
            await STATE_WRITE_QUEUE.stop()

    asyncio.run(_drive())
    assert _saved_step(tmp_path / "step0" / "state.json") == 0
    assert _saved_step(tmp_path / "step1" / "state.json") == 1


def _completed_step0(tmp_path) -> ConditioningStep:
    """Aggregate step 0, done: chr2 has a single lead (so it passes), chr1 two (so it carries on)."""
    cfg = _make_aou_cfg(tmp_path)
    cfg.chromosomes_to_run = ["1", "2"]
    step = ConditioningStep.from_discovery_config(cfg, cfg.pheno_names[0], CondType.AGGREGATE)
    step.step_path.mkdir(parents=True)
    ids = ["SET_A.mask.singleton", "SET_B.mask.singleton", "SET_C.mask.singleton"]
    n = len(ids)
    pl.DataFrame(
        {
            "ID": ids,
            "LOG10P_IN": [20.0] * n,
            "A1FREQ_IN": [0.001] * n,
            "BETA_IN": [1.0] * n,
            "SE_IN": [0.1] * n,
            "TEST": ["BURDEN"] * n,
            "CHROM": ["1", "1", "2"],
            "START": [1] * n,
            "END": [2] * n,
            "A1FREQ_OUT": [0.001] * n,
            "BETA_OUT": [1.0] * n,
            "SE_OUT": [0.1] * n,
            "LOG10P_OUT": [20.0] * n,
            "STEP": [0] * n,
            "JUMPER": [False] * n,
        }
    ).write_csv(step.step_path / "summary_output", separator="\t")
    step.testing_ids = {"1": ["SET_A.mask", "SET_B.mask"], "2": ["SET_C.mask"]}
    step.state = StepState.COMPLETED
    return step


def _unchanged(step: ConditioningStep) -> bool:
    """Return whether a step is still the completed step 0 the test built."""
    return (
        step.step_number == 0
        and step.state == StepState.COMPLETED
        and step.passing_ids == {"1": [], "2": []}
        and step.testing_ids == {"1": ["SET_A.mask", "SET_B.mask"], "2": ["SET_C.mask"]}
    )


def test_increment_into_an_existing_step_changes_nothing(tmp_path):
    step = _completed_step0(tmp_path)
    (step.run_path / "agg_cond_step1").mkdir()

    with pytest.raises(FileExistsError):
        step.increment()
    assert _unchanged(step)


def test_increment_that_fails_to_write_the_next_step_changes_nothing(tmp_path, monkeypatch):
    step = _completed_step0(tmp_path)

    def _fail(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(pl.DataFrame, "write_csv", _fail)
    with pytest.raises(OSError):
        step.increment()
    assert _unchanged(step)
    assert not (step.run_path / "agg_cond_step1").exists()


def test_increment_advances_and_records_the_passing_lead(tmp_path):
    step = _completed_step0(tmp_path)
    step.increment()
    assert step.step_number == 1 and step.state == StepState.READY
    assert step.passing_ids == {"1": [], "2": ["SET_C.mask"]}
    assert (step.step_path / "tophits_input").read_text().count("\n") == 2


def test_resume_refuses_a_conditioning_state_file_for_another_step(tmp_path):
    cfg = _make_aou_cfg(tmp_path)
    disc = _make_disc_run(cfg)
    disc.cond_step = 0
    pheno = cfg.pheno_names[0]
    later = ConditioningStep.from_discovery_config(cfg, pheno, CondType.SINGLEVAR, step_number=1)
    step0_file = ConditioningStep.from_discovery_config(cfg, pheno, CondType.SINGLEVAR).state_path
    write_text_atomic(step0_file, make_conv().dumps(later, indent=4))

    async def _drive() -> CondOutcome:
        STATE_WRITE_QUEUE.start()
        try:
            return await disc._run_conditioning(CondType.SINGLEVAR)
        finally:
            await STATE_WRITE_QUEUE.stop()

    assert asyncio.run(_drive()) is CondOutcome.FAILED
    assert disc.cond_status[CondType.SINGLEVAR] == {pheno: CondOutcome.FAILED}


def test_rewriting_final_results_moves_the_previous_ones_to_old(tmp_path, monkeypatch):
    import fedcat.core.conditioning as cond_mod

    cfg = _make_aou_cfg(tmp_path)
    step = ConditioningStep.from_discovery_config(cfg, cfg.pheno_names[0], CondType.AGGREGATE, step_number=2)
    step.step_path.mkdir(parents=True)
    (step.step_path / "summary_output").write_text("")
    end_steps: list[int] = []
    uploads = iter([False, True, True])

    def _fake_generate(*, run_path, cond_type, end_step, **_kwargs):
        end_steps.append(end_step)
        results_dir = run_path / f"{cond_type.value}_cond_results"
        results_dir.mkdir()  # fails if the previous results were left in place
        (results_dir / "attempt").write_text(str(len(end_steps)))
        return results_dir

    async def _fake_upload(self, results_dir):
        return next(uploads)

    monkeypatch.setattr(cond_mod, "generate_final_output", _fake_generate)
    monkeypatch.setattr(ConditioningStep, "_upload_final_results_to_cloud", _fake_upload)

    assert asyncio.run(step.write_final_results()) is False  # upload failed, so a resume retries
    assert asyncio.run(step.write_final_results()) is True
    assert asyncio.run(step.write_final_results()) is True
    assert end_steps == [2, 2, 2]
    assert (step.run_path / "agg_cond_results" / "attempt").read_text() == "3"
    assert (step.run_path / "agg_cond_results_old" / "attempt").read_text() == "2"
