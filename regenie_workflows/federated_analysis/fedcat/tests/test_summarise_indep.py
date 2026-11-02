"""Tests for the ``summarise-indep`` CLI and its backing helpers.

Covers `generate_final_output`, `summarise_independent_results` with detected and
overridden inputs, state files of older versions (only `passing_ids` is read), and
the command line.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_summarise_indep.py -v
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import patch

import polars as pl
import pytest

import fedcat.__main__ as cli
from fedcat.core.conditioning import generate_final_output
from fedcat.core.io import SUMMARY_OUTPUT_DTYPES
from fedcat.core.structures import CondType
from fedcat.stats.summarise_indep import (
    _locate_passing_ids_json,
    summarise_independent_results,
)

# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _row(
    *,
    id: str,
    chrom: str,
    step: int,
    test: str = "ADD",
    log10p_in: float = 7.0,
    log10p_out: float | None = None,
    jumper: bool = False,
    start: int = 100,
    end: int = 200,
) -> dict:
    """Build one summary_output row matching SUMMARY_OUTPUT_DTYPES."""
    return {
        "ID": id,
        "LOG10P_IN": log10p_in,
        "A1FREQ_IN": 0.01,
        "BETA_IN": 0.5,
        "SE_IN": 0.1,
        "TEST": test,
        "CHROM": chrom,
        "START": start,
        "END": end,
        "A1FREQ_OUT": None if log10p_out is None else 0.01,
        "BETA_OUT": None if log10p_out is None else 0.4,
        "SE_OUT": None if log10p_out is None else 0.12,
        "LOG10P_OUT": log10p_out,
        "STEP": step,
        "JUMPER": jumper,
    }


def _write_step(step_dir: Path, rows: list[dict]) -> None:
    """Write a step's `summary_output`."""
    step_dir.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows, schema=SUMMARY_OUTPUT_DTYPES).write_csv(
        step_dir / "summary_output", separator="\t", null_value="NA"
    )


def _write_step_p(step_dir: Path, rows: list[dict]) -> None:
    """Like _write_step but with legacy ``P_IN``/``P_OUT`` instead of ``LOG10P_*``."""
    step_dir.mkdir(parents=True, exist_ok=True)
    p_schema = pl.Schema(
        {
            ("P_IN" if k == "LOG10P_IN" else "P_OUT" if k == "LOG10P_OUT" else k): v
            for k, v in SUMMARY_OUTPUT_DTYPES.items()
        }
    )
    p_rows = [
        {("P_IN" if k == "LOG10P_IN" else "P_OUT" if k == "LOG10P_OUT" else k): v for k, v in r.items()}
        for r in rows
    ]
    pl.DataFrame(p_rows, schema=p_schema).write_csv(step_dir / "summary_output", separator="\t", null_value="NA")


def _build_singles_run(
    base: Path,
    *,
    extra_json_keys: dict | None = None,
    write_state_json: bool = True,
) -> Path:
    """Build a SINGLEVAR conditioning run dir with 3 step dirs + a state JSON.

    Passing units (per the conditioning STEP-rank invariant max_step == rank):
      chr1: ["1:100:A:T" (rank1, max STEP 1), "1:200:G:C" (rank2, max STEP 2)]
      chr2: ["2:300:C:G" (rank1, max STEP 1)]

    "1:200:G:C" has two rows at its max step (ADD + DOM) so all_top (max-LOG10P
    per unit) differs from all_sig.  "1:100:A:T" carries a JUMPER row.
    """
    run_dir = base / "height_singlevar_conditioning_test"
    run_dir.mkdir(parents=True)

    # step0: unconditioned (OUT cols null)
    _write_step(
        run_dir / "singles_cond_step0",
        [
            _row(id="1:100:A:T", chrom="1", step=0),
            _row(id="1:200:G:C", chrom="1", step=0),
            _row(id="2:300:C:G", chrom="2", step=0),
        ],
    )
    # step1: first selection conditioned in
    _write_step(
        run_dir / "singles_cond_step1",
        [
            _row(id="1:100:A:T", chrom="1", step=1, log10p_out=6.5, jumper=True),
            _row(id="1:200:G:C", chrom="1", step=1, log10p_out=4.0),
            _row(id="2:300:C:G", chrom="2", step=1, log10p_out=5.5),
        ],
    )
    # step2: only the rank-2 unit on chr1 survives, with two tests
    _write_step(
        run_dir / "singles_cond_step2",
        [
            _row(id="1:200:G:C", chrom="1", step=2, test="ADD", log10p_out=5.0),
            _row(id="1:200:G:C", chrom="1", step=2, test="DOM", log10p_out=9.0),
        ],
    )

    if write_state_json:
        state = {
            "passing_ids": {
                "1": ["1:100:A:T", "1:200:G:C"],
                "2": ["2:300:C:G"],
            },
        }
        if extra_json_keys:
            state.update(extra_json_keys)
        (run_dir / "singles_cond_step2" / "conditioning_state_singles_step2.json").write_text(
            json.dumps(state, indent=2)
        )

    return run_dir


def _read_tsv(path: Path) -> pl.DataFrame:
    """Read a results TSV."""
    return pl.read_csv(path, separator="\t", null_values=["NA"])


# ---------------------------------------------------------------------------
# Auto-detect happy path
# ---------------------------------------------------------------------------


class TestSummariseAutoDetect:
    """Summaries with every input detected."""

    def test_generates_expected_files(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        assert results_dir == run_dir / "singles_cond_results"
        assert results_dir.is_dir()
        assert (results_dir / "all_sig_conditioned_singles_height.tsv").is_file()
        assert (results_dir / "all_top_sig_conditioned_singles_height.tsv").is_file()
        # one trace file per passing unit
        assert (results_dir / "trace_chr1_1-100-A-T.tsv").is_file()
        assert (results_dir / "trace_chr1_1-200-G-C.tsv").is_file()
        assert (results_dir / "trace_chr2_2-300-C-G.tsv").is_file()

    def test_all_sig_row_counts_and_ids(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        all_sig = _read_tsv(results_dir / "all_sig_conditioned_singles_height.tsv")
        # 1:100 (1 row) + 1:200 (ADD+DOM = 2 rows) + 2:300 (1 row) = 4
        assert all_sig.height == 4
        assert set(all_sig["ID"].to_list()) == {"1:100:A:T", "1:200:G:C", "2:300:C:G"}

    def test_all_top_keeps_max_log10p_per_unit(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        all_top = _read_tsv(results_dir / "all_top_sig_conditioned_singles_height.tsv")
        # one row per unit → 3 rows
        assert all_top.height == 3
        assert all_top["ID"].to_list().count("1:200:G:C") == 1
        # the DOM row (LOG10P_OUT=9.0) wins over ADD (5.0)
        row = all_top.filter(pl.col("ID") == "1:200:G:C")
        assert row["LOG10P"][0] == pytest.approx(9.0)
        assert row["TEST"][0] == "DOM"

    def test_log10p_coalesces_out_over_in(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        all_sig = _read_tsv(results_dir / "all_sig_conditioned_singles_height.tsv")
        # 1:100 final row is at STEP1 with LOG10P_OUT=6.5
        r = all_sig.filter(pl.col("ID") == "1:100:A:T")
        assert r["LOG10P"][0] == pytest.approx(6.5)

    def test_jumper_note_populated(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        all_sig = _read_tsv(results_dir / "all_sig_conditioned_singles_height.tsv")
        note = all_sig.filter(pl.col("ID") == "1:100:A:T")["NOTE"][0]
        assert note == "jumped_in_step:1"
        # non-jumper unit has empty note
        note2 = all_sig.filter(pl.col("ID") == "2:300:C:G")["NOTE"][0]
        assert note2 in ("", None)


# ---------------------------------------------------------------------------
# Back-compatibility: only passing_ids is read
# ---------------------------------------------------------------------------


class TestBackCompat:
    """State files of older versions."""

    def test_extra_and_renamed_sibling_keys_ignored(self, tmp_path):
        # Simulate an old/foreign state JSON with unrelated sibling keys.
        run_dir = _build_singles_run(
            tmp_path,
            extra_json_keys={
                "some_removed_field": 123,
                "renamed_thing": {"nested": [1, 2, 3]},
                "step_number": 2,
                "config": {"this": "would break cattrs"},
            },
        )
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        assert (results_dir / "all_sig_conditioned_singles_height.tsv").is_file()

    def test_missing_passing_ids_key_raises(self, tmp_path):
        run_dir = _build_singles_run(tmp_path, write_state_json=False)
        bad_json = run_dir / "singles_cond_step2" / "conditioning_state_singles_step2.json"
        bad_json.write_text(json.dumps({"not_passing_ids": {}}))
        with pytest.raises(ValueError, match="passing_ids"):
            summarise_independent_results(run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height")


# ---------------------------------------------------------------------------
# Overrides & detection edge cases
# ---------------------------------------------------------------------------


class TestOverridesAndDetection:
    """Overridden and detected end step, state file and output directory."""

    def test_end_step_override(self, tmp_path):
        # Limit to end_step=1: 1:200:G:C never reaches step2, so its max STEP is 1,
        # matching its rank-2 last-unit allowance (i == len-1 and max_step == i).
        run_dir = _build_singles_run(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height", end_step=1
        )
        all_sig = _read_tsv(results_dir / "all_sig_conditioned_singles_height.tsv")
        # 1:200 now resolves to its STEP1 row only (one row), so total = 3
        assert all_sig.height == 3
        assert all_sig.filter(pl.col("ID") == "1:200:G:C").height == 1

    def test_end_step_skips_a_step_that_ended_before_its_meta_analysis(self, tmp_path):
        # Step 3 holds state but no summary_output, so the run ends on step 2
        run_dir = _build_singles_run(tmp_path)
        step3 = run_dir / "singles_cond_step3"
        step3.mkdir()
        (run_dir / "singles_cond_step2" / "conditioning_state_singles_step2.json").rename(
            step3 / "conditioning_state_singles_step3.json"
        )
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        assert (results_dir / "all_sig_conditioned_singles_height.tsv").is_file()

    def test_passing_ids_json_override(self, tmp_path):
        run_dir = _build_singles_run(tmp_path, write_state_json=False)
        # Put the JSON somewhere non-standard and point at it explicitly.
        alt = tmp_path / "elsewhere" / "old_state.json"
        alt.parent.mkdir(parents=True)
        alt.write_text(json.dumps({"passing_ids": {"1": ["1:100:A:T", "1:200:G:C"], "2": ["2:300:C:G"]}}))
        results_dir = summarise_independent_results(
            run_dir=run_dir,
            cond_type=CondType.SINGLEVAR,
            pheno_name="height",
            passing_ids_json=alt,
        )
        assert (results_dir / "all_sig_conditioned_singles_height.tsv").is_file()

    def test_output_dir_override(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        custom = tmp_path / "my_results"
        results_dir = summarise_independent_results(
            run_dir=run_dir,
            cond_type=CondType.SINGLEVAR,
            pheno_name="height",
            output_dir=custom,
        )
        assert results_dir == custom
        assert (custom / "all_sig_conditioned_singles_height.tsv").is_file()

    def test_no_step_dirs_raises(self, tmp_path):
        empty = tmp_path / "empty_run"
        empty.mkdir()
        with pytest.raises(FileNotFoundError, match="cond_step"):
            summarise_independent_results(run_dir=empty, cond_type=CondType.SINGLEVAR, pheno_name="height")

    def test_missing_run_dir_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            summarise_independent_results(
                run_dir=tmp_path / "does_not_exist",
                cond_type=CondType.SINGLEVAR,
                pheno_name="height",
            )


# ---------------------------------------------------------------------------
# Aggregate GRP_ID stripping
# ---------------------------------------------------------------------------


class TestAggregateGrpId:
    """Aggregate IDs grouped without their frequency bin suffix."""

    def test_grpid_suffix_stripped_for_aggregate(self, tmp_path):
        run_dir = tmp_path / "height_aggregate_conditioning_test"
        run_dir.mkdir(parents=True)
        # ID carries a mask/maf suffix that GRP_ID strips: GENE1.M1.0.001 -> GENE1.M1
        _write_step(
            run_dir / "agg_cond_step0",
            [_row(id="GENE1.M1.0.001", chrom="1", step=0)],
        )
        _write_step(
            run_dir / "agg_cond_step1",
            [_row(id="GENE1.M1.0.001", chrom="1", step=1, log10p_out=8.0)],
        )
        (run_dir / "agg_cond_step1" / "conditioning_state_agg_step1.json").write_text(
            json.dumps({"passing_ids": {"1": ["GENE1.M1"]}})
        )

        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.AGGREGATE, pheno_name="height"
        )
        all_sig = _read_tsv(results_dir / "all_sig_conditioned_agg_height.tsv")
        assert all_sig.height == 1
        assert all_sig["ID"][0] == "GENE1.M1.0.001"
        assert (results_dir / "trace_chr1_GENE1.M1.tsv").is_file()


# ---------------------------------------------------------------------------
# Helper-level behaviour (required vs final dir naming)
# ---------------------------------------------------------------------------


class TestGenerateFinalOutputHelper:
    """generate_final_output directly."""

    def test_step0_with_nothing_passed_copies_only(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        results_dir = generate_final_output(
            run_path=run_dir,
            cond_type=CondType.SINGLEVAR,
            passing_ids={"1": [], "2": []},
            pheno_name="height",
            required=False,
            end_step=0,
            output_dir=tmp_path / "step0_out",
        )
        assert results_dir is None
        assert (tmp_path / "step0_out" / "step0_summary_output").is_file()

    def test_step0_passes_are_summarised(self, tmp_path):
        # One unit per chrom passes at step 0 and conditioning ends there
        run_dir = _build_singles_run(tmp_path)
        results_dir = generate_final_output(
            run_path=run_dir,
            cond_type=CondType.SINGLEVAR,
            passing_ids={"1": ["1:100:A:T"], "2": ["2:300:C:G"]},
            pheno_name="height",
            required=False,
            end_step=0,
            output_dir=tmp_path / "step0_out",
        )
        assert results_dir == tmp_path / "step0_out"
        all_sig = _read_tsv(results_dir / "all_sig_conditioned_singles_height.tsv")
        assert sorted(all_sig["ID"].to_list()) == ["1:100:A:T", "2:300:C:G"]
        assert (results_dir / "trace_chr2_2-300-C-G.tsv").is_file()

    def test_required_uses_partial_dir_name(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        results_dir = generate_final_output(
            run_path=run_dir,
            cond_type=CondType.SINGLEVAR,
            passing_ids={"1": ["1:100:A:T", "1:200:G:C"], "2": ["2:300:C:G"]},
            pheno_name="height",
            required=True,
            end_step=2,
        )
        assert results_dir == run_dir / "singles_cond_partial_results_step2"

    def test_existing_results_dir_returns_none(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        (run_dir / "singles_cond_results").mkdir()
        results_dir = generate_final_output(
            run_path=run_dir,
            cond_type=CondType.SINGLEVAR,
            passing_ids={"1": ["1:100:A:T"]},
            pheno_name="height",
            required=False,
            end_step=2,
        )
        assert results_dir is None

    def test_invalid_end_step_returns_none(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        assert (
            generate_final_output(
                run_path=run_dir,
                cond_type=CondType.SINGLEVAR,
                passing_ids={},
                pheno_name="height",
                required=False,
                end_step=-1,
            )
            is None
        )


# ---------------------------------------------------------------------------
# CLI wiring
# ---------------------------------------------------------------------------


class TestCli:
    """The `fedcat summarise-indep` command."""

    def test_cli_dispatch_generates_summaries(self, tmp_path):
        run_dir = _build_singles_run(tmp_path)
        argv = [
            "fedcat",
            "summarise-indep",
            "--run-dir",
            str(run_dir),
            "--cond-type",
            "singles",
            "--pheno",
            "height",
        ]
        with patch.object(sys, "argv", argv):
            cli.cli()
        assert (run_dir / "singles_cond_results" / "all_sig_conditioned_singles_height.tsv").is_file()

    def test_cli_missing_run_dir_exits_nonzero(self, tmp_path):
        argv = [
            "fedcat",
            "summarise-indep",
            "--run-dir",
            str(tmp_path / "nope"),
            "--cond-type",
            "agg",
            "--pheno",
            "height",
        ]
        with patch.object(sys, "argv", argv), pytest.raises(SystemExit) as exc:
            cli.cli()
        assert exc.value.code == 1


# ---------------------------------------------------------------------------
# Legacy P (instead of LOG10P) columns
# ---------------------------------------------------------------------------


def _build_singles_run_p(base: Path) -> Path:
    """Same structure as _build_singles_run but with legacy P_IN/P_OUT columns.

    Also uses the legacy ``step{N}_config.json`` filename for passing_ids.
    """
    run_dir = base / "height_singlevar_conditioning_ptest"
    run_dir.mkdir(parents=True)
    _write_step_p(
        run_dir / "singles_cond_step0",
        [
            _row(id="1:100:A:T", chrom="1", step=0),
            _row(id="1:200:G:C", chrom="1", step=0),
        ],
    )
    _write_step_p(
        run_dir / "singles_cond_step1",
        [
            _row(id="1:100:A:T", chrom="1", step=1, log10p_out=6.5),
            _row(id="1:200:G:C", chrom="1", step=1, log10p_out=4.0),
        ],
    )
    _write_step_p(
        run_dir / "singles_cond_step2",
        [_row(id="1:200:G:C", chrom="1", step=2, log10p_out=5.0)],
    )
    # Legacy passing_ids filename: step{N}_config.json
    (run_dir / "singles_cond_step2" / "step2_config.json").write_text(
        json.dumps({"passing_ids": {"1": ["1:100:A:T", "1:200:G:C"]}})
    )
    return run_dir


class TestLegacyPColumns:
    """Older runs with P rather than LOG10P columns."""

    def test_autodetect_p_and_keep_p_in_output(self, tmp_path):
        run_dir = _build_singles_run_p(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        all_sig = pl.read_csv(
            results_dir / "all_sig_conditioned_singles_height.tsv",
            separator="\t",
            null_values=["NA"],
        )
        # Output keeps the "P" column name, not "LOG10P".
        assert "P" in all_sig.columns
        assert "LOG10P" not in all_sig.columns
        # 1:100 final at step1 → P=6.5 (coalesce P_OUT over P_IN)
        assert all_sig.filter(pl.col("ID") == "1:100:A:T")["P"][0] == pytest.approx(6.5)

    def test_trace_files_use_p_columns(self, tmp_path):
        run_dir = _build_singles_run_p(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir, cond_type=CondType.SINGLEVAR, pheno_name="height"
        )
        trace = pl.read_csv(results_dir / "trace_chr1_1-100-A-T.tsv", separator="\t", null_values=["NA"])
        assert "P_IN" in trace.columns
        assert "P_OUT" in trace.columns
        assert "LOG10P_IN" not in trace.columns

    def test_step_config_json_located(self, tmp_path):
        # Legacy step{N}_config.json must be found by auto-location.
        run_dir = _build_singles_run_p(tmp_path)
        located = _locate_passing_ids_json(run_dir, CondType.SINGLEVAR, end_step=2)
        assert located.name == "step2_config.json"

    def test_explicit_stat_col_override(self, tmp_path):
        # Default LOG10P run, but force the helper via explicit stat_col through the API.
        run_dir = _build_singles_run(tmp_path)
        results_dir = summarise_independent_results(
            run_dir=run_dir,
            cond_type=CondType.SINGLEVAR,
            pheno_name="height",
            stat_col="LOG10P",
        )
        all_sig = pl.read_csv(
            results_dir / "all_sig_conditioned_singles_height.tsv",
            separator="\t",
            null_values=["NA"],
        )
        assert "LOG10P" in all_sig.columns

    def test_cli_stat_col_flag(self, tmp_path):
        run_dir = _build_singles_run_p(tmp_path)
        argv = [
            "fedcat",
            "summarise-indep",
            "--run-dir",
            str(run_dir),
            "--cond-type",
            "singles",
            "--pheno",
            "height",
            "--stat-col",
            "P",
        ]
        with patch.object(sys, "argv", argv):
            cli.cli()
        all_sig = pl.read_csv(
            run_dir / "singles_cond_results" / "all_sig_conditioned_singles_height.tsv",
            separator="\t",
            null_values=["NA"],
        )
        assert "P" in all_sig.columns
