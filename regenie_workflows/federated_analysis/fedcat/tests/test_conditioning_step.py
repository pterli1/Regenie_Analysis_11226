"""A conditioning step's job building, result checks and meta-analysis, on small real files.

The meta-analysis and ACAT run for real; only the cloud is faked.
"""

from __future__ import annotations

import asyncio
import gzip
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from fedcat.core.conditioning import ConditioningStep
from fedcat.core.io import COLLINEARITY_DTYPES
from fedcat.core.structures import BioBank, CondType, StepState
from fedcat.stats.meta_analysis import z_to_log10p
from fedcat.tests.conftest import FakePlatform, cohort_raw, make_config

REGENIE_HEADER = "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA"
CMA_HEADER = "Chr\tSNP\tbp\tA1\tA2\tfreq\tb\tse\tp\tn\tbC\tbC_se\tpC"


@pytest.fixture(autouse=True)
def _offline(monkeypatch):
    """Stop conditioning steps saving state or uploading."""

    async def _noop(*_args, **_kwargs):
        return None

    monkeypatch.setattr(ConditioningStep, "persist_state", lambda *_a, **_k: None)
    monkeypatch.setattr(ConditioningStep, "_upload_meta_outputs_to_cloud", _noop)


def _step(
    tmp_path, cond_type: CondType, *, step_number: int = 1, aou_only: bool = True, **config
) -> ConditioningStep:
    """Return a conditioning step for an AoU (or two-cohort) config, with its directory created."""
    cohorts = {"aou": cohort_raw("aou", ["eur", "afr"])} if aou_only else None
    cfg = make_config(cohorts, local_run_path=str(tmp_path), **config)
    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", cond_type, step_number=step_number)
    step.step_path.mkdir(parents=True)
    return step


def _regenie(path: Path, rows: list[tuple[str, float, float, float]]) -> Path:
    """A regenie result file of (ID, A1FREQ, BETA, SE) rows."""
    lines = [REGENIE_HEADER]
    for vid, freq, beta, se in rows:
        chrom, pos = vid.split(":")[:2] if ":" in vid else ("1", "100")
        lines.append(f"{chrom} {pos} {vid} A G {freq} 1000 ADD {beta} {se} {(beta / se) ** 2} 0 NA")
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def _tophits(step: ConditioningStep, rows: list[tuple[str, float, str, str]]) -> None:
    """tophits_input of (ID, LOG10P_IN, TEST, CHROM) rows."""
    text = "".join(f"{i}\t{lp}\t0.001\t0.5\t0.1\t{t}\t{c}\t100\t300\n" for i, lp, t, c in rows)
    (step.step_path / "tophits_input").write_text(text)


def _summary(step: ConditioningStep) -> dict[tuple[str, str], dict]:
    """Return the step's `summary_output` rows by (ID, TEST)."""
    df = pl.read_csv(step.step_path / "summary_output", separator="\t", null_values="NA", infer_schema_length=0)
    return {(r["ID"], r["TEST"]): r for r in df.iter_rows(named=True)}


def _log10p(beta: float, se: float) -> float:
    """Return -log10(P) of `beta / se`."""
    return float(z_to_log10p(np.array([beta / se]))[0])


# ── Job generation ────────────────────────────────────────────────────────────


def test_jobs_cover_the_chromosomes_each_stratum_can_condition(tmp_path):
    """AoU runs one workflow over every stratum; DNAnexus one job per stratum and chromosome."""
    step = _step(tmp_path, CondType.SINGLEVAR, aou_only=False, chromosomes_to_run=["1", "2", "X"])
    step.platforms = {BioBank.AOU: FakePlatform(), BioBank.UKB: FakePlatform()}
    step.passing_ids = {"1": ["1:5:A:G"], "2": ["2:5:A:G"], "X": ["X:5:A:G"]}
    step.testing_ids = {"1": ["1:6:A:G"], "2": ["2:6:A:G"], "X": []}
    step.chroms_cannot_condition[BioBank.AOU]["afr"].add("2")
    step.chroms_cannot_condition[BioBank.UKB]["eur"].add("2")

    asyncio.run(step._generate_jobs())

    assert step.state is StepState.RUN_JOBS
    ids_file = step.step_path / "ids_to_test_step1"
    assert ids_file.read_text() == "1\t1:6:A:G\n2\t2:6:A:G\n"
    for bb, platform in step.platforms.items():
        assert platform.uploaded == [
            (ids_file, f"{step.config.output_location(bb)}/meta_analysis/{step.dir_name}/testing_ids_lists")
        ]
    by_bb = {bb: [j for j in step.jobs if j.biobank is bb] for bb in step.platforms}
    (aou,) = by_bb[BioBank.AOU]
    assert aou.description == "aou_singles_step1_chr1,2"
    assert aou.recipe["input_list"] == [("1", "1:5:A:G"), ("2", "2:5:A:G")]
    assert aou.recipe["n_tasks"] == 3  # chr1 in eur and afr, chr2 in eur only
    (ukb,) = by_bb[BioBank.UKB]
    assert (ukb.description, ukb.recipe["input_list"]) == ("eur_ukb_singles_step1_chr1", [("1", "1:5:A:G")])


def test_no_jobs_to_build_is_an_error(tmp_path):
    step = _step(tmp_path, CondType.SINGLEVAR)
    step.platforms = {BioBank.AOU: FakePlatform()}
    asyncio.run(step._generate_jobs())
    assert step.state is StepState.ERROR


# ── Result checks ─────────────────────────────────────────────────────────────


def test_results_are_checked_only_where_jobs_ran(tmp_path):
    step = _step(tmp_path, CondType.SINGLEVAR, chromosomes_to_run=["1", "2"])
    step.testing_ids = {"1": ["1:6:A:G"], "2": []}
    step.chroms_cannot_condition[BioBank.AOU]["afr"].add("1")
    results = step.step_path / "regenie_output_aou"
    results.mkdir()
    _regenie(results / "x_eur_chr1_singlevar_TEST_PHENO.regenie.gz", [("1:6:A:G", 0.01, 0.2, 0.05)])
    (results / "x_eur_chr1_variant_collinearity").write_text("ID\tcand_rsq\tmax_rsq\n")

    asyncio.run(step._check_results())
    assert step.state is StepState.META_ANALYSE
    assert step.chroms_ever_conditioned[BioBank.AOU]["eur"] == {"1"}

    (results / "x_eur_chr1_singlevar_TEST_PHENO.regenie.gz").unlink()
    step.state = StepState.CHECK_RESULTS
    asyncio.run(step._check_results())
    assert step.state is StepState.ERROR


# ── Singles meta-analysis ─────────────────────────────────────────────────────


def _singles_step(tmp_path, eur_chr1_rows) -> ConditioningStep:
    """Step 1 on chr1 and chr2, each with a passing lead, eur and afr results on chr1 and none on chr2."""
    step = _step(tmp_path, CondType.SINGLEVAR, chromosomes_to_run=["1", "2"])
    step.passing_ids = {"1": ["1:100:A:G"], "2": ["2:100:A:G"]}
    _tophits(
        step,
        [("1:100:A:G", 12, "ADD", "1"), ("1:200:A:G", 8, "ADD", "1"), ("1:300:A:G", 7, "ADD", "1")]
        + [("2:100:A:G", 9, "ADD", "2"), ("2:200:A:G", 7, "ADD", "2")],
    )
    results = step.step_path / "regenie_output_aou"
    results.mkdir()
    _regenie(results / "x_eur_chr1_singlevar_TEST_PHENO.regenie.gz", eur_chr1_rows)
    _regenie(results / "x_afr_chr1_singlevar_TEST_PHENO.regenie.gz", [("1:200:A:G", 0.01, 0.3, 0.1)])
    (results / "x_eur_chr1_variant_collinearity").write_text("ID\tcand_rsq\tmax_rsq\n1:200:A:G\t0.1\t0.2\n")
    (results / "x_afr_chr1_variant_collinearity").write_text("ID\tcand_rsq\tmax_rsq\n1:200:A:G\t0.3\t0.15\n")
    return step


def test_singles_meta_analysis_summarises_the_step(tmp_path):
    step = _singles_step(tmp_path, [("1:200:A:G", 0.01, 0.2, 0.05), ("1:300:A:G", 0.02, 0.4, 0.1)])
    asyncio.run(step._meta_analyse_singles())

    assert step.state is StepState.COMPLETED
    rows = _summary(step)
    assert float(rows[("1:300:A:G", "ADD")]["LOG10P_OUT"]) == pytest.approx(_log10p(0.4, 0.1), rel=1e-4)
    w = np.array([1 / 0.05**2, 1 / 0.1**2])
    beta, se = (w * [0.2, 0.3]).sum() / w.sum(), 1 / np.sqrt(w.sum())
    assert float(rows[("1:200:A:G", "ADD")]["LOG10P_OUT"]) == pytest.approx(_log10p(beta, se), rel=1e-4)
    # the passing leads and chr2, which no stratum returned, have no result
    for vid in ("1:100:A:G", "2:100:A:G", "2:200:A:G"):
        assert rows[(vid, "ADD")]["LOG10P_OUT"] is None
    assert {r["STEP"] for r in rows.values()} == {"1"}
    # with nothing left on chr2, no stratum conditions there again
    assert all(
        "2" in chroms and "1" not in chroms for chroms in step.chroms_cannot_condition[BioBank.AOU].values()
    )
    # collinearity is the worst over the cohorts
    coll = pl.read_csv(
        step.step_path / "collinearity" / "all_cohorts_max_collinearity",
        separator="\t",
        schema_overrides=COLLINEARITY_DTYPES,
    )
    assert coll.rows() == [("1:200:A:G", 0.3, 0.2)]


def test_singles_result_for_a_passing_lead_is_an_error(tmp_path):
    step = _singles_step(tmp_path, [("1:100:A:G", 0.01, 0.2, 0.05), ("1:200:A:G", 0.01, 0.2, 0.05)])
    asyncio.run(step._meta_analyse_singles())
    assert step.state is StepState.ERROR


# ── Aggregate meta-analysis ───────────────────────────────────────────────────


def test_mancojo_aggregate_step_uses_the_conditional_results(tmp_path):
    """From step 1 the manc-cojo engine's .cma.cojo are the results; its regenie runs are unconditional."""
    step = _step(tmp_path, CondType.AGGREGATE, cond_engine="mancojo")
    step.passing_ids = {"1": ["H.lof"]}
    _tophits(
        step,
        [("H.lof.0.001", 12, "BURDEN", "1"), ("G.lof.0.001", 8, "BURDEN", "1"), ("G.lof.0.001", 9, "ACAT", "1")],
    )
    results = step.step_path / "regenie_output_aou"
    results.mkdir()
    cma = {
        "burden": "1\tG.lof.0.001\t200\tG\tA\t0.0008\t0.9\t0.1\t1e-19\t5000\t0.5\t0.1\t5.7e-7",
        "singlevars": "1\t1:200:A:G\t200\tG\tA\t0.0005\t0.7\t0.2\t4e-4\t5000\t0.6\t0.2\t2.7e-3",
    }
    for build, row in cma.items():
        (results / f"x_eur_chr1_annot1_{build}_mancojo.cma.cojo").write_text(f"{CMA_HEADER}\n{row}\n")
    # the engine's own unconditional regenie run, which must not be read
    _regenie(results / "x_eur_chr1_annot1_burden_TEST_PHENO.regenie.gz", [("G.lof.0.001", 0.0008, 0.9, 0.1)])
    step0 = step.run_path / "agg_cond_step0"
    step0.mkdir()
    (step0 / "annotations").write_text("1:200:A:G\tG\tpLoF\n")
    (step0 / "masks").write_text("lof\tpLoF\n")

    asyncio.run(step._meta_analyse_aggs())

    assert step.state is StepState.COMPLETED
    rows = _summary(step)
    burden = float(rows[("G.lof.0.001", "BURDEN")]["LOG10P_OUT"])
    assert burden == pytest.approx(_log10p(0.5, 0.1), rel=1e-4)  # bC/bC_se, not the marginal b/se
    assert rows[("G.lof.0.001", "ACAT")]["LOG10P_OUT"] is not None
    assert rows[("H.lof.0.001", "BURDEN")]["LOG10P_OUT"] is None


def test_step0_extracts_its_subset_once(tmp_path):
    import tarfile

    tophits = tmp_path / "tophits_input"
    tophits.write_text("1:100:A:G\t7\t0.01\t0.1\t0.02\tADD\t1\t100\t100\n")
    tarball = tmp_path / "subset.tar.gz"
    with tarfile.open(tarball, "w:gz") as tar:
        tar.add(tophits, arcname="tophits_input")
    cfg = make_config({"aou": cohort_raw("aou", ["eur"])}, local_run_path=str(tmp_path))
    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", CondType.SINGLEVAR, tarball_path=tarball)

    asyncio.run(step._ready())
    assert (step.step_path / "tophits_input").read_text() == tophits.read_text()
    assert step.state == StepState.SORT_INPUT

    # A step 0 directory already there is kept, not extracted over
    (step.step_path / "tophits_input").write_text("kept\n")
    step.state = StepState.READY
    asyncio.run(step._ready())
    assert (step.step_path / "tophits_input").read_text() == "kept\n"
