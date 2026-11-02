"""Tests for the variant collinearity checker (fedcat.stats.collinearity).

The core quantities are validated against an *independent* least-squares oracle
(``numpy.linalg.lstsq``) rather than a frozen copy of the implementation, so the
tests pin the mathematical meaning of the outputs:

``cand_rsq``
    R² of a candidate variant regressed on the selected set.
``max_rsq``
    Largest VIF R² across the augmented set ``selected ∪ {candidate}``.

Coverage:
  1. cand_rsq / max_rsq vs the lstsq oracle across a range of collinearity.
  2. Realistic {0,1,2}+missing genotypes vs the oracle.
  3. Degenerate cases → the collinearity sentinel (2.0).
  4. Invariance to candidate ordering.

These require only numpy + scipy (no pgenlib).
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from fedcat.stats.collinearity import (
    COLLINEAR_SENTINEL,
    candidate_collinearity,
    compute_variances,
    mean_centre,
)


# ---------------------------------------------------------------------------
# lstsq oracle
# ---------------------------------------------------------------------------
def _rsq_on(design: np.ndarray, y: np.ndarray) -> float:
    """R² of ``y`` regressed on the columns of ``design`` (both centred)."""
    if design.shape[1] == 0:
        return 0.0
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    return float(1.0 - np.sum(resid**2) / np.sum(y**2))


def _oracle_vif_rsq(g: np.ndarray, idx: list[int]) -> np.ndarray:
    """VIF R² of each variant in ``idx`` regressed on the others (centred ``g``)."""
    sub = g[idx]
    out = np.empty(len(idx))
    for i in range(len(idx)):
        others = np.delete(sub, i, axis=0).T
        out[i] = _rsq_on(others, sub[i])
    return out


def _oracle_candidate(g: np.ndarray, slct_idx: list[int], cand_i: int) -> tuple[float, float]:
    """Return `(cand_rsq, max_rsq)` of one candidate by direct regression."""
    cand_rsq = _rsq_on(g[slct_idx].T, g[cand_i])
    max_rsq = float(_oracle_vif_rsq(g, list(slct_idx) + [cand_i]).max())
    return cand_rsq, max_rsq


# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------
def _make_candidate(selected: np.ndarray, target_r2: float, rng: np.random.Generator) -> np.ndarray:
    """A candidate row with (approximately) ``target_r2`` against the selected span."""
    k, n = selected.shape
    signal = rng.standard_normal(k) @ selected
    noise = rng.standard_normal(n)
    if target_r2 <= 0.0:
        return noise
    var_sig = np.var(signal)
    scale = np.sqrt(var_sig * (1.0 - target_r2) / (target_r2 * np.var(noise)))
    return signal + scale * noise


def _centred_blocks(selected: np.ndarray, candidates: np.ndarray):
    """Stack, mean-centre together, split back; return blocks + variances + full g."""
    g = np.vstack([selected, candidates]).astype(np.float64)
    mean_centre(g)
    k = selected.shape[0]
    g_slct, g_cand = g[:k], g[k:]
    return g, g_slct, compute_variances(g_slct), g_cand, compute_variances(g_cand)


# ---------------------------------------------------------------------------
# 1. Correctness across a range of collinearity, vs the lstsq oracle
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("target_r2", [0.0, 0.1, 0.3, 0.5, 0.7, 0.85, 0.95])
def test_cand_rsq_matches_oracle(target_r2):
    rng = np.random.default_rng(20240722 + int(target_r2 * 100))
    n_samples, k = 800, 5
    selected = rng.standard_normal((k, n_samples))
    candidates = np.vstack([_make_candidate(selected, target_r2, rng) for _ in range(8)])

    g, g_slct, var_slct, g_cand, var_cand = _centred_blocks(selected, candidates)
    cand_rsq, max_rsq = candidate_collinearity(g_slct, var_slct, g_cand, var_cand)

    for j in range(candidates.shape[0]):
        exp_cand, exp_max = _oracle_candidate(g, list(range(k)), k + j)
        assert cand_rsq[j] == pytest.approx(exp_cand, abs=1e-6)
        assert max_rsq[j] == pytest.approx(exp_max, abs=1e-6)


def test_cand_rsq_monotonic_in_collinearity():
    """Higher designed collinearity yields higher reported cand_rsq."""
    rng = np.random.default_rng(7)
    n_samples, k = 600, 4
    selected = rng.standard_normal((k, n_samples))
    targets = [0.05, 0.25, 0.5, 0.75, 0.9]
    candidates = np.vstack([_make_candidate(selected, t, rng) for t in targets])

    _, g_slct, var_slct, g_cand, var_cand = _centred_blocks(selected, candidates)
    cand_rsq, _ = candidate_collinearity(g_slct, var_slct, g_cand, var_cand)
    assert np.all(np.diff(cand_rsq) > 0), cand_rsq
    assert np.all((cand_rsq > 0) & (cand_rsq < 1))


# ---------------------------------------------------------------------------
# 2. Realistic {0,1,2}+missing genotypes vs the oracle
# ---------------------------------------------------------------------------
def test_genotypes_with_missing_match_oracle():
    """Discrete hardcalls with missing (-9→nan) entries agree with the oracle.

    Both the checker and the oracle operate on the same mean-centred data (where
    missing entries have been imputed to 0), so they must agree exactly.
    """
    rng = np.random.default_rng(202)
    n_samples, k, m = 400, 6, 25
    selected = rng.integers(0, 3, size=(k, n_samples)).astype(np.float64)
    candidates = rng.integers(0, 3, size=(m, n_samples)).astype(np.float64)
    for block in (selected, candidates):
        block[rng.random(block.shape) < 0.02] = np.nan

    g = np.vstack([selected, candidates])
    mean_centre(g)
    g_slct, g_cand = g[:k], g[k:]
    cand_rsq, max_rsq = candidate_collinearity(
        g_slct, compute_variances(g_slct), g_cand, compute_variances(g_cand)
    )

    for j in range(m):
        exp_cand, exp_max = _oracle_candidate(g, list(range(k)), k + j)
        assert cand_rsq[j] == pytest.approx(exp_cand, abs=1e-6)
        assert max_rsq[j] == pytest.approx(exp_max, abs=1e-6)


# ---------------------------------------------------------------------------
# 3. Degenerate cases → sentinel
# ---------------------------------------------------------------------------
def test_perfect_pairwise_collinearity_sentinel():
    rng = np.random.default_rng(303)
    n_samples, k = 300, 4
    selected = rng.standard_normal((k, n_samples))
    # Candidate 0 duplicates a selected variant (pairwise r²=1); candidate 1 is clean.
    candidates = np.vstack([selected[2].copy(), rng.standard_normal(n_samples)])

    _, g_slct, var_slct, g_cand, var_cand = _centred_blocks(selected, candidates)
    cand_rsq, max_rsq = candidate_collinearity(g_slct, var_slct, g_cand, var_cand)

    assert cand_rsq[0] == COLLINEAR_SENTINEL
    assert max_rsq[0] == COLLINEAR_SENTINEL
    assert cand_rsq[1] < 0.9  # the clean candidate is unaffected


def test_exact_linear_combination_sentinel():
    """A candidate in the exact span of the selected set is degenerate.

    This span-collinearity has no single pairwise r² over the threshold, so it is
    caught by the Schur complement being within rounding error of zero, rather than by
    the pairwise check.
    """
    rng = np.random.default_rng(404)
    n_samples, k = 300, 5
    selected = rng.standard_normal((k, n_samples))
    candidates = (selected[0] + selected[1] - selected[3])[None, :]

    _, g_slct, var_slct, g_cand, var_cand = _centred_blocks(selected, candidates)
    cand_rsq, max_rsq = candidate_collinearity(g_slct, var_slct, g_cand, var_cand)
    assert cand_rsq[0] == COLLINEAR_SENTINEL
    assert max_rsq[0] == COLLINEAR_SENTINEL


# ---------------------------------------------------------------------------
# 4. Ordering invariance
# ---------------------------------------------------------------------------
def test_candidate_order_invariance():
    rng = np.random.default_rng(505)
    n_samples, k = 500, 5
    selected = rng.standard_normal((k, n_samples))
    candidates = np.vstack([_make_candidate(selected, t, rng) for t in [0.2, 0.5, 0.8, 0.1]])
    _, g_slct, var_slct, g_cand, var_cand = _centred_blocks(selected, candidates)

    base_rsq, base_max = candidate_collinearity(g_slct, var_slct, g_cand, var_cand)
    perm = np.array([2, 0, 3, 1])
    perm_rsq, perm_max = candidate_collinearity(g_slct, var_slct, g_cand[perm], var_cand[perm])
    np.testing.assert_allclose(perm_rsq, base_rsq[perm], atol=1e-10)
    np.testing.assert_allclose(perm_max, base_max[perm], atol=1e-10)


def test_collinearity_cli_threads_reach_blas_env(monkeypatch):
    """--threads is applied to the BLAS thread variables (before numpy loads in the CLI)."""
    from fedcat._blas_threads import BLAS_THREAD_VARS, preset_blas_threads

    for var in BLAS_THREAD_VARS:
        monkeypatch.delenv(var, raising=False)
    preset_blas_threads(["fedcat", "meta", "--threads", "3"])
    assert all(var not in os.environ for var in BLAS_THREAD_VARS)
    preset_blas_threads(["fedcat", "collinearity", "--pfile", "x", "-t", "3"])
    assert all(os.environ[var] == "3" for var in BLAS_THREAD_VARS)


def test_run_collinearity_check_reads_the_pfile_in_batches(tmp_path):
    """The file-level check, as tasks run it: IDs resolved through the pvar, candidates streamed in batches."""
    pg = pytest.importorskip("pgenlib")
    from fedcat.stats.collinearity import run_collinearity_check

    rng = np.random.default_rng(3)
    n_samples, ids = 60, [f"1:{100 * (i + 1)}:A:G" for i in range(6)]
    geno = rng.integers(0, 3, size=(len(ids), n_samples)).astype(np.int8)
    geno[4] = geno[0]  # a candidate identical to a selected variant
    geno[5, :3] = -9  # with some missing calls
    with pg.PgenWriter(str(tmp_path / "x.pgen").encode(), n_samples, len(ids), False) as writer:
        for row in geno:
            writer.append_biallelic(np.ascontiguousarray(row))
    (tmp_path / "x.pvar").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n" + "".join(f"1\t{v.split(':')[1]}\t{v}\tA\tG\n" for v in ids)
    )
    (tmp_path / "selected").write_text(f"{ids[0]}\n{ids[1]}\n")
    (tmp_path / "candidates").write_text("\n".join([ids[2], ids[3], "1:999:A:G", ids[4], ids[5]]) + "\n")

    out = tmp_path / "out.tsv"
    run_collinearity_check(
        str(tmp_path / "x"), str(tmp_path / "selected"), str(tmp_path / "candidates"), str(out), batch_size=2
    )

    got = {
        row[0]: (float(row[1]), float(row[2]))
        for row in (line.split("\t") for line in out.read_text().splitlines()[1:])
    }
    assert list(got) == [ids[2], ids[3], ids[4], ids[5]]  # the ID missing from the pvar is left out
    g = geno.astype(np.float64)
    g[g == -9] = np.nan
    slct, cand = g[:2].copy(), g[2:].copy()
    for block in (slct, cand):
        mean_centre(block)
    cand_rsq, max_rsq = candidate_collinearity(slct, compute_variances(slct), cand, compute_variances(cand))
    for j, vid in enumerate(ids[2:]):
        assert got[vid] == pytest.approx((cand_rsq[j], max_rsq[j]), rel=1e-9), vid
    assert got[ids[4]][1] >= 0.99  # flagged as collinear with the selected copy
