"""Collinearity of candidate variants with the selected (conditioned-on) variants.

For each candidate:

- `cand_rsq`: R² of the candidate regressed on the selected set (its VIF R²), how
  much of it the selected variants explain.
- `max_rsq`: the largest VIF R² in the set `selected ∪ {candidate}`, each member
  regressed on the others, which also catches a candidate making a selected variant
  collinear.

A candidate whose pairwise r² with a selected variant exceeds `rsq_threshold`, or
that makes the set singular, gets `COLLINEAR_SENTINEL` (2.0) for both. `pgenlib` is
only needed to read genotypes, by
[`run_collinearity_check`][fedcat.stats.collinearity.run_collinearity_check].
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from loguru import logger
from scipy.linalg import LinAlgError, cho_factor, cho_solve

# Missing genotype in pgenlib hardcalls
_MISSING_CODE = -9
# Reported instead of an R² for a candidate collinear with the selected set
COLLINEAR_SENTINEL = 2.0
# Residual variance, relative to the candidate's, at or below which it is in the span
# of the selected set; above rounding error, which differs between BLAS libraries
SPAN_TOLERANCE = 1e-10


# ---------------------------------------------------------------------------
# Genotype preprocessing
# ---------------------------------------------------------------------------
def mean_centre(g: np.ndarray) -> None:
    """Mean-centre each variant in place, setting missing genotypes to the mean (0 once centred).

    Args:
        g: `(n_variants, n_samples)` genotypes, missing as `nan`.
    """
    counts = np.sum(~np.isnan(g), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        row_mean = np.nansum(g, axis=1) / counts
    centred = g - row_mean[:, None]
    np.copyto(g, np.where(np.isnan(g), 0.0, centred))


def compute_variances(g: np.ndarray) -> np.ndarray:
    """Return each variant's variance, `sum(g**2) / n_samples`, of centred genotypes.

    Args:
        g: Centred `(n_variants, n_samples)` genotypes.

    Returns:
        `(n_variants,)` variances: the diagonal of the 1/n-scaled Gram matrix.
    """
    n = g.shape[1]
    return np.einsum("ij,ij->i", g, g) / n


def selected_gram(g_slct: np.ndarray) -> np.ndarray:
    """Return the 1/n-scaled `(k, k)` Gram matrix of centred `(k, n_samples)` genotypes."""
    n = g_slct.shape[1]
    return (g_slct @ g_slct.T) / n


def rsqs_from_gram(B: np.ndarray) -> np.ndarray:
    """Return each variable's VIF R², `1 - 1 / (B_ii · (B⁻¹)_ii)`, from a Gram matrix.

    Inverts `B` by Cholesky, or by pseudo-inverse when it is singular.

    Args:
        B: Gram matrix.

    Returns:
        VIF R² per variable.
    """
    try:
        c, low = cho_factor(B, lower=True)
        B_inv = cho_solve((c, low), np.identity(B.shape[0]))
    except LinAlgError:
        logger.warning("Cholesky decomposition failed; falling back to pseudo-inverse.")
        B_inv = np.linalg.pinv(B)
    return 1.0 - 1.0 / (np.diag(B) * np.diag(B_inv))


# ---------------------------------------------------------------------------
# Candidate collinearity
# ---------------------------------------------------------------------------
def candidate_collinearity(
    g_slct: np.ndarray,
    var_slct: np.ndarray,
    g_cand: np.ndarray,
    var_cand: np.ndarray,
    rsq_threshold: float = 0.99,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute each candidate's collinearity with the selected set.

    The selected Gram matrix `B` is factorised once, and every candidate is handled as
    a block-inverse update of it, in batched matrix operations. For a candidate with
    cross-covariances `b` to the selected set and variance `d`, the explained variance
    is `q = bᵀB⁻¹b` and the Schur complement `s = d - q`; then `cand_rsq = q / d`, and
    selected variant `i`'s VIF R² in the augmented set uses
    `(B_aug⁻¹)_ii = (B⁻¹)_ii + (B⁻¹b)_i² / s`.

    Args:
        g_slct: Centred selected genotypes, `(k, n_samples)`.
        var_slct: Selected variances, `(k,)`.
        g_cand: Centred candidate genotypes, `(m, n_samples)`.
        var_cand: Candidate variances, `(m,)`.
        rsq_threshold: Pairwise r² above which a candidate gets `COLLINEAR_SENTINEL`.

    Returns:
        `(cand_rsq, max_rsq)`, each of length `m`.
    """
    k, n = g_slct.shape
    m = g_cand.shape[0]

    # Cross-covariances, (k, m)
    Bmat = (g_slct @ g_cand.T) / n

    B = selected_gram(g_slct)
    try:
        cho = cho_factor(B, lower=True)
        diag_Binv = np.diag(cho_solve(cho, np.identity(k)))
        W = cho_solve(cho, Bmat)  # B⁻¹ @ Bmat, (k, m)
    except LinAlgError:
        logger.warning("Cholesky decomposition of selected Gram failed; using pseudo-inverse.")
        B_inv = np.linalg.pinv(B)
        diag_Binv = np.diag(B_inv)
        W = B_inv @ Bmat

    q = np.einsum("ij,ij->j", Bmat, W)  # explained variance per candidate
    s = var_cand - q  # residual variance

    cand_rsq = np.divide(q, var_cand, out=np.zeros_like(q), where=var_cand > 0)

    with np.errstate(divide="ignore", invalid="ignore"):
        inv_diag = diag_Binv[:, None] + (W * W) / s[None, :]  # (k, m)
        rsq_sel = 1.0 - 1.0 / (var_slct[:, None] * inv_diag)  # (k, m)
    max_sel = np.max(rsq_sel, axis=0) if k else np.zeros(m)
    max_rsq = np.maximum(cand_rsq, max_sel)

    # A residual variance within rounding error of zero means the candidate is in the
    # span of the selected set; a high pairwise r² catches the common case directly
    pairwise_rsq = (Bmat * Bmat) / (var_cand[None, :] * var_slct[:, None])  # (k, m)
    degenerate = (s <= SPAN_TOLERANCE * var_cand) | np.any(pairwise_rsq > rsq_threshold, axis=0)
    cand_rsq[degenerate] = COLLINEAR_SENTINEL
    max_rsq[degenerate] = COLLINEAR_SENTINEL
    return cand_rsq, max_rsq


# ---------------------------------------------------------------------------
# Reading genotypes
# ---------------------------------------------------------------------------
def _read_variants(pgen_path: Path, pvar_path: Path, indices: np.ndarray) -> np.ndarray:
    """Read hardcall genotypes of pvar rows as a `(n_variants, n_samples)` float array, missing as `nan`."""
    import pgenlib as pg

    with (
        pg.PvarReader(str(pvar_path).encode()) as pvar_reader,
        pg.PgenReader(str(pgen_path).encode(), pvar=pvar_reader) as pgen_reader,
    ):
        n_samples = pgen_reader.get_raw_sample_ct()
        out = np.empty((len(indices), n_samples), dtype=np.int8)
        pgen_reader.read_list(indices.astype(np.uint32), out, sample_maj=False)
    g = out.astype(np.float64)
    g[g == _MISSING_CODE] = np.nan
    return g


def _prepare_block(pgen_path: Path, pvar_path: Path, indices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Read and mean-centre genotypes of pvar rows, returning them with their variances."""
    g = _read_variants(pgen_path, pvar_path, indices)
    mean_centre(g)
    return g, compute_variances(g)


def run_collinearity_check(
    pfile: str,
    selected: str,
    candidate: str,
    out: str,
    batch_size: int = 1000,
    rsq_threshold: float = 0.99,
) -> None:
    """Write each candidate variant's collinearity with the selected set to a TSV.

    Variants not in the pvar are skipped.

    Args:
        pfile: PLINK2 file set prefix, with `.pgen` and `.pvar` files.
        selected: File of selected variant IDs, one per line.
        candidate: File of candidate variant IDs, one per line.
        out: Output TSV with columns `ID`, `cand_rsq`, `max_rsq`.
        batch_size: Candidates processed per batch.
        rsq_threshold: Pairwise r² above which a candidate gets `COLLINEAR_SENTINEL`.
    """
    start = time.time()

    selected_ids = _read_id_list(Path(selected))
    candidate_ids = _read_id_list(Path(candidate))

    pgen_path = Path(f"{pfile}.pgen")
    pvar_path = Path(f"{pfile}.pvar")

    logger.info("Reading pvar variant IDs...")
    id_to_index = _pvar_index(pvar_path)

    slct_idx = _resolve_indices(selected_ids, id_to_index)
    logger.info(f"Found {len(slct_idx)}/{len(selected_ids)} selected variants in pvar")
    g_slct, var_slct = _prepare_block(pgen_path, pvar_path, slct_idx)

    if len(slct_idx) > 1:
        rsqs_slct = rsqs_from_gram(selected_gram(g_slct))
        logger.info(f"Maximum R² among selected variants: {rsqs_slct.max():.4f}")

    present_cand_ids = [cid for cid in candidate_ids if cid in id_to_index]
    n_cand = len(present_cand_ids)
    logger.info(f"Found {n_cand}/{len(candidate_ids)} candidate variants in pvar")

    with open(out, "w") as fh:
        fh.write("ID\tcand_rsq\tmax_rsq\n")
        for i in range(0, n_cand, batch_size):
            batch_ids = present_cand_ids[i : i + batch_size]
            logger.info(f"Processing batch {i // batch_size + 1}/{(n_cand + batch_size - 1) // batch_size}...")
            cand_idx = np.array([id_to_index[c] for c in batch_ids], dtype=np.int64)
            g_cand, var_cand = _prepare_block(pgen_path, pvar_path, cand_idx)
            cand_rsq, max_rsq = candidate_collinearity(
                g_slct, var_slct, g_cand, var_cand, rsq_threshold=rsq_threshold
            )
            fh.writelines(f"{cid}\t{cand_rsq[j]}\t{max_rsq[j]}\n" for j, cid in enumerate(batch_ids))

    logger.success(f"Collinearity check complete in {time.time() - start:.2f}s → {out}")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _read_id_list(path: Path) -> list[str]:
    """Read non-empty lines of a file."""
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def _pvar_index(pvar_path: Path) -> dict[str, int]:
    """Return each variant ID's pvar row index."""
    import pgenlib as pg

    with pg.PvarReader(str(pvar_path).encode()) as reader:
        return {reader.get_variant_id(i).decode(): i for i in range(reader.get_variant_ct())}


def _resolve_indices(ids: list[str], id_to_index: dict[str, int]) -> np.ndarray:
    """Return the pvar row indices of the IDs found in the pvar."""
    return np.array([id_to_index[i] for i in ids if i in id_to_index], dtype=np.int64)
