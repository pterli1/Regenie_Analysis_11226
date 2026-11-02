"""Tests for fedcat.stats.meta_analysis: Z → -log10(P) and the IVW meta-analysis on both backends.

Expected values come from an independent IVW and scipy, not from the code under test.
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from loguru import logger
from scipy.stats import chi2, norm

from fedcat.__main__ import _meta_main
from fedcat.stats.meta_analysis import log10p_to_p, meta_analyse_regenie, z_to_log10p

META_OUTPUT_COLUMNS = [
    "CHROM", "POS", "ID", "A0", "A1", "A1FREQ", "N", "TEST", "BETA", "SE", "Z", "LOG10P",
    "Q", "Q_P", "I2", "N_STUDIES", "DIRECTIONS",
]  # fmt: skip
REGENIE_COLUMNS = [
    "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA",
]  # fmt: skip


def _scipy_log10p(z: float) -> float:
    """Reference -log10(two-sided P) of a Z-score, by scipy."""
    return -(norm.logsf(abs(z)) + math.log(2.0)) / math.log(10.0)


def _row(pos: int = 100, chrom: str = "1", **overrides) -> dict:
    """A regenie result row; the ID defaults to chrom:pos:A:T."""
    row = {
        "CHROM": chrom, "GENPOS": pos, "ID": f"{chrom}:{pos}:A:T", "ALLELE0": "A", "ALLELE1": "T",
        "A1FREQ": 0.1, "N": 1000, "TEST": "ADD", "BETA": 0.5, "SE": 0.1, "CHISQ": 25, "LOG10P": 5.0, "EXTRA": "NA",
    }  # fmt: skip
    row.update(overrides)
    return row


def _chrom_key(chrom: str) -> int:
    """Sort key of a chromosome: 1-22, X, Y, MT, then others."""
    return int(chrom) if chrom.isdigit() else {"X": 23, "Y": 24, "MT": 25}.get(chrom.upper(), 99)


def _make_regenie_file(
    tmp_dir: Path, label: str, rows: list[dict], *, compression: str = "gz", add_comments: bool = False
) -> Path:
    """Write *rows* as a regenie file, sorted by (CHROM, GENPOS) as the merge backend needs."""
    rows = sorted(rows, key=lambda r: (_chrom_key(str(r["CHROM"])), int(r["GENPOS"])))
    lines = [" ".join(REGENIE_COLUMNS)]
    for row in rows:
        if add_comments:
            lines.append("## synthetic comment before data row")
        lines.append(" ".join(str(row[c]) for c in REGENIE_COLUMNS))
    text = ("## synthetic comment before header\n" if add_comments else "") + "\n".join(lines) + "\n"
    match compression:
        case "gz":
            path = tmp_dir / f"{label}.regenie.gz"
            with gzip.open(path, "wt") as fh:
                fh.write(text)
        case "zst":
            zstd = pytest.importorskip("zstandard")
            path = tmp_dir / f"{label}.regenie.zst"
            path.write_bytes(zstd.ZstdCompressor(level=3).compress(text.encode()))
        case _:
            path = tmp_dir / f"{label}.regenie"
            path.write_text(text)
    return path


def _meta(tmp_path: Path, cohorts: dict[str, list[dict]], backend: str = "streaming", **kwargs) -> pl.DataFrame:
    """Meta-analyse cohorts' rows with a backend, returning the output."""
    files = {label: [_make_regenie_file(tmp_path, f"{label}_{backend}", rows)] for label, rows in cohorts.items()}
    out = tmp_path / f"meta_{backend}.tsv"
    meta_analyse_regenie(files, out, kwargs.pop("res_type", "single_variants"), backend=backend, **kwargs)
    return _read_meta_output(out)


def _read_meta_output(path: Path) -> pl.DataFrame:
    """Read meta-analysis output."""
    return pl.read_csv(
        path, separator="\t", comment_prefix="##", schema_overrides={"CHROM": pl.Utf8, "POS": pl.Int64}
    )


@pytest.fixture(params=["streaming", "merge"])
def backend(request):
    """Each meta-analysis backend."""
    return request.param


# ── Z → -log10(P) ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize("use_rust", [True, False])
def test_z_to_log10p_matches_scipy(use_rust):
    z = np.array([0.0, 0.5, 1.96, -1.96, 3.0, 5.45, -10.0, 20.0, 37.0])
    expected = [_scipy_log10p(v) for v in z]
    np.testing.assert_allclose(z_to_log10p(z, use_rust=use_rust), expected, rtol=1e-9, atol=1e-12)
    assert len(z_to_log10p(np.array([]), use_rust=use_rust)) == 0


def test_z_to_log10p_stays_finite_past_float_underflow():
    """Past |Z| ~ 37 the P-value underflows; the asymptotic expansion keeps LOG10P finite and increasing."""
    z = np.array([38.0, 40.0, 100.0, -100.0])
    out = z_to_log10p(z)
    # log10 of the leading term of the normal tail, phi(z)/z, which the expansion refines
    leading = [
        -(-0.5 * v * v - math.log(abs(v) * math.sqrt(2 * math.pi)) + math.log(2.0)) / math.log(10) for v in z
    ]
    np.testing.assert_allclose(out, leading, rtol=1e-3)
    assert np.all(np.diff(out[:3]) > 0) and out[2] == out[3]


def test_log10p_to_p():
    np.testing.assert_allclose(log10p_to_p(np.array([0.0, 1.0, -math.log10(0.05), 400.0])), [1.0, 0.1, 0.05, 0.0])


# ── IVW against an independent calculation ────────────────────────────────────


def _manual_ivw(betas, ses, a1freqs, ns) -> dict:
    """Reference IVW meta-analysis of one variant's per-cohort results."""
    betas, ses, a1freqs, ns = (np.asarray(x, dtype=np.float64) for x in (betas, ses, a1freqs, ns))
    w = 1.0 / ses**2
    beta = (w * betas).sum() / w.sum()
    se = 1.0 / math.sqrt(w.sum())
    q = float((w * (betas - beta) ** 2).sum())
    df = len(betas) - 1
    return {
        "BETA": beta,
        "SE": se,
        "Z": beta / se,
        "LOG10P": _scipy_log10p(beta / se),
        "A1FREQ": (a1freqs * ns).sum() / ns.sum(),
        "N": int(ns.sum()),
        "Q": q,
        "Q_P": chi2.sf(q, df),
        "I2": max(0.0, (q - df) / q),
        "N_STUDIES": len(betas),
    }


def test_three_cohort_ivw_matches_manual(tmp_path, backend):
    data = [(0.30, 0.10, 0.15, 1000), (0.25, 0.12, 0.18, 1500), (0.65, 0.08, 0.12, 2000)]
    cohorts = {f"c{i}": [_row(BETA=b, SE=s, A1FREQ=f, N=n)] for i, (b, s, f, n) in enumerate(data)}
    df = _meta(tmp_path, cohorts, backend)
    assert list(df.columns) == META_OUTPUT_COLUMNS
    assert df.height == 1
    got = df.row(0, named=True)
    expected = _manual_ivw(*zip(*data))
    assert expected["I2"] > 0  # heterogeneous enough to test I2 and Q_P away from their bounds
    for col, value in expected.items():
        assert got[col] == pytest.approx(value, rel=1e-3), col
    assert (got["CHROM"], got["POS"], got["A0"], got["A1"], got["TEST"]) == ("1", 100, "A", "T", "ADD")


def test_many_variants_match_manual_ivw(tmp_path, backend):
    """50 variants x 5 cohorts, each cohort missing some: every row matches the IVW of the cohorts it has."""
    rng = np.random.default_rng(123)
    per_variant: dict[str, list[tuple]] = {}
    cohorts: dict[str, list[dict]] = {}
    for ci in range(5):
        rows = []
        for vi in range(50):
            if rng.random() < 0.15:
                continue
            values = (round(rng.normal(0, 0.5), 6), round(rng.uniform(0.05, 0.3), 6))
            values += (round(rng.uniform(0.01, 0.49), 6), int(rng.integers(500, 5000)))
            row = _row(1000 + vi * 10, BETA=values[0], SE=values[1], A1FREQ=values[2], N=values[3])
            rows.append(row)
            per_variant.setdefault(row["ID"], []).append(values)
        cohorts[f"c{ci}"] = rows
    df = _meta(tmp_path, cohorts, backend)
    assert sorted(df["ID"].to_list()) == sorted(per_variant)
    for got in df.iter_rows(named=True):
        expected = _manual_ivw(*zip(*per_variant[got["ID"]]))
        for col in ("BETA", "SE", "A1FREQ", "N", "N_STUDIES"):
            assert got[col] == pytest.approx(expected[col], rel=1e-4), (got["ID"], col)


def test_single_study_passes_through(tmp_path, backend):
    """One study: its estimate unchanged, with no heterogeneity whatever float noise Q picks up."""
    got = _meta(tmp_path, {"c0": [_row(500, "5", BETA=0.75, SE=0.15, A1FREQ=0.2, N=3000)]}, backend).row(
        0, named=True
    )
    assert (got["BETA"], got["SE"], got["A1FREQ"], got["N"]) == (
        pytest.approx(0.75),
        pytest.approx(0.15),
        0.2,
        3000,
    )
    assert (got["N_STUDIES"], got["I2"], got["DIRECTIONS"]) == (1, 0.0, "+")
    # No Q p-value from one study
    assert got["Q_P"] in (None, "NA")


def test_backends_agree(tmp_path):
    rng = np.random.default_rng(42)
    cohorts = {
        f"c{ci}": [
            _row(
                100 + vi * 10,
                str(vi % 3 + 1),
                BETA=round(rng.normal(0, 0.5), 6),
                SE=round(rng.uniform(0.05, 0.3), 6),
            )
            for vi in range(20)
            if rng.random() < 0.85
        ]
        for ci in range(4)
    }
    df_s = _meta(tmp_path, cohorts, "streaming").sort("ID")
    df_m = _meta(tmp_path, cohorts, "merge").sort("ID")
    for col in ("ID", "N", "N_STUDIES", "DIRECTIONS"):
        assert df_s[col].to_list() == df_m[col].to_list(), col
    # merge writes fewer decimal places
    for col in ("BETA", "SE", "Z", "A1FREQ", "Q", "I2", "LOG10P"):
        np.testing.assert_allclose(df_s[col].to_numpy(), df_m[col].to_numpy(), rtol=5e-4, atol=1e-4, err_msg=col)


# ── Input handling ────────────────────────────────────────────────────────────


def test_empty_cohort_list_raises(tmp_path, backend):
    with pytest.raises(ValueError, match="empty"):
        meta_analyse_regenie({}, tmp_path / "nope.tsv", "single_variants", backend=backend)


@pytest.mark.parametrize(
    "column, sentinel",
    [("BETA", s) for s in ("NA", "nan", "-nan", "NaN", "inf", "-inf")]
    + [("SE", "-nan"), ("A1FREQ", "nan"), ("N", "NA"), ("N", "inf")],
)
def test_non_finite_values_drop_the_row(tmp_path, backend, column, sentinel):
    """regenie writes NA/nan/-nan for stats it could not compute; that cohort's row is dropped, not parsed."""
    got = _meta(tmp_path, {"dirty": [_row(**{column: sentinel})], "clean": [_row(BETA=0.4, SE=0.08)]}, backend)
    assert got.height == 1
    assert (got["N_STUDIES"][0], got["BETA"][0], got["SE"][0]) == (1, pytest.approx(0.4), pytest.approx(0.08))


def test_variant_non_finite_in_every_cohort_is_dropped(tmp_path, backend):
    assert _meta(tmp_path, {"c1": [_row(BETA="NA")], "c2": [_row(BETA="-nan")]}, backend).height == 0


def test_burden_ids_are_kept(tmp_path, backend):
    df = _meta(
        tmp_path,
        {"c1": [_row(0, ID="GENE1.mask1.0.001", ALLELE0="ref", ALLELE1="alt")]},
        backend,
        res_type="burden",
    )
    assert df["ID"].to_list() == ["GENE1.mask1.0.001"]


def test_output_in_natural_chromosome_order(tmp_path, backend):
    rows = [_row(100, c, ID=f"{c}:100") for c in ("2", "23", "1", "10")]
    assert _meta(tmp_path, {"c1": rows}, backend)["CHROM"].to_list() == ["1", "2", "10", "X"]


@pytest.mark.parametrize("compression", ["plain", "gz", "zst"])
def test_compressed_and_commented_inputs(tmp_path, backend, compression):
    files = {
        label: [
            _make_regenie_file(tmp_path, f"{label}_{backend}", [row], compression=compression, add_comments=True)
        ]
        for label, row in (("c1", _row()), ("c2", _row(BETA=0.4, SE=0.11)))
    }
    out = tmp_path / "meta.tsv"
    meta_analyse_regenie(files, out, "single_variants", backend=backend)
    df = _read_meta_output(out)
    assert (df.height, df["N_STUDIES"][0]) == (1, 2)


# ── Output annotations ────────────────────────────────────────────────────────


def test_position_column_can_be_renamed(tmp_path, backend):
    files = {}
    for label, beta in (("c1", 0.4), ("c2", 0.6)):
        path = _make_regenie_file(tmp_path, f"{label}_{backend}", [_row(100, BETA=beta), _row(250, BETA=beta)])
        renamed = tmp_path / f"{label}_{backend}_pos.regenie.gz"
        with gzip.open(path, "rt") as src, gzip.open(renamed, "wt") as dst:
            dst.write(src.read().replace("GENPOS", "BP", 1))
        files[label] = [renamed]
    out = tmp_path / "meta.tsv"
    meta_analyse_regenie(files, out, "single_variants", column_config={"genpos": "BP"}, backend=backend)
    df = _read_meta_output(out)
    assert df["POS"].to_list() == [100, 250]
    assert df["N_STUDIES"].to_list() == [2, 2]


def test_matching_duplicate_rows_count_once(tmp_path, backend):
    """A cohort's variant in two overlapping files, or twice in one, counts once."""
    files = {
        "c1": [
            _make_regenie_file(tmp_path, f"c1a_{backend}", [_row(100, BETA=0.4), _row(200, BETA=0.2)]),
            _make_regenie_file(tmp_path, f"c1b_{backend}", [_row(200, BETA=0.2), _row(300)]),
        ],
        "c2": [_make_regenie_file(tmp_path, f"c2_{backend}", [_row(100, BETA=0.6), _row(100, BETA=0.6)])],
    }
    out = tmp_path / "meta.tsv"
    meta_analyse_regenie(files, out, "single_variants", backend=backend)
    df = _read_meta_output(out)
    assert dict(zip(df["ID"], df["N_STUDIES"])) == {"1:100:A:T": 2, "1:200:A:T": 1, "1:300:A:T": 1}
    assert df.filter(pl.col("ID") == "1:100:A:T")["BETA"][0] == pytest.approx(0.5, rel=1e-4)


@pytest.mark.parametrize("change", [{"BETA": 0.3}, {"SE": 0.2}])
def test_conflicting_duplicate_rows_raise(tmp_path, backend, change):
    files = {
        "c1": [
            _make_regenie_file(tmp_path, f"c1a_{backend}", [_row(100)]),
            _make_regenie_file(tmp_path, f"c1b_{backend}", [_row(100, **change)]),
        ]
    }
    with pytest.raises(ValueError, match="Conflicting results for variant '1:100:A:T' .* in cohort 'c1'"):
        meta_analyse_regenie(files, tmp_path / "meta.tsv", "single_variants", backend=backend)


def test_merge_warns_of_unsorted_input_with_both_positions(tmp_path, capfd):
    path = tmp_path / "unsorted.regenie"
    rows = [_row(300), _row(100)]
    path.write_text(
        " ".join(REGENIE_COLUMNS)
        + "\n"
        + "".join(" ".join(str(r[c]) for c in REGENIE_COLUMNS) + "\n" for r in rows)
    )
    meta_analyse_regenie({"c1": [path]}, tmp_path / "meta.tsv", "single_variants", backend="merge")
    assert "not sorted at pos 1:100 (prev 1:300)" in capfd.readouterr().err


def test_directions_give_each_cohort_sign_or_absence(tmp_path, backend):
    cohorts = {
        "c1": [_row(100, BETA=0.3), _row(200, BETA=0.1)],
        "c2": [_row(100, BETA=-0.2)],
        "c3": [_row(200, BETA=-0.1)],
    }
    df = _meta(tmp_path, cohorts, backend)
    assert dict(zip(df["ID"], df["DIRECTIONS"])) == {"1:100:A:T": "+-?", "1:200:A:T": "+?-"}


def test_output_names_cohorts_in_order(tmp_path, backend):
    files = {
        name: [_make_regenie_file(tmp_path, f"{name}_{backend}", [_row()])]
        for name in ("zebra", "apple", "mango")
    }
    out = tmp_path / "meta.tsv"
    meta_analyse_regenie(files, out, "single_variants", backend=backend)
    assert out.read_text().splitlines()[0] == "## Cohorts: zebra, apple, mango"


# ── Per-stratum AAF bins (single-variant input to ACAT) ───────────────────────


class TestAafBins:
    """A cohort contributes to bin t only where its own A1FREQ <= t."""

    @staticmethod
    def _cohorts() -> dict[str, list[dict]]:
        """Two cohorts where one variant is rare in only one, and another too common for any bin."""
        return {
            "A": [_row(100, A1FREQ=0.0005, N=10000, SE=0.2), _row(200, A1FREQ=0.02, N=10000, BETA=0.1, SE=0.05)],
            "B": [_row(100, A1FREQ=0.005, N=5000, BETA=0.3, SE=0.1)],
        }

    def test_bins_use_only_strata_rare_for_the_bin(self, tmp_path, backend):
        df = _meta(tmp_path, self._cohorts(), backend, aaf_bins=["0.001", "0.01"])
        assert list(df.columns) == META_OUTPUT_COLUMNS + ["AAF_BIN"]
        assert set(df["ID"].to_list()) == {"1:100:A:T"}
        rows = {str(r["AAF_BIN"]): r for r in df.iter_rows(named=True)}
        # bin 0.001: stratum A alone, so its estimate passes through unchanged
        assert (rows["0.001"]["N_STUDIES"], rows["0.001"]["DIRECTIONS"]) == (1, "+?")
        assert rows["0.001"]["BETA"] == pytest.approx(0.5, abs=1e-6)
        assert rows["0.001"]["A1FREQ"] == pytest.approx(0.0005, abs=1e-6)
        # bin 0.01: both strata, IVW-combined
        w_a, w_b = 1 / 0.2**2, 1 / 0.1**2
        assert rows["0.01"]["N_STUDIES"] == 2
        assert rows["0.01"]["BETA"] == pytest.approx((w_a * 0.5 + w_b * 0.3) / (w_a + w_b), abs=1e-6)
        assert rows["0.01"]["A1FREQ"] == pytest.approx((0.0005 * 10000 + 0.005 * 5000) / 15000, abs=1e-6)

    def test_without_bins_every_variant_is_kept(self, tmp_path, backend):
        df = _meta(tmp_path, self._cohorts(), backend)
        assert list(df.columns) == META_OUTPUT_COLUMNS
        assert df["ID"].to_list() == ["1:100:A:T", "1:200:A:T"]

    def test_no_rare_rows_writes_header_only(self, tmp_path, backend):
        df = _meta(tmp_path, self._cohorts(), backend, aaf_bins=["0.0001"])
        assert df.height == 0 and "AAF_BIN" in df.columns

    def test_invalid_bin_rejected(self, tmp_path, backend):
        with pytest.raises(ValueError, match="aaf bin"):
            _meta(tmp_path, self._cohorts(), backend, aaf_bins=["0"])


# ── CLI and warnings ──────────────────────────────────────────────────────────


def test_meta_cli_with_no_listed_files_writes_nothing(tmp_path):
    """A cohort file listing no files (e.g. no mac10 results) leaves no output for the workflow to find."""
    cohort_file = tmp_path / "cohorts.txt"
    cohort_file.write_text("# eur_aou\n# afr_aou\n")
    out = tmp_path / "meta.tsv.zst"
    _meta_main(argparse.Namespace(cohort_file=str(cohort_file), output_file=str(out), res_type="mac10"))
    assert not out.exists()


@pytest.mark.parametrize("compression", ["gz", "zst"])
@pytest.mark.parametrize(
    "ids, warned", [(("chr1:100:A:G", "1:100:A:G"), True), (("1:100:A:G", "1:100:A:G"), False)]
)
def test_meta_warns_when_cohorts_mix_chr_prefixed_ids(tmp_path, ids, warned, compression):
    """The meta matches IDs exactly, so cohorts disagreeing on 'chr' would never combine."""
    files = {
        label: [_make_regenie_file(tmp_path, label, [_row(ID=vid)], compression=compression)]
        for label, vid in zip(("aou", "ukb"), ids)
    }
    messages: list[str] = []
    sink = logger.add(lambda m: messages.append(m.record["message"]), level="WARNING")
    try:
        meta_analyse_regenie(files, tmp_path / "meta.tsv", "single_variants")
    finally:
        logger.remove(sink)
    hits = [m for m in messages if "mix 'chr'-prefixed and bare" in m]
    assert bool(hits) is warned
    if warned:
        assert "'chr' in aou; bare in ukb" in hits[0]
