"""Tests for ACAT: the Rust kernel against independent references, and compute_acat end to end.

The kernel is checked against an exact P-space ACAT where P is representable, and against
its analytic tail, P ~ 1 / sum(w / p), where it is not.
"""

from __future__ import annotations

import gzip
import math
from collections import defaultdict as _ddict
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from fedcat.stats.acat import compute_acat
from fedcat_rs import acat_log10p_rs, compute_acat_streaming_rs, convert_tsv_to_parquet_rs

ACAT_COLUMNS = [
    "CHR", "START", "END", "ID", "UNIT", "MASK", "THRESHOLD", "N_AGG_COMBINED", "N_SINGLEVAR_COMBINED", "LOG10P_ACAT",
]  # fmt: skip


def _ref_acat_log10p(log10p, weights=None) -> float:
    """Exact P-space ACAT: T = sum w tan((0.5 - p) pi) = sum w / tan(p pi), P = atan(1/T) / pi for T > 0."""
    p = np.minimum(10.0 ** -np.asarray(log10p, dtype=np.float64), 0.999)
    w = np.full(len(p), 1.0 / len(p)) if weights is None else np.asarray(weights) / np.sum(weights)
    t = float(np.sum(w / np.tan(p * np.pi)))
    return -math.log10(math.atan(1.0 / t) / math.pi if t > 0 else 0.5 - math.atan(t) / math.pi)


def _tail_acat_log10p(log10p, weights=None) -> float:
    """ACAT where P underflows: T ~ sum w / (pi p), so -log10 P ~ log10(sum w 10^v)."""
    v = np.asarray(log10p, dtype=np.float64)
    w = np.full(len(v), 1.0 / len(v)) if weights is None else np.asarray(weights) / np.sum(weights)
    terms = v + np.log10(w)
    top = terms.max()
    return float(top + math.log10(np.sum(10.0 ** (terms - top))))


def _weight_standard(af: float) -> float:
    """Reference ACAT weight: maf × (1 - maf)^49."""
    maf = min(af, 1 - af)
    return maf * (1 - maf) ** 49


# ── Kernel ────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "log10p, weights",
    [
        ([2.0, 3.0, 1.5], None),
        ([2.0, 2.0, 2.0], None),
        ([0.5, 1.0, 1.5, 2.0, 3.0], None),
        ([10.0, 20.0, 30.0], None),
        ([100.0, 200.0, 300.0], None),
        ([14.9, 15.1, 3.0], None),  # either side of the asymptotic switch at 15
        ([0.01, 0.1, 0.5], None),
        ([0.05, 0.1, 20.0], None),  # p > 0.5 terms pull T down
        ([0.01, 0.02], None),  # every p > 0.5, T < 0
        ([2.0, 3.0, 1.5], [1.0, 2.0, 0.5]),
        ([5.0, 10.0, 15.0], [0.1, 0.3, 0.6]),
    ],
)
def test_kernel_matches_exact_acat(log10p, weights):
    w = None if weights is None else np.array(weights)
    assert acat_log10p_rs(np.array(log10p), weights=w) == pytest.approx(
        _ref_acat_log10p(log10p, weights), rel=1e-9
    )


@pytest.mark.parametrize(
    "log10p, weights",
    [
        ([500.0, 1.0, 2.0], None),
        ([1000.0, 500.0, 100.0], None),
        ([3000.0, 0.5], None),
        ([308.0, 308.0, 308.0], None),
        ([400.0, 300.0, 200.0, 100.0, 50.0], None),
        ([1000.0, 0.5], [0.9, 0.1]),
        ([3000.0, 500.0, 50.0], [2.0, 3.0, 5.0]),
    ],
)
def test_kernel_matches_the_analytic_tail_past_underflow(log10p, weights):
    w = None if weights is None else np.array(weights)
    assert acat_log10p_rs(np.array(log10p), weights=w) == pytest.approx(
        _tail_acat_log10p(log10p, weights), rel=1e-9
    )


@pytest.mark.parametrize("v", [0.5, 3.0, 12.0, 14.9, 15.1, 300.0, 3000.0])
def test_kernel_single_p_value_combines_to_itself(v):
    assert acat_log10p_rs(np.array([v])) == pytest.approx(v, rel=1e-9)


def test_kernel_weights_are_normalised_and_order_free():
    v = np.array([3.0, 7.0, 1.0])
    w = np.array([1.0, 2.0, 3.0])
    assert acat_log10p_rs(v, weights=w) == pytest.approx(acat_log10p_rs(v, weights=10 * w), rel=1e-12)
    assert acat_log10p_rs(v[::-1].copy(), weights=w[::-1].copy()) == pytest.approx(
        acat_log10p_rs(v, weights=w), rel=1e-12
    )
    assert math.isnan(acat_log10p_rs(np.array([])))


# ── compute_acat on generated inputs ─────────────────────────────────────────

# Reproducible RNG
_RNG = np.random.default_rng(12345)

# Masks and annotations
MASK_ANNOTS = {
    "mask1": ["pLoF"],
    "mask2": ["pLoF", "missense"],
    "mask3": ["pLoF", "missense", "synonymous"],
}
ALL_ANNOTS = sorted({a for annots in MASK_ANNOTS.values() for a in annots})

# 3 genes across 2 chromosomes
GENES = [
    ("GENE_A", "1"),
    ("GENE_B", "1"),
    ("GENE_C", "2"),
]

# Generate variant IDs — 10 variants per gene, positions 1000–10000
VARIANTS: dict[str, list[tuple[str, int]]] = {}
for gene, chrom in GENES:
    positions = sorted(_RNG.integers(1000, 10000, size=10))
    bases = ["A", "C", "G", "T"]
    vids = []
    for pos in positions:
        ref, alt = _RNG.choice(bases, size=2, replace=False)
        vid = f"{chrom}:{pos}:{ref}:{alt}"
        vids.append((vid, int(pos)))
    VARIANTS[gene] = vids


def _write_singles_tsv(path: Path) -> None:
    """Write random single-variant results for every test variant."""
    rows = []
    for gene, chrom in GENES:
        for vid, pos in VARIANTS[gene]:
            af = float(_RNG.uniform(0.0001, 0.0005))
            log10p = float(_RNG.uniform(0.5, 50.0))
            rows.append({"ID": vid, "LOG10P": round(log10p, 4), "A1FREQ": round(af, 6)})
    pl.DataFrame(rows).write_csv(path, separator="\t")


def _write_burden_tsv(path: Path) -> None:
    """Write random burden results for every test gene and mask."""
    rows = []
    for gene, chrom in GENES:
        for mask_name in MASK_ANNOTS:
            burden_id = f"{gene}.{mask_name}.0.001"
            af = float(_RNG.uniform(0.00005, 0.0005))
            log10p = float(_RNG.uniform(0.5, 50.0))
            rows.append(
                {
                    "ID": burden_id,
                    "LOG10P": round(log10p, 4),
                    "A1FREQ": round(af, 6),
                    "TEST": "ADD",
                    "CHROM": chrom,
                }
            )
    pl.DataFrame(rows).write_csv(path, separator="\t")


def _with_comments_and_compression(src: Path, dst: Path, compression: str) -> Path:
    """Copy `src` to `dst` with a `##` comment line, compressed as `"plain"`, `"gz"` or `"zst"`."""
    text = src.read_text()
    payload = f"## comment: should be skipped\n{text}"
    if compression == "plain":
        dst.write_text(payload)
    elif compression == "gz":
        with gzip.open(dst, "wt") as fh:
            fh.write(payload)
    elif compression == "zst":
        zstd = pytest.importorskip("zstandard")
        cctx = zstd.ZstdCompressor(level=3)
        with open(dst, "wb") as fh:
            fh.write(cctx.compress(payload.encode("utf-8")))
    else:
        raise ValueError(f"Unsupported compression: {compression}")
    return dst


def _build_annotations() -> tuple[list[str], list[str], list[str]]:
    """Return random annotations of the test variants as parallel (variant, unit, annotation) lists."""
    vids, units, annots = [], [], []
    for gene, chrom in GENES:
        for vid, _pos in VARIANTS[gene]:
            n_annots = _RNG.integers(1, 4)
            chosen = _RNG.choice(ALL_ANNOTS, size=n_annots, replace=False)
            for annot in chosen:
                vids.append(vid)
                units.append(gene)
                annots.append(annot)
    return vids, units, annots


def _write_tsv_annotations(vids, units, annots, out_dir: Path) -> list[Path]:
    """Write annotations as REGENIE annotation files, one per chromosome."""
    out_dir.mkdir(parents=True, exist_ok=True)
    chrs = [v.split(":")[0] for v in vids]
    paths = []
    for c in sorted(set(chrs)):
        path = out_dir / f"annotations_chr{c}.tsv"
        with open(path, "w") as fh:
            for vid, unit, annot, ch in zip(vids, units, annots, chrs):
                if ch == c:
                    fh.write(f"{vid}\t{unit}\t{annot}\n")
        paths.append(path)
    return sorted(paths)


def _write_optimised_annotations(
    vids,
    units,
    annots,
    out_dir: Path,
    sort_order: str = "unit_start",
) -> list[Path]:
    """Write optimised-schema annotation directory."""
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pl.DataFrame({"VARIANT_ID": vids, "UNIT": units, "ANNOT": annots})
    df = df.with_columns(
        pl.col("VARIANT_ID").str.split(":").list.get(0).alias("CHR"),
        pl.col("VARIANT_ID").str.split(":").list.get(1).cast(pl.Int64).alias("POS"),
    )

    unique_vars = sorted(df["VARIANT_ID"].unique().to_list())
    unique_units = sorted(df["UNIT"].unique().to_list())
    unique_annots = sorted(df["ANNOT"].unique().to_list())

    var_to_code = {v: i for i, v in enumerate(unique_vars)}
    unit_to_code = {u: i for i, u in enumerate(unique_units)}
    annot_to_code = {a: i for i, a in enumerate(unique_annots)}

    pl.DataFrame({"VARIANT_CODE": list(range(len(unique_vars))), "VARIANT_ID": unique_vars}).write_parquet(
        out_dir / "dim_variants.parquet"
    )
    pl.DataFrame({"UNIT_CODE": list(range(len(unique_units))), "UNIT": unique_units}).write_parquet(
        out_dir / "dim_units.parquet"
    )
    pl.DataFrame({"ANNOT_CODE": list(range(len(unique_annots))), "ANNOT": unique_annots}).write_parquet(
        out_dir / "dim_annots.parquet"
    )

    df = df.with_columns(
        pl.col("VARIANT_ID").replace(var_to_code).cast(pl.Int32).alias("VARIANT_CODE"),
        pl.col("UNIT").replace(unit_to_code).cast(pl.Int32).alias("UNIT_CODE"),
        pl.col("ANNOT").replace(annot_to_code).cast(pl.Int32).alias("ANNOT_CODE"),
        pl.col("CHR")
        .replace({"X": "23", "23": "23", "Y": "24", "MT": "25"})
        .cast(pl.Int32, strict=False)
        .fill_null(26)
        .alias("CHR_RANK"),
    )

    if sort_order == "position":
        df = df.sort(["CHR_RANK", "POS", "UNIT_CODE"])
    elif sort_order == "unit_start":
        df = df.with_columns(pl.col("POS").min().over(["UNIT_CODE", "CHR_RANK"]).alias("UNIT_MIN_POS")).sort(
            ["CHR_RANK", "UNIT_MIN_POS", "UNIT_CODE", "POS"]
        )
    else:
        df = df.sort(["UNIT_CODE", "CHR_RANK", "POS"])

    df.select(["VARIANT_CODE", "UNIT_CODE", "ANNOT_CODE", "CHR", "POS"]).write_parquet(
        out_dir / "annotations.parquet"
    )

    # unit_mask_bounds.parquet
    annot_to_masks: dict[str, list[str]] = _ddict(list)
    for mask_name, mask_annots_list in MASK_ANNOTS.items():
        for a in mask_annots_list:
            annot_to_masks[a].append(mask_name)
    umb_rows: list[dict] = []
    for (unit, chr_val), grp in df.group_by(["UNIT", "CHR"]):
        annots_in_group = grp["ANNOT"].to_list()
        positions_in_group = grp["POS"].to_list()
        mask_positions: dict[str, list[int]] = _ddict(list)
        for a, p in zip(annots_in_group, positions_in_group):
            for m in annot_to_masks.get(a, []):
                mask_positions[m].append(p)
        for m_name, positions in mask_positions.items():
            umb_rows.append(
                {
                    "UNIT": unit,
                    "CHR": chr_val,
                    "MASK": m_name,
                    "MIN_POS": min(positions),
                    "MAX_POS": max(positions),
                }
            )
        umb_rows.append(
            {
                "UNIT": unit,
                "CHR": chr_val,
                "MASK": "__UNIT__",
                "MIN_POS": min(positions_in_group),
                "MAX_POS": max(positions_in_group),
            }
        )
    if umb_rows:
        pl.DataFrame(umb_rows).write_parquet(out_dir / "unit_mask_bounds.parquet")

    # mask_annotations.parquet — must match schema written by convert_tsv_to_parquet_rs:
    # columns MASK (str), ANNOT_CODE (int64), ANNOT (str), one row per mask×annotation pair
    ma_rows = [
        {"MASK": mask_name, "ANNOT_CODE": annot_to_code[a], "ANNOT": a}
        for mask_name, mask_annots_list in MASK_ANNOTS.items()
        for a in mask_annots_list
        if a in annot_to_code
    ]
    pl.DataFrame(ma_rows, schema={"MASK": pl.Utf8, "ANNOT_CODE": pl.Int64, "ANNOT": pl.Utf8}).write_parquet(
        out_dir / "mask_annotations.parquet"
    )

    return [out_dir]


def _write_masks_file(path: Path) -> None:
    """Write the test masks as a REGENIE mask file."""
    with open(path, "w") as fh:
        fh.writelines(f"{mask_name}\t{','.join(annot_list)}\n" for mask_name, annot_list in MASK_ANNOTS.items())


@pytest.fixture(scope="module")
def pipeline_data(tmp_path_factory):
    """Shared ACAT inputs: results, masks, and annotations as TSVs and as a parquet schema."""
    base = tmp_path_factory.mktemp("acat_pipeline")
    _write_singles_tsv(base / "singles.tsv")
    _write_burden_tsv(base / "burden.tsv")
    _write_masks_file(base / "masks.tsv")
    vids, units, annots = _build_annotations()
    return {
        "singles": base / "singles.tsv",
        "burden": base / "burden.tsv",
        "masks": base / "masks.tsv",
        "tsv_annots": _write_tsv_annotations(vids, units, annots, base / "annots_tsv"),
        "schema": _write_optimised_annotations(vids, units, annots, base / "annots_optimised")[0],
        "base": base,
    }


def _run(data, out_name: str, *, singles=None, burden=None, tsv_annots=False, masks=None, annots=None, **kwargs):
    """compute_acat on the pipeline data, with the Python-built schema or the TSV annotations."""
    out = data["base"] / out_name
    annotations = (
        {"masks_tsvs": [masks or data["masks"]], "annotations_tsvs": annots or data["tsv_annots"]}
        if tsv_annots or annots
        else {"schema_dir": data["schema"]}
    )
    compute_acat(
        single_variants_file=singles or data["singles"],
        burden_file=burden or data["burden"],
        output_file=out,
        aggregate_freqs=["0.001"],
        n_threads=1,
        **annotations,
        **kwargs,
    )
    return pl.read_csv(str(out), separator="\t", null_values=["NA"]).sort("ID")


@pytest.fixture(scope="module")
def baseline(pipeline_data) -> pl.DataFrame:
    """ACAT output of the default settings, to compare other settings against."""
    return _run(pipeline_data, "baseline.tsv")


def _assert_same_results(df: pl.DataFrame, expected: pl.DataFrame, rel: float = 1e-9) -> None:
    """Assert two ACAT outputs match: IDs, counts and positions exactly, -log10(P) to `rel`."""
    assert df["ID"].to_list() == expected["ID"].to_list()
    for col in ("N_AGG_COMBINED", "N_SINGLEVAR_COMBINED", "START", "END"):
        assert df[col].to_list() == expected[col].to_list(), col
    np.testing.assert_allclose(df["LOG10P_ACAT"].to_numpy(), expected["LOG10P_ACAT"].to_numpy(), rtol=rel)


def test_baseline_output(baseline):
    assert baseline.columns == ACAT_COLUMNS
    assert baseline.height == len(GENES) * len(MASK_ANNOTS)
    assert (baseline["LOG10P_ACAT"] >= 0).all() and baseline["LOG10P_ACAT"].is_finite().all()


def test_tsv_annotations_match_the_prebuilt_schema(pipeline_data, baseline):
    """The TSV path converts with the Rust converter; the baseline uses a schema built here in Python."""
    _assert_same_results(_run(pipeline_data, "tsv.tsv", tsv_annots=True), baseline)


def test_converter_writes_the_schema(pipeline_data):
    singles = pl.read_csv(pipeline_data["singles"], separator="\t")
    valid_ids = singles.filter(pl.col("A1FREQ") <= 0.001)["ID"].to_list()
    out_dir = pipeline_data["base"] / "converted"
    info = convert_tsv_to_parquet_rs(
        [str(p) for p in pipeline_data["tsv_annots"]], [str(pipeline_data["masks"])], str(out_dir), valid_ids
    )
    assert (info["n_variants"], info["n_units"], info["n_annots"]) == (
        len(valid_ids),
        len(GENES),
        len(ALL_ANNOTS),
    )
    for name in (
        "dim_variants",
        "dim_units",
        "dim_annots",
        "annotations",
        "unit_mask_bounds",
        "mask_annotations",
    ):
        assert (out_dir / f"{name}.parquet").exists(), name


def test_sweep_and_arrow_backends_agree(pipeline_data, baseline):
    _assert_same_results(_run(pipeline_data, "arrow.tsv", backend="arrow"), baseline, rel=1e-12)


@pytest.mark.parametrize("compression, suffix", [("plain", ".tsv"), ("gz", ".tsv.gz"), ("zst", ".tsv.zst")])
@pytest.mark.parametrize("tsv_annots", [False, True], ids=["schema", "tsv_annots"])
def test_commented_and_compressed_results(pipeline_data, baseline, compression, suffix, tsv_annots):
    """A leading ## line and compression in the singles and burden files, on both annotation paths."""
    base = pipeline_data["base"]
    singles = _with_comments_and_compression(
        pipeline_data["singles"], base / f"s_{tsv_annots}{suffix}", compression
    )
    burden = _with_comments_and_compression(
        pipeline_data["burden"], base / f"b_{tsv_annots}{suffix}", compression
    )
    df = _run(
        pipeline_data, f"c_{compression}_{tsv_annots}.tsv", singles=singles, burden=burden, tsv_annots=tsv_annots
    )
    _assert_same_results(df, baseline)


@pytest.mark.parametrize("compression, suffix", [("plain", ".tsv"), ("gz", ".tsv.gz"), ("zst", ".tsv.zst")])
def test_commented_and_compressed_annotations(pipeline_data, baseline, compression, suffix):
    base = pipeline_data["base"]
    annots = [
        _with_comments_and_compression(src, base / f"a{i}_{compression}{suffix}", compression)
        for i, src in enumerate(pipeline_data["tsv_annots"])
    ]
    masks = _with_comments_and_compression(pipeline_data["masks"], base / f"m_{compression}{suffix}", compression)
    _assert_same_results(_run(pipeline_data, f"ca_{compression}.tsv", masks=masks, annots=annots), baseline)


# ── compute_acat on hand-built inputs ─────────────────────────────────────────


@pytest.mark.parametrize("scheme", ["standard", "equal"])
def test_aggregate_combines_its_mask_variants_and_burden(tmp_path, scheme):
    """Each aggregate is the ACAT of the bin's single variants in its mask plus its own burden."""
    singles = [  # (id, A1FREQ, LOG10P, annotation)
        ("1:100:A:C", 0.0004, 3.0, "pLoF"),
        ("1:200:A:G", 0.0009, 5.5, "missense"),
        ("1:300:A:T", 0.0002, 1.2, "pLoF"),
        ("1:400:C:T", 0.004, 9.0, "pLoF"),  # above the bin
    ]
    burdens = {"lof": (0.0006, 2.0), "all": (0.0015, 4.0)}  # mask: (A1FREQ, LOG10P)
    pl.DataFrame(
        {"ID": [s[0] for s in singles], "A1FREQ": [s[1] for s in singles], "LOG10P": [s[2] for s in singles]}
    ).write_csv(tmp_path / "singles.tsv", separator="\t")
    pl.DataFrame(
        [
            {"ID": f"G.{m}.0.001", "A1FREQ": af, "LOG10P": lp, "TEST": "ADD", "CHROM": "1"}
            for m, (af, lp) in burdens.items()
        ]
    ).write_csv(tmp_path / "burden.tsv", separator="\t")
    (tmp_path / "masks.tsv").write_text("lof\tpLoF\nall\tpLoF,missense\n")
    (tmp_path / "annots.tsv").write_text("".join(f"{s[0]}\tG\t{s[3]}\n" for s in singles))
    compute_acat(
        single_variants_file=tmp_path / "singles.tsv",
        burden_file=tmp_path / "burden.tsv",
        masks_tsvs=[tmp_path / "masks.tsv"],
        annotations_tsvs=[tmp_path / "annots.tsv"],
        output_file=tmp_path / "acat.tsv",
        weighting_scheme=scheme,
        n_threads=1,
    )
    rows = {r["MASK"]: r for r in pl.read_csv(tmp_path / "acat.tsv", separator="\t").iter_rows(named=True)}
    weight = _weight_standard if scheme == "standard" else (lambda af: 1.0)
    for mask, annots in {"lof": {"pLoF"}, "all": {"pLoF", "missense"}}.items():
        parts = [(af, lp) for _, af, lp, a in singles if a in annots and af <= 0.001] + [burdens[mask]]
        expected = _ref_acat_log10p([lp for _, lp in parts], [weight(af) for af, _ in parts])
        assert rows[mask]["LOG10P_ACAT"] == pytest.approx(expected, rel=1e-6), mask
        assert (rows[mask]["N_SINGLEVAR_COMBINED"], rows[mask]["N_AGG_COMBINED"]) == (len(parts) - 1, 1)
        assert (rows[mask]["START"], rows[mask]["END"]) == (100, 400)  # full bounds include the common variant


class TestMultipleThresholds:
    """ACAT with more than one aggregate-frequency bin (regenie --aaf-bins).

    Bins are cumulative: a variant with AAF <= t belongs to every bin whose
    bound is >= its AAF. Each threshold produces its own set of aggregates.
    """

    def _write_inputs(self, base: Path):
        """One gene, one mask (pLoF); two variants in the 0.001 bin and two more
        that only enter the 0.01 bin. Burden entries exist at both thresholds."""
        base.mkdir(parents=True, exist_ok=True)
        gene, chrom = "GENE_X", "1"
        # (vid, af, log10p)
        variants = [
            (f"{chrom}:1000:A:C", 0.0004, 3.0),  # 0.001 bin
            (f"{chrom}:2000:A:G", 0.0009, 4.0),  # 0.001 bin
            (f"{chrom}:3000:A:T", 0.0050, 5.0),  # 0.01 bin only
            (f"{chrom}:4000:C:G", 0.0080, 6.0),  # 0.01 bin only
        ]
        singles = base / "singles.tsv"
        pl.DataFrame([{"ID": v, "LOG10P": lp, "A1FREQ": af} for v, af, lp in variants]).write_csv(
            singles, separator="\t"
        )

        burden = base / "burden.tsv"
        pl.DataFrame(
            [
                {"ID": f"{gene}.mask1.{t}", "LOG10P": lp, "A1FREQ": 0.0003, "TEST": "ADD", "CHROM": chrom}
                for t, lp in [("0.001", 2.0), ("0.01", 2.5)]
            ]
        ).write_csv(burden, separator="\t")

        masks = base / "masks.tsv"
        with open(masks, "w") as fh:
            fh.write("mask1\tpLoF\n")
        annots = base / "annots_chr1.tsv"
        with open(annots, "w") as fh:
            fh.writelines(f"{v}\t{gene}\tpLoF\n" for v, _af, _lp in variants)
        return singles, burden, masks, annots

    def test_two_thresholds_produce_nested_bins(self, tmp_path):
        singles, burden, masks, annots = self._write_inputs(tmp_path)
        out = tmp_path / "multi.tsv"
        compute_acat(
            single_variants_file=singles,
            burden_file=burden,
            masks_tsvs=[masks],
            annotations_tsvs=[annots],
            output_file=out,
            aggregate_freqs=["0.001", "0.01"],
            log10p_col="LOG10P",
            weighting_scheme="standard",
            n_threads=1,
        )
        df = pl.read_csv(str(out), separator="\t", null_values=["NA"], schema_overrides={"THRESHOLD": pl.Utf8})

        # One aggregate row per (unit, mask, threshold)
        ids = set(df["ID"].to_list())
        assert "GENE_X.mask1.0.001" in ids
        assert "GENE_X.mask1.0.01" in ids

        by_thresh = {r["THRESHOLD"]: r for r in df.iter_rows(named=True)}
        # Both thresholds present, each with a THRESHOLD column matching its label
        assert set(by_thresh) == {"0.001", "0.01"}
        # 0.001 bin: only the two low-AF variants; 0.01 bin: all four (cumulative)
        assert by_thresh["0.001"]["N_SINGLEVAR_COMBINED"] == 2
        assert by_thresh["0.01"]["N_SINGLEVAR_COMBINED"] == 4
        # Burden component present at each threshold (its own .THRESHOLD suffix)
        assert by_thresh["0.001"]["N_AGG_COMBINED"] == 1
        assert by_thresh["0.01"]["N_AGG_COMBINED"] == 1

    def test_multi_threshold_first_bin_matches_single_run(self, tmp_path):
        """The 0.001 aggregate must be identical whether run alone or alongside 0.01."""
        singles, burden, masks, annots = self._write_inputs(tmp_path)

        def _run(freqs, name):
            out = tmp_path / name
            compute_acat(
                single_variants_file=singles,
                burden_file=burden,
                masks_tsvs=[masks],
                annotations_tsvs=[annots],
                output_file=out,
                aggregate_freqs=freqs,
                log10p_col="LOG10P",
                weighting_scheme="standard",
                n_threads=1,
            )
            return pl.read_csv(
                str(out), separator="\t", null_values=["NA"], schema_overrides={"THRESHOLD": pl.Utf8}
            )

        single = _run(["0.001"], "single.tsv")
        multi = _run(["0.001", "0.01"], "multi2.tsv")

        s_row = single.filter(pl.col("THRESHOLD") == "0.001").row(0, named=True)
        m_row = multi.filter(pl.col("THRESHOLD") == "0.001").row(0, named=True)
        assert m_row["LOG10P_ACAT"] == pytest.approx(s_row["LOG10P_ACAT"], rel=1e-9, abs=1e-12)
        assert m_row["N_SINGLEVAR_COMBINED"] == s_row["N_SINGLEVAR_COMBINED"]

    def test_burden_only_unit_with_non_add_test_emitted(self, tmp_path):
        """A burden-only unit (no annotation group) is emitted regardless of its
        TEST value — the burden-only path keys on (chr, unit, test), not a
        hardcoded 'ADD'."""
        base = tmp_path
        base.mkdir(parents=True, exist_ok=True)
        chrom = "1"
        # GENE_X has singles + annotations (a real annotation group, TEST=ADD).
        # GENE_Z has only a burden entry with a non-ADD TEST -> burden-only path.
        pl.DataFrame([{"ID": f"{chrom}:1000:A:C", "LOG10P": 3.0, "A1FREQ": 0.0004}]).write_csv(
            base / "singles.tsv", separator="\t"
        )
        pl.DataFrame(
            [
                {"ID": "GENE_X.mask1.0.001", "LOG10P": 2.0, "A1FREQ": 0.0003, "TEST": "ADD", "CHROM": chrom},
                {
                    "ID": "GENE_Z.mask1.0.001",
                    "LOG10P": 2.5,
                    "A1FREQ": 0.0003,
                    "TEST": "ADD-EXTRA",
                    "CHROM": chrom,
                },
            ]
        ).write_csv(base / "burden.tsv", separator="\t")
        with open(base / "masks.tsv", "w") as fh:
            fh.write("mask1\tpLoF\n")
        with open(base / "annots.tsv", "w") as fh:
            fh.write(f"{chrom}:1000:A:C\tGENE_X\tpLoF\n")

        out = base / "burden_only.tsv"
        compute_acat(
            single_variants_file=base / "singles.tsv",
            burden_file=base / "burden.tsv",
            masks_tsvs=[base / "masks.tsv"],
            annotations_tsvs=[base / "annots.tsv"],
            output_file=out,
            aggregate_freqs=["0.001"],
            log10p_col="LOG10P",
            weighting_scheme="standard",
            n_threads=1,
        )
        ids = set(pl.read_csv(str(out), separator="\t", null_values=["NA"])["ID"].to_list())
        # The non-ADD burden-only unit must still produce an aggregate row.
        assert "GENE_Z.mask1.0.001" in ids
        assert "GENE_X.mask1.0.001" in ids


class TestPerBinSingles:
    """Singles with an AAF_BIN column (``fedcat meta --aaf-bins``): each row is used
    only for its own bin, whatever its A1FREQ."""

    def _run(self, tmp_path, freqs):
        """Run ACAT on binned single-variant results for the given frequency bins."""
        gene, chrom = "GENE_X", "1"
        v1, v2 = f"{chrom}:1000:A:C", f"{chrom}:2000:A:G"
        singles = tmp_path / "singles.tsv"
        pl.DataFrame(
            [
                # v1: rare for 0.001 in some strata (lp 3), more strata enter 0.01 (lp 5)
                {"ID": v1, "LOG10P": 3.0, "A1FREQ": 0.0004, "AAF_BIN": "0.001"},
                {"ID": v1, "LOG10P": 5.0, "A1FREQ": 0.002, "AAF_BIN": "0.01"},
                # v2: no stratum rare enough for 0.001, though its 0.01 meta A1FREQ is
                {"ID": v2, "LOG10P": 4.0, "A1FREQ": 0.0005, "AAF_BIN": "0.01"},
                # a bin not requested is ignored
                {"ID": v2, "LOG10P": 9.0, "A1FREQ": 0.03, "AAF_BIN": "0.05"},
            ]
        ).write_csv(singles, separator="\t")
        masks = tmp_path / "masks.tsv"
        masks.write_text("mask1\tpLoF\n")
        annots = tmp_path / "annots_chr1.tsv"
        annots.write_text(f"{v1}\t{gene}\tpLoF\n{v2}\t{gene}\tpLoF\n")
        out = tmp_path / "acat.tsv"
        compute_acat(
            single_variants_file=singles,
            masks_tsvs=[masks],
            annotations_tsvs=[annots],
            output_file=out,
            aggregate_freqs=freqs,
            n_threads=1,
        )
        df = pl.read_csv(str(out), separator="\t", null_values=["NA"], schema_overrides={"THRESHOLD": pl.Utf8})
        return {r["THRESHOLD"]: r for r in df.iter_rows(named=True)}

    def test_each_row_only_feeds_its_bin(self, tmp_path):
        rows = self._run(tmp_path, ["0.001", "0.01"])
        assert rows["0.001"]["N_SINGLEVAR_COMBINED"] == 1
        # a single component's ACAT is its own p-value: v1's 0.001-bin result, not its 0.01 one
        assert rows["0.001"]["LOG10P_ACAT"] == pytest.approx(3.0, rel=1e-6)
        assert rows["0.01"]["N_SINGLEVAR_COMBINED"] == 2

    def test_bin_labels_match_on_value(self, tmp_path):
        """A config spelling the bin differently from the meta still picks up its rows."""
        rows = self._run(tmp_path, ["1e-3"])
        assert rows["1e-3"]["N_SINGLEVAR_COMBINED"] == 1


class TestMaskBounds:
    """START/END of ACAT rows when TSV annotations are filtered to variants with
    single-variant results. A mask made only of variants without results (e.g.
    MAC <= 10) still gets a burden-only ACAT row and needs positions."""

    def _write_inputs(self, base: Path):
        """Write inputs where one mask has variants without results and another has none with results."""
        base.mkdir(parents=True, exist_ok=True)
        # pLoF: two variants with results, one without; missense: none with results
        singles = base / "singles.tsv"
        pl.DataFrame(
            [{"ID": v, "LOG10P": lp, "A1FREQ": 0.0005} for v, lp in [("1:2000:A:C", 3.0), ("1:3000:A:G", 4.0)]]
        ).write_csv(singles, separator="\t")
        burden = base / "burden.tsv"
        pl.DataFrame(
            [
                {"ID": f"{unit}.{mask}.0.001", "LOG10P": lp, "A1FREQ": 0.0003, "TEST": "ADD", "CHROM": "1"}
                for unit, mask, lp in [("GENE_X", "lof", 2.0), ("GENE_X", "mis", 5.0), ("GENE_Z", "lof", 1.0)]
            ]
        ).write_csv(burden, separator="\t")
        masks = base / "masks.tsv"
        masks.write_text("lof\tpLoF\nmis\tmissense\n")
        annots = base / "annots.tsv"
        annots.write_text(
            "1:1000:A:T\tGENE_X\tpLoF\n"
            "1:2000:A:C\tGENE_X\tpLoF\n"
            "1:3000:A:G\tGENE_X\tpLoF\n"
            "1:5000:C:T\tGENE_X\tmissense\n"
            "1:9000:G:A\tGENE_X\tmissense\n"
        )
        return singles, burden, masks, annots

    def _run(self, tmp_path, full_bounds: bool) -> dict[str, dict]:
        """Run ACAT from TSV annotations, returning each output row by ID."""
        singles, burden, masks, annots = self._write_inputs(tmp_path / "in")
        out = tmp_path / f"acat_{full_bounds}.tsv"
        compute_acat(
            single_variants_file=singles,
            burden_file=burden,
            masks_tsvs=[masks],
            annotations_tsvs=[annots],
            output_file=out,
            log10p_col="LOG10P",
            n_threads=1,
            full_bounds=full_bounds,
        )
        df = pl.read_csv(str(out), separator="\t", null_values=["NA"])
        return {r["ID"]: r for r in df.iter_rows(named=True)}

    def test_full_bounds_cover_all_annotated_variants(self, tmp_path):
        rows = self._run(tmp_path, full_bounds=True)
        assert (rows["GENE_X.lof.0.001"]["START"], rows["GENE_X.lof.0.001"]["END"]) == (1000, 3000)
        assert (rows["GENE_X.mis.0.001"]["START"], rows["GENE_X.mis.0.001"]["END"]) == (5000, 9000)

    def test_filtered_bounds_leave_burden_only_masks_na(self, tmp_path):
        rows = self._run(tmp_path, full_bounds=False)
        assert (rows["GENE_X.lof.0.001"]["START"], rows["GENE_X.lof.0.001"]["END"]) == (2000, 3000)
        assert rows["GENE_X.mis.0.001"]["START"] is None
        assert rows["GENE_X.mis.0.001"]["END"] is None

    def test_unannotated_unit_is_na(self, tmp_path):
        rows = self._run(tmp_path, full_bounds=True)
        assert rows["GENE_Z.lof.0.001"]["START"] is None
        assert rows["GENE_Z.lof.0.001"]["END"] is None

    def test_bounds_scan_does_not_change_pvalues(self, tmp_path):
        full = self._run(tmp_path / "full", full_bounds=True)
        filtered = self._run(tmp_path / "filtered", full_bounds=False)
        assert set(full) == set(filtered)
        for k in full:
            assert full[k]["LOG10P_ACAT"] == pytest.approx(filtered[k]["LOG10P_ACAT"], rel=1e-12)
            assert full[k]["N_SINGLEVAR_COMBINED"] == filtered[k]["N_SINGLEVAR_COMBINED"]


def _write_tiny_schema(out_dir: Path, rows: list[tuple[str, str, str]], masks: dict[str, list[str]]) -> None:
    """Parquet schema with *rows* (variant, unit, annot) written in the order given."""
    out_dir.mkdir(parents=True, exist_ok=True)
    vids = sorted({r[0] for r in rows})
    units = sorted({r[1] for r in rows})
    annots = sorted({r[2] for r in rows} | {a for annots in masks.values() for a in annots})
    code = {name: i for names in (vids, units, annots) for i, name in enumerate(names)}
    pl.DataFrame({"VARIANT_CODE": range(len(vids)), "VARIANT_ID": vids}).write_parquet(
        out_dir / "dim_variants.parquet"
    )
    pl.DataFrame({"UNIT_CODE": range(len(units)), "UNIT": units}).write_parquet(out_dir / "dim_units.parquet")
    pl.DataFrame({"ANNOT_CODE": range(len(annots)), "ANNOT": annots}).write_parquet(
        out_dir / "dim_annots.parquet"
    )
    pl.DataFrame(
        {
            "VARIANT_CODE": [code[v] for v, _, _ in rows],
            "UNIT_CODE": [code[u] for _, u, _ in rows],
            "ANNOT_CODE": [code[a] for _, _, a in rows],
            "CHR": [v.split(":")[0] for v, _, _ in rows],
            "POS": [int(v.split(":")[1]) for v, _, _ in rows],
        }
    ).write_parquet(out_dir / "annotations.parquet")
    pl.DataFrame(
        {"UNIT": [], "CHR": [], "MASK": [], "MIN_POS": [], "MAX_POS": []},
        schema={"UNIT": pl.Utf8, "CHR": pl.Utf8, "MASK": pl.Utf8, "MIN_POS": pl.Int64, "MAX_POS": pl.Int64},
    ).write_parquet(out_dir / "unit_mask_bounds.parquet")
    pl.DataFrame(
        [{"MASK": m, "ANNOT_CODE": code[a], "ANNOT": a} for m, ma in masks.items() for a in ma],
        schema={"MASK": pl.Utf8, "ANNOT_CODE": pl.Int64, "ANNOT": pl.Utf8},
    ).write_parquet(out_dir / "mask_annotations.parquet")


def _tiny_singles(path: Path, vids: list[str]) -> Path:
    """Write single-variant results of `vids`, all rare and nominally significant."""
    pl.DataFrame({"ID": vids, "A1FREQ": [0.0005] * len(vids), "LOG10P": [3.0] * len(vids)}).write_csv(
        path, separator="\t"
    )
    return path


def test_sweep_refuses_a_unit_split_across_the_annotations(tmp_path):
    """A unit that reappears after another would be ACAT'd as two partial units; it must fail instead."""
    rows = [("1:100:A:G", "U1", "lof"), ("1:200:A:G", "U2", "lof"), ("1:300:A:G", "U1", "lof")]
    _write_tiny_schema(tmp_path / "schema", rows, {"M1": ["lof"]})
    singles = _tiny_singles(tmp_path / "singles.tsv", [r[0] for r in rows])
    with pytest.raises(BaseException, match="not correctly sorted"):
        compute_acat_streaming_rs(
            str(singles), "", str(tmp_path / "schema"), ["0.001"], output_path="", backend="sweep"
        )


def test_arrow_backend_refuses_more_than_64_masks(tmp_path):
    """Masks are bits of a u64 in the arrow backend; mask 64 would alias mask 0."""
    rows = [("1:100:A:G", "U1", "lof")]
    _write_tiny_schema(tmp_path / "schema", rows, {f"M{i}": ["lof"] for i in range(65)})
    singles = _tiny_singles(tmp_path / "singles.tsv", ["1:100:A:G"])
    with pytest.raises(BaseException, match="at most 64"):
        compute_acat_streaming_rs(
            str(singles), "", str(tmp_path / "schema"), ["0.001"], output_path="", backend="arrow"
        )
    # the sweep backend handles it
    compute_acat_streaming_rs(
        str(singles), "", str(tmp_path / "schema"), ["0.001"], output_path="", backend="sweep"
    )
