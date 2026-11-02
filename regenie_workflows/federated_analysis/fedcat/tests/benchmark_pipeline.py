#!/usr/bin/env python3
"""Standalone benchmark: meta-analysis -> ACAT pipeline.

Generates dummy REGENIE files for N cohorts × (1M singles, 50K aggregates),
then runs meta-analysis and ACAT with/without Rust across TSV and
optimised annotation backends.  Each benchmark case runs in an isolated
forked subprocess so RSS memory tracking is accurate.

Usage::

    python fedcat/tests/benchmark_pipeline.py
    python fedcat/tests/benchmark_pipeline.py --regenerate   # force re-create data
    python fedcat/tests/benchmark_pipeline.py --skip-datagen # skip data gen entirely

Files are cached in /tmp/fedcat_benchmark_1m/ and reused across runs.
"""

from __future__ import annotations

import argparse
import functools
import gc
import gzip
import io
import math
import multiprocessing as _mp
import os
import random
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.parquet as pq

import fedcat_rs

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

os.environ.setdefault("WORKSPACE_BUCKET", "gs://dummy-bucket")

import fedcat.stats.acat as _acat_mod
import fedcat.stats.meta_analysis as _meta_mod
from fedcat.stats.acat import compute_acat
from fedcat.stats.meta_analysis import meta_analyse_regenie

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BENCH_DIR = Path("/tmp/fedcat_benchmark_5m")
N_SINGLE_VARIANTS = 5_000_000
N_AGGREGATES = 1_000_000
N_COHORTS = 18
COHORT_NAMES = [f"cohort_{chr(ord('A') + i)}" for i in range(N_COHORTS)]
N_MASKS = 12
MASK_NAMES = [f"mask{i}" for i in range(1, N_MASKS + 1)]
N_ANNOTS = 48
ANNOT_NAMES = [f"annot{i}" for i in range(1, N_ANNOTS + 1)]
CHRS = [str(c) for c in range(1, 23)] + ["X"]
SEED = 42
N_CHUNKS_PER_COHORT = 24  # genome chunks per cohort
CHUNK_OVERLAP_FRAC = 0.05  # fraction of rows duplicated between adjacent chunks

# Input compression formats to benchmark (same data, different encodings).
INPUT_FORMATS = [
    ("gz", ".regenie.gz"),
    ("zst", ".regenie.zst"),
]

# Randomly assign 1-5 annotations to each mask (seeded for reproducibility)
_mask_rng = random.Random(SEED)
MASK_ANNOTS = {mask: _mask_rng.sample(ANNOT_NAMES, _mask_rng.randint(1, min(5, N_ANNOTS))) for mask in MASK_NAMES}


# ---------------------------------------------------------------------------
# Timing infrastructure
# ---------------------------------------------------------------------------


class TimingTracker:
    """Accumulates wall-clock time for wrapped functions."""

    def __init__(self):
        """Start with no timings."""
        self.timings: dict[str, list[float]] = {}

    def wrap(self, func, name: str | None = None):
        """Return `func` wrapped to record the time of each call under `name`."""
        fname = name or func.__qualname__
        self.timings.setdefault(fname, [])

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            t0 = time.perf_counter()
            result = func(*args, **kwargs)
            self.timings[fname].append(time.perf_counter() - t0)
            return result

        return wrapper

    def reset(self):
        """Clear the timings."""
        self.timings.clear()

    def total(self, name: str) -> float:
        """Return the total seconds recorded for `name`."""
        return sum(self.timings.get(name, []))

    def calls(self, name: str) -> int:
        """Return the calls recorded for `name`."""
        return len(self.timings.get(name, []))

    def summary_table(self) -> list[dict]:
        """Return the calls, total and mean time of each timed function."""
        rows = []
        for name, times in sorted(self.timings.items()):
            if not times:
                continue
            rows.append(
                {
                    "function": name,
                    "calls": len(times),
                    "total_s": sum(times),
                    "mean_ms": 1000 * sum(times) / len(times),
                }
            )
        return rows


import ctypes
import ctypes.util
import threading

# ---------------------------------------------------------------------------
# Portable RSS measurement
# ---------------------------------------------------------------------------

if sys.platform == "darwin":
    _libc = ctypes.CDLL(ctypes.util.find_library("c"))
    _libc.task_info.argtypes = [
        ctypes.c_uint32,
        ctypes.c_int32,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_uint32),
    ]
    _libc.task_info.restype = ctypes.c_int32

    _mach_task_self = ctypes.c_uint32.in_dll(_libc, "mach_task_self_")
    _MACH_TASK_BASIC_INFO = 20

    class _time_value_t(ctypes.Structure):
        _fields_ = [("seconds", ctypes.c_int32), ("microseconds", ctypes.c_int32)]

    class _mach_task_basic_info(ctypes.Structure):
        _fields_ = [
            ("virtual_size", ctypes.c_uint64),
            ("resident_size", ctypes.c_uint64),
            ("resident_size_max", ctypes.c_uint64),
            ("user_time", _time_value_t),
            ("system_time", _time_value_t),
            ("policy", ctypes.c_int32),
            ("suspend_count", ctypes.c_int32),
        ]

    def _get_rss_bytes() -> int:
        """Get current RSS in bytes via mach MACH_TASK_BASIC_INFO (flavor 20)."""
        info = _mach_task_basic_info()
        count = ctypes.c_uint32(ctypes.sizeof(info) // 4)
        ret = _libc.task_info(
            _mach_task_self.value,
            _MACH_TASK_BASIC_INFO,
            ctypes.byref(info),
            ctypes.byref(count),
        )
        return info.resident_size if ret == 0 else 0

elif sys.platform == "linux":
    _page_size = os.sysconf("SC_PAGE_SIZE")

    def _get_rss_bytes() -> int:
        """Get current RSS in bytes from /proc/self/statm."""
        try:
            with open("/proc/self/statm") as f:
                return int(f.read().split()[1]) * _page_size
        except (OSError, ValueError):
            return 0

else:

    def _get_rss_bytes() -> int:
        return 0


class RSSSampler:
    """Track peak RSS delta during a code section using a background sampling thread.

    Usage::

        with RSSSampler() as sampler:
            do_work()
        print(f"Peak delta: {sampler.peak_delta_mb:.0f} MB")
    """

    def __init__(self, interval_s: float = 0.02):
        """Set the sampling interval."""
        self._interval = interval_s
        self._baseline = 0
        self._peak = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def __enter__(self):
        """Start sampling RSS in a background thread."""
        gc.collect()
        self._baseline = _get_rss_bytes()
        self._peak = self._baseline
        self._stop.clear()
        self._thread = threading.Thread(target=self._sample, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc):
        """Stop sampling."""
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2.0)
        final = _get_rss_bytes()
        self._peak = max(self._peak, final)

    def _sample(self):
        """Record the peak RSS until stopped."""
        while not self._stop.is_set():
            rss = _get_rss_bytes()
            self._peak = max(self._peak, rss)
            self._stop.wait(self._interval)

    @property
    def peak_rss_mb(self) -> float:
        """Absolute peak RSS in MB during the tracked section."""
        return self._peak / (1024 * 1024)

    @property
    def peak_delta_mb(self) -> float:
        """Peak RSS increase from baseline in MB."""
        return max(0, self._peak - self._baseline) / (1024 * 1024)


def _fork_child(conn, func, args, kwargs):
    """Child entry point for :func:`_fork_run`.

    Uses ``resource.getrusage(RUSAGE_SELF).ru_maxrss`` for memory tracking.
    This is a kernel-maintained high-water mark that is accurate even when
    Rust FFI code holds the GIL (which blocks Python sampling threads).
    """
    gc.collect()
    # Baseline: ru_maxrss after fresh process startup (spawn context)
    baseline_maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result = func(*args, **kwargs)
    peak_maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS: ru_maxrss is bytes; Linux: ru_maxrss is KB
    if sys.platform == "darwin":
        peak_mb = peak_maxrss / (1024 * 1024)
        baseline_mb = baseline_maxrss / (1024 * 1024)
    else:
        peak_mb = peak_maxrss / 1024
        baseline_mb = baseline_maxrss / 1024
    conn.send(
        {
            "result": result,
            "peak_rss_mb": peak_mb,
            "delta_rss_mb": max(0.0, peak_mb - baseline_mb),
        }
    )
    conn.close()


def _fork_run(func, *args, timeout_s=None, **kwargs):
    """Run *func* in a spawned subprocess with kernel-tracked RSS.

    Uses ``spawn`` (not ``fork``) to avoid macOS deadlocks caused by Polars'
    Rayon thread pool holding system locks at fork time.

    Returns ``(func_return_value, peak_rss_mb, delta_rss_mb)``.

    Args:
        timeout_s: If set, kill the child after this many seconds and raise
            :class:`TimeoutError`.
    """
    ctx = _mp.get_context("spawn")
    parent_conn, child_conn = ctx.Pipe(duplex=False)
    p = ctx.Process(target=_fork_child, args=(child_conn, func, args, kwargs))
    p.start()
    child_conn.close()  # parent only reads
    if timeout_s is not None:
        if parent_conn.poll(timeout_s):
            data = parent_conn.recv()
        else:
            p.kill()
            p.join()
            raise TimeoutError(f"Subprocess timed out after {timeout_s}s")
    else:
        try:
            data = parent_conn.recv()
        except EOFError:
            p.join()
            raise RuntimeError(f"Subprocess crashed (exit code {p.exitcode})")
    p.join()
    if p.exitcode != 0:
        raise RuntimeError(f"Subprocess exited with code {p.exitcode}")
    return data["result"], data["peak_rss_mb"], data["delta_rss_mb"]


# ---------------------------------------------------------------------------
# Data generation helpers
# ---------------------------------------------------------------------------


def _chrom_rank_expr() -> pl.Expr:
    """Map CHROM to a numeric sort rank: 1..22, X=23, Y=24, MT=25, else=26."""
    numeric = pl.col("CHROM").cast(pl.UInt32, strict=False)
    return (
        pl.when(numeric.is_not_null())
        .then(numeric)
        .when(pl.col("CHROM") == "X")
        .then(pl.lit(23, dtype=pl.UInt32))
        .when(pl.col("CHROM") == "Y")
        .then(pl.lit(24, dtype=pl.UInt32))
        .when(pl.col("CHROM") == "MT")
        .then(pl.lit(25, dtype=pl.UInt32))
        .otherwise(pl.lit(26, dtype=pl.UInt32))
    )


def _write_chunked_regenie_gz(
    base_path: Path,
    df: pl.DataFrame,
    n_chunks: int = N_CHUNKS_PER_COHORT,
    overlap_frac: float = CHUNK_OVERLAP_FRAC,
) -> list[Path]:
    """Split a sorted DataFrame into overlapping chunk files.

    Splits by row count into *n_chunks* pieces.  Adjacent chunks share
    *overlap_frac* of rows at their boundaries (duplicated).  Each chunk
    file is independently sorted.

    Returns list of written chunk file paths.
    """
    df = df.with_columns(_chrom_rank_expr().alias("_cr")).sort(["_cr", "GENPOS", "ID"]).drop("_cr")
    n = df.height
    base_size = n // n_chunks
    overlap = max(1, int(base_size * overlap_frac))

    paths = []
    for i in range(n_chunks):
        start = i * base_size
        if i == n_chunks - 1:
            end = n
        else:
            end = (i + 1) * base_size + overlap  # overlap with next chunk
        chunk_df = df.slice(start, end - start)
        stem = base_path.stem.replace(".regenie", "")
        chunk_path = base_path.parent / f"{stem}_chunk_{i}.regenie.gz"
        buf = io.BytesIO()
        chunk_df.write_csv(buf, separator=" ", include_header=True)
        with gzip.open(chunk_path, "wb", compresslevel=1) as f:
            f.write(buf.getvalue())
        paths.append(chunk_path)
    return paths


def _generate_singles_df(n: int, rng: np.random.Generator) -> pl.DataFrame:
    """Generate *n* single-variant REGENIE rows (vectorised)."""
    chrs = rng.choice(CHRS, n)
    positions = rng.integers(1, 250_000_000, n)
    bases = np.array(["A", "C", "G", "T"])
    ref_idx = rng.integers(0, 4, n)
    alt_idx = (ref_idx + rng.integers(1, 4, n)) % 4
    afs = rng.uniform(0.00001, 0.001, n)
    betas = rng.normal(0, 0.05, n)
    ses = rng.uniform(0.005, 0.05, n)
    z = betas / ses

    df = pl.DataFrame(
        {
            "CHROM": chrs,
            "GENPOS": positions.astype(np.int64),
            "ALLELE0": bases[ref_idx],
            "ALLELE1": bases[alt_idx],
            "A1FREQ": np.round(afs, 6),
            "N": np.full(n, 50_000, dtype=np.int64),
            "TEST": np.full(n, "ADD", dtype=object),
            "BETA": np.round(betas, 6),
            "SE": np.round(ses, 6),
            "CHISQ": np.round(z**2, 4),
            "LOG10P": np.round(np.clip(np.abs(z), 0.01, 50.0), 4),
            "EXTRA": np.full(n, "NA", dtype=object),
        }
    )

    # Build ID column in polars (fast Rust strings)
    df = df.with_columns(
        (
            pl.col("CHROM")
            + pl.lit(":")
            + pl.col("GENPOS").cast(pl.Utf8)
            + pl.lit(":")
            + pl.col("ALLELE0")
            + pl.lit(":")
            + pl.col("ALLELE1")
        ).alias("ID")
    )

    return df.select(
        [
            "CHROM",
            "GENPOS",
            "ID",
            "ALLELE0",
            "ALLELE1",
            "A1FREQ",
            "N",
            "TEST",
            "BETA",
            "SE",
            "CHISQ",
            "LOG10P",
            "EXTRA",
        ]
    )


def _generate_burden_df(
    units: np.ndarray,
    unit_chrs: np.ndarray,
    unit_genpos: np.ndarray,
    rng: np.random.Generator,
    valid_unit_masks: set[tuple[str, str]] | None = None,
) -> pl.DataFrame:
    """Generate burden/MAC10 rows for annotated units x masks only (vectorised)."""
    n_masks = len(MASK_NAMES)
    n_total = len(units) * n_masks

    chr_arr = np.repeat(unit_chrs, n_masks)
    unit_arr = np.repeat(units, n_masks)
    mask_arr = np.tile(MASK_NAMES, len(units))
    genpos_arr = np.repeat(unit_genpos, n_masks)

    afs = rng.uniform(0.00005, 0.0005, n_total)
    betas = rng.normal(0, 0.1, n_total)
    ses = rng.uniform(0.01, 0.1, n_total)
    z = betas / ses

    # Build burden IDs inline with polars later
    df = pl.DataFrame(
        {
            "CHROM": chr_arr,
            "GENPOS": genpos_arr,
            "UNIT": unit_arr,
            "MASK": mask_arr,
            "ALLELE0": np.full(n_total, "A", dtype=object),
            "ALLELE1": np.full(n_total, "T", dtype=object),
            "A1FREQ": np.round(afs, 6),
            "N": np.full(n_total, 50_000, dtype=np.int64),
            "TEST": np.full(n_total, "ADD", dtype=object),
            "BETA": np.round(betas, 6),
            "SE": np.round(ses, 6),
            "CHISQ": np.round(z**2, 4),
            "LOG10P": np.round(np.clip(np.abs(z), 0.01, 50.0), 4),
            "EXTRA": np.full(n_total, "NA", dtype=object),
        }
    )

    df = df.with_columns((pl.col("UNIT") + pl.lit(".") + pl.col("MASK") + pl.lit(".0.001")).alias("ID"))

    if valid_unit_masks is not None:
        vm_df = pl.DataFrame({"UNIT": [u for u, m in valid_unit_masks], "MASK": [m for u, m in valid_unit_masks]})
        df = df.join(vm_df, on=["UNIT", "MASK"], how="inner")

    return df.select(
        [
            "CHROM",
            "GENPOS",
            "ID",
            "ALLELE0",
            "ALLELE1",
            "A1FREQ",
            "N",
            "TEST",
            "BETA",
            "SE",
            "CHISQ",
            "LOG10P",
            "EXTRA",
        ]
    )


def _generate_annotations(
    singles_df: pl.DataFrame,
    units: np.ndarray,
    unit_chrs: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Assign singles to units: 0-1000 per unit (geometric ~100 mean).

    Returns ``(variant_ids, unit_names, annotation_names)`` aligned arrays.
    A single variant can appear in multiple units through shared chromosome.
    Within each unit, variants share annotations that allow them to appear in
    multiple masks (e.g. pLoF appears in mask1, mask2, mask3).
    """
    chr_col = singles_df["CHROM"].to_numpy()
    id_col = singles_df["ID"].to_numpy()

    # Index singles by chromosome
    chr_to_ids: dict[str, np.ndarray] = {}
    for c in CHRS:
        mask = chr_col == c
        chr_to_ids[c] = id_col[mask]

    # Per-unit count: geometric(p=0.01) → mean 100, cap at 1000
    n_per_unit = np.minimum(rng.geometric(p=0.01, size=len(units)), 1000)
    total_expected = int(n_per_unit.sum())
    print(f"    annotation tuples expected: ~{total_expected:,}")

    annot_choices = np.array(ANNOT_NAMES)

    # Pre-allocate
    all_vids: list[np.ndarray] = []
    all_units: list[np.ndarray] = []
    all_annots: list[np.ndarray] = []

    for i in range(len(units)):
        chrom = unit_chrs[i]
        available = chr_to_ids.get(chrom)
        if available is None or len(available) == 0:
            continue
        n_pick = min(int(n_per_unit[i]), len(available))
        if n_pick == 0:
            continue
        chosen_idx = rng.choice(len(available), n_pick, replace=False)
        all_vids.append(available[chosen_idx])
        all_units.append(np.full(n_pick, units[i]))
        all_annots.append(rng.choice(annot_choices, n_pick))

    return (
        np.concatenate(all_vids),
        np.concatenate(all_units),
        np.concatenate(all_annots),
    )


# ---------------------------------------------------------------------------
# Annotation writers
# ---------------------------------------------------------------------------


def _write_tsv_annotations(
    vids: np.ndarray,
    unit_names: np.ndarray,
    annot_names: np.ndarray,
    out_dir: Path,
) -> list[Path]:
    """Write TSV annotation files (one per chromosome)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    # Extract chromosome from variant ID  (chr:pos:ref:alt)
    chrs = np.array([v.split(":")[0] for v in vids])
    paths = []
    for c in sorted(set(chrs)):
        mask = chrs == c
        path = out_dir / f"annotations_chr{c}.tsv"
        with open(path, "w") as fh:
            fh.writelines(
                f"{vid}\t{unit}\t{annot}\n"
                for vid, unit, annot in zip(vids[mask], unit_names[mask], annot_names[mask])
            )
        paths.append(path)
    return sorted(paths)


def _write_masks_file(path: Path) -> None:
    """Write the benchmark masks as a REGENIE mask file."""
    with open(path, "w") as fh:
        fh.writelines(f"{mask_name}\t{','.join(annot_list)}\n" for mask_name, annot_list in MASK_ANNOTS.items())


def _create_format_variants(regenie_dir: Path) -> None:
    """Create zstd copies of every .regenie.gz file.

    Allows benchmarking the same data across gz and zstd input compressions
    to separate decompression overhead from computation.
    """
    gz_files = sorted(regenie_dir.glob("*.regenie.gz"))
    n_files = len(gz_files)
    for i, gz_path in enumerate(gz_files, 1):
        with gzip.open(gz_path, "rb") as f:
            data = f.read()

        # Zstd .regenie.zst
        zst_path = Path(str(gz_path).replace(".regenie.gz", ".regenie.zst"))
        subprocess.run(
            ["zstd", "-3", "-f", "--no-progress", "-o", str(zst_path)],
            input=data,
            capture_output=True,
            check=True,
        )

        if i % 50 == 0 or i == n_files:
            print(f"         {i}/{n_files} files")


# ---------------------------------------------------------------------------
# Data generation orchestration
# ---------------------------------------------------------------------------


def generate_data(force: bool = False) -> dict:
    """Generate (or reuse) all benchmark input files.  Returns path dict."""
    bench_dir = BENCH_DIR
    bench_dir.mkdir(parents=True, exist_ok=True)
    marker = bench_dir / ".generated"

    if marker.exists() and not force:
        print(f"Reusing cached data in {bench_dir}")
        return _collect_paths(bench_dir)

    print(f"\n{'=' * 70}")
    print(
        f"  GENERATING BENCHMARK DATA  ({N_SINGLE_VARIANTS:,} singles, "
        f"{N_AGGREGATES:,} aggregates, {N_COHORTS} cohorts)"
    )
    print(f"{'=' * 70}")
    t0 = time.perf_counter()

    rng_master = np.random.default_rng(SEED)
    cohort_rngs = [np.random.default_rng(SEED + i + 1) for i in range(N_COHORTS)]

    # 1. Generate a single reference set of variant IDs (shared across cohorts
    #    with different effect sizes) — avoids generating 6x 10M ID strings
    print("  [1/6] Generating single-variant IDs ...")
    ref_singles = _generate_singles_df(N_SINGLE_VARIANTS, rng_master)
    # Dedup by ID to avoid collisions (random CHROM:POS:REF:ALT can collide)
    ref_singles = ref_singles.unique(subset=["ID"], keep="first")
    print(f"         done ({ref_singles.height:,} unique IDs, {time.perf_counter() - t0:.1f}s)")

    # 2. Gene units
    print("  [2/6] Generating gene units ...")
    unit_chrs = rng_master.choice(CHRS, N_AGGREGATES)
    units = np.array([f"GENE{i:06d}" for i in range(N_AGGREGATES)])
    # Fixed GENPOS per unit — must be shared across cohorts so meta grouping works
    unit_genpos = rng_master.integers(1, 250_000_000, N_AGGREGATES).astype(np.int64)
    print(f"         {N_AGGREGATES:,} units across {len(set(unit_chrs))} chromosomes")

    # 3. Annotations (0–1000 singles/unit, geometric mean ~100) — generated
    #    before burden writing so we can filter burden to annotated unit×mask pairs only.
    print("  [3/6] Generating annotations ...")
    ann_vids, ann_units, ann_annots = _generate_annotations(
        ref_singles,
        units,
        unit_chrs,
        rng_master,
    )
    print(f"         {len(ann_vids):,} tuples ({time.perf_counter() - t0:.1f}s)")

    # Build set of (unit, mask) pairs that have at least one annotation variant.
    # A burden result cannot exist in reality without annotated variants.
    valid_unit_masks: set[tuple[str, str]] = set()
    for _u, _a in zip(ann_units, ann_annots):
        for _mask_name, _mask_annots in MASK_ANNOTS.items():
            if _a in _mask_annots:
                valid_unit_masks.add((_u, _mask_name))

    # 4. Write REGENIE files per cohort (burden filtered to valid unit×mask pairs)
    regenie_dir = bench_dir / "regenie"
    regenie_dir.mkdir(exist_ok=True)

    print("  [4/6] Writing REGENIE files ...")
    n_singles_actual = ref_singles.height
    ref_afs = ref_singles["A1FREQ"].to_numpy()
    for ci, cohort in enumerate(COHORT_NAMES):
        crng = cohort_rngs[ci]
        # Per-cohort variant missingness: each cohort drops 5-15% of variants
        miss_frac = crng.uniform(0.05, 0.15)
        keep_mask = crng.random(n_singles_actual) >= miss_frac
        n_kept = int(keep_mask.sum())
        # Per-cohort MAF variation: proportional noise around reference AF
        kept_afs = ref_afs[keep_mask]
        cohort_afs = np.round(
            np.clip(kept_afs + crng.normal(0, kept_afs * 0.1, n_kept), 1e-6, 1.0 - 1e-6),
            6,
        )
        cdf = (
            ref_singles.filter(pl.lit(keep_mask))
            .with_columns(
                [
                    pl.lit(cohort_afs).alias("A1FREQ"),
                    pl.lit(np.round(crng.normal(0, 0.05, n_kept), 6)).alias("BETA"),
                    pl.lit(np.round(crng.uniform(0.005, 0.05, n_kept), 6)).alias("SE"),
                ]
            )
            .with_columns(
                [
                    (pl.col("BETA") / pl.col("SE")).pow(2).round(4).alias("CHISQ"),
                    (pl.col("BETA") / pl.col("SE")).abs().clip(0.01, 50.0).round(4).alias("LOG10P"),
                ]
            )
        )
        _write_chunked_regenie_gz(
            regenie_dir / f"{cohort}_singlevars_chr_all.regenie.gz",
            cdf,
        )
        # Burden + MAC10 (independent effect sizes, filtered to annotated unit×mask pairs)
        bdf = _generate_burden_df(units, unit_chrs, unit_genpos, crng, valid_unit_masks)
        b_miss = crng.uniform(0.05, 0.15)
        bdf = bdf.filter(pl.lit(crng.random(bdf.height) >= b_miss))
        _write_chunked_regenie_gz(
            regenie_dir / f"{cohort}_burden_chr_all.regenie.gz",
            bdf,
        )
        mdf = _generate_burden_df(units, unit_chrs, unit_genpos, crng, valid_unit_masks)
        m_miss = crng.uniform(0.05, 0.15)
        mdf = mdf.filter(pl.lit(crng.random(mdf.height) >= m_miss))
        _write_chunked_regenie_gz(
            regenie_dir / f"{cohort}_mac10_chr_all.regenie.gz",
            mdf,
        )
        print(f"         {cohort} written ({time.perf_counter() - t0:.1f}s)")

    # 4b. Create zstd variants of all .gz files
    print("         Creating zstd format variants ...")
    _create_format_variants(regenie_dir)
    print(f"         done ({time.perf_counter() - t0:.1f}s)")

    # 5. Write TSV annotation files
    print("  [5/6] Writing TSV annotations ...")
    _write_tsv_annotations(ann_vids, ann_units, ann_annots, bench_dir / "annots_tsv")
    print(f"         done ({time.perf_counter() - t0:.1f}s)")

    # 6. Masks file
    print("  [6/6] Writing masks file ...")
    _write_masks_file(bench_dir / "masks.tsv")

    marker.write_text(f"singles={N_SINGLE_VARIANTS} aggregates={N_AGGREGATES} cohorts={N_COHORTS} seed={SEED}\n")
    total = time.perf_counter() - t0
    print(f"\n  Data generation complete: {total:.1f}s\n")

    return _collect_paths(bench_dir)


def _collect_paths(bench_dir: Path) -> dict:
    """Return the benchmark's input and output paths."""
    return {
        "bench_dir": bench_dir,
        "regenie_dir": bench_dir / "regenie",
        "masks": bench_dir / "masks.tsv",
        "tsv_annots": sorted((bench_dir / "annots_tsv").glob("*.tsv")),
        "output_dir": bench_dir / "acat_output",
    }


# ---------------------------------------------------------------------------
# Benchmark runner
# ---------------------------------------------------------------------------


def _install_timing_hooks(tracker: TimingTracker) -> dict:
    """Monkey-patch key functions with timing wrappers.  Returns originals."""
    targets = {
        # acat.py
        "_convert_to_acat_input_rs": (_acat_mod, "rust_convert_to_acat_input"),
        "_compute_acat_streaming_rs": (_acat_mod, "rust_streaming_acat"),
        # meta_analysis.py
        "z_to_log10p": (_meta_mod, "z_to_log10p"),
    }
    originals = {}
    for attr, (mod, label) in targets.items():
        orig = getattr(mod, attr, None)
        if orig is None:
            continue
        originals[(mod, attr)] = orig
        setattr(mod, attr, tracker.wrap(orig, label))
    return originals


def _uninstall_timing_hooks(originals: dict) -> None:
    """Restore the functions timing hooks replaced."""
    for (mod, attr), orig in originals.items():
        setattr(mod, attr, orig)


def _acat_case(paths, backend, label, schema_dir=None):
    """Run ACAT pipeline (designed for subprocess isolation).

    If *schema_dir* is a non-empty string, it is used directly as the
    pre-converted parquet schema directory (output of the conversion
    benchmark).  Otherwise, TSVs are passed directly and ACAT converts
    internally.

    Returns a picklable metrics dict.
    """
    tracker = TimingTracker()
    originals = _install_timing_hooks(tracker)

    out_dir = paths["output_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    output_file = out_dir / f"acat_{label}.tsv"

    try:
        t0 = time.perf_counter()
        if schema_dir:
            compute_acat(
                single_variants_file=paths["meta_singles"],
                burden_file=paths["meta_mac10"],
                schema_dir=Path(schema_dir),
                output_file=output_file,
                aggregate_freqs=["0.001"],
                log10p_col="LOG10P",
                weighting_scheme="standard",
                n_threads=4,
                backend=backend,
            )
        else:
            compute_acat(
                single_variants_file=paths["meta_singles"],
                burden_file=paths["meta_mac10"],
                masks_tsvs=[paths["masks"]],
                annotations_tsvs=paths["tsv_annots"],
                output_file=output_file,
                aggregate_freqs=["0.001"],
                log10p_col="LOG10P",
                weighting_scheme="standard",
                n_threads=4,
                backend=backend,
            )
        elapsed = time.perf_counter() - t0
    finally:
        _uninstall_timing_hooks(originals)

    internal = tracker.summary_table()
    load_time = sum(r["total_s"] for r in internal if r["function"].startswith("load_"))
    index_time = sum(r["total_s"] for r in internal if r["function"].startswith("build_"))

    # Count result rows (subtract 1 for header line)
    try:
        with open(output_file) as fh:
            n_results = sum(1 for _ in fh) - 1
    except OSError:
        n_results = 0

    return {
        "elapsed": elapsed,
        "load_s": load_time,
        "index_s": index_time,
        "compute_s": elapsed - load_time - index_time,
        "n_results": n_results,
        "internal": internal,
    }


def _conversion_case(paths, method, label):
    """Run TSV→parquet conversion in subprocess isolation.

    *method*: 'python_polars' | 'rust_v1'
    Output is written to a persistent directory under bench_dir so it can
    be reused by the ACAT benchmark.  Returns a picklable metrics dict
    including ``out_dir`` (str) for downstream consumers.
    """
    out_dir = str(Path(paths["bench_dir"]) / f"conv_output_{label}")
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    tsv_files = [str(p) for p in paths["tsv_annots"]]
    masks_tsv = str(paths["masks"])

    t0 = time.perf_counter()
    match method:
        case "python_polars":
            from fedcat.stats.acat import (
                convert_tsv_to_parquet_polars,
            )

            convert_tsv_to_parquet_polars([Path(p) for p in tsv_files], [Path(masks_tsv)], Path(out_dir))
        case "rust_v1":
            fedcat_rs.convert_tsv_to_parquet_rs(tsv_files, [masks_tsv], out_dir)
        case _:
            raise ValueError(f"Unknown method: {method}")
    elapsed = time.perf_counter() - t0

    # Count output rows from annotations.parquet
    ann_pq = Path(out_dir) / "annotations.parquet"
    n_rows = pq.read_metadata(ann_pq).num_rows if ann_pq.exists() else 0

    return {"elapsed": elapsed, "n_rows": n_rows, "out_dir": out_dir}


def _meta_case(regenie_dir, meta_dir, backend, res_type, suffix, input_ext=".regenie.gz"):
    """Run a single meta-analysis type (designed for subprocess isolation).

    Returns a picklable metrics dict including the output file path.
    """
    tracker = TimingTracker()
    originals = _install_timing_hooks(tracker)
    meta_dir = Path(meta_dir)
    meta_dir.mkdir(exist_ok=True)

    try:
        cohort_files = {}
        for c in COHORT_NAMES:
            cohort_files[c] = sorted(Path(regenie_dir).glob(f"{c}_{suffix}_chr_all_chunk_*{input_ext}"))
        t0 = time.perf_counter()
        meta_analyse_regenie(
            per_cohort_files=cohort_files,
            output_file=meta_dir / f"META_{res_type}.tsv",
            res_type=res_type,
            backend=backend,
        )
        elapsed = time.perf_counter() - t0
    finally:
        _uninstall_timing_hooks(originals)

    return {
        "elapsed": elapsed,
        "internal": tracker.summary_table(),
        "output_path": str(meta_dir / f"META_{res_type}.tsv"),
    }


_META_TYPES = [
    ("single_variants", "singlevars"),
    ("burden", "burden"),
    ("mac10", "mac10"),
]

# Timeout for a single meta-analysis subprocess (seconds).
META_TIMEOUT_S = 1800


def _run_meta_suite(
    label: str,
    regenie_dir: Path,
    meta_dir: Path,
    backend: str,
    input_ext: str = ".regenie.gz",
) -> tuple[dict, dict]:
    """Run all 3 meta types in separate subprocesses.

    Returns ``(timings_dict, per_type_mem_dict)`` where each maps
    res_type → seconds / peak_rss_mb.  Also prints live progress.
    Meta output paths are stored in ``meta_dir``.
    """
    timings = {}
    mem = {}
    for res_type, suffix in _META_TYPES:
        print(f"    {res_type} ...", end=" ", flush=True)
        try:
            metrics, peak_mb, _delta_mb = _fork_run(
                _meta_case,
                regenie_dir,
                meta_dir,
                backend,
                res_type,
                suffix,
                input_ext=input_ext,
                timeout_s=META_TIMEOUT_S,
            )
            timings[res_type] = metrics["elapsed"]
            mem[res_type] = peak_mb
            print(f"{metrics['elapsed']:.1f}s | {peak_mb:.0f}MB peak")
        except TimeoutError:
            timings[res_type] = float("nan")
            mem[res_type] = float("nan")
            print(f"TIMEOUT (>{META_TIMEOUT_S}s)")
        except RuntimeError as e:
            timings[res_type] = float("nan")
            mem[res_type] = float("nan")
            print(f"FAILED: {e}")
    return timings, mem


def run_benchmark(paths: dict, low_memory_only: bool = False) -> None:
    """Execute the full benchmark suite and print results."""
    regenie_dir = paths["regenie_dir"]

    # Ensure format variants exist (needed for compression benchmarks)
    if not next(regenie_dir.glob("*_singlevars_chr_all_chunk_*.regenie.zst"), None):
        print("\n  Creating format variants (zstd) ...")
        t0 = time.perf_counter()
        _create_format_variants(regenie_dir)
        print(f"  Format variants created in {time.perf_counter() - t0:.1f}s")

    # ── Meta-analysis benchmark ──────────────────────────────────
    print(f"\n{'=' * 70}")
    print("  META-ANALYSIS BENCHMARK")
    print("  Each result type runs in a separate subprocess for accurate memory.")
    print(f"  Input formats: {', '.join(f for f, _ in INPUT_FORMATS)}")
    print(f"{'=' * 70}")

    meta_runs = []
    meta_output_paths = {}

    # ── Backends × input formats ─────────────────────────────────
    for backend in ["streaming"]:  # , "merge"]:
        for fmt_label, input_ext in INPUT_FORMATS:
            meta_dir = paths["bench_dir"] / f"meta_{backend}_{fmt_label}"
            print(f"\n  {backend} meta-analysis (input: {fmt_label}):")
            timings, mem_usage = _run_meta_suite(
                f"{backend}/{fmt_label}",
                regenie_dir,
                meta_dir,
                backend=backend,
                input_ext=input_ext,
            )
            meta_runs.append((f"{backend}/{fmt_label}", timings, mem_usage))
            total = sum(v for v in timings.values() if not math.isnan(v))
            print(f"  total: {total:.2f}s")

            # Take first complete output for ACAT
            if not meta_output_paths and all(p in timings and not math.isnan(timings[p]) for p, _ in _META_TYPES):
                meta_output_paths = {
                    "meta_singles": meta_dir / "META_single_variants.tsv",
                    "meta_mac10": meta_dir / "META_mac10.tsv",
                    "meta_burden": meta_dir / "META_burden.tsv",
                }
    if not meta_output_paths:
        print("\n  WARNING: No complete meta-analysis output available for ACAT.")

    # ── Conversion benchmark ─────────────────────────────────────
    CONV_TIMEOUT_S = 600
    conv_results: list[dict] = []

    conv_runs = [
        ("python_polars", "python_polars", "Python Polars"),
        ("rust_v1", "rust_v1", "Rust v1 (csv crate)"),
    ]

    print(f"\n{'=' * 70}")
    print("  TSV → PARQUET CONVERSION BENCHMARK")
    print("  Each converter runs in a separate subprocess for accurate memory.")
    print("  Output is written to bench_dir and reused by the ACAT benchmark.")
    print(f"{'=' * 70}")

    conv_paths = dict(paths)

    conv_schema_dirs: dict[str, str] = {}  # label -> persistent output dir

    for label, method, desc in conv_runs:
        print(f"\n  {desc} ...", end=" ", flush=True)
        try:
            metrics, peak_mb, _delta_mb = _fork_run(
                _conversion_case, conv_paths, method, label, timeout_s=CONV_TIMEOUT_S
            )
            conv_results.append(
                {
                    "label": label,
                    "desc": desc,
                    "elapsed": metrics["elapsed"],
                    "n_rows": metrics["n_rows"],
                    "peak_mem_mb": peak_mb,
                    "out_dir": metrics.get("out_dir", ""),
                }
            )
            conv_schema_dirs[label] = metrics.get("out_dir", "")
            print(f"{metrics['elapsed']:.2f}s | {peak_mb:.0f}MB peak | {metrics['n_rows']:,} rows")
        except TimeoutError:
            conv_results.append(
                {
                    "label": label,
                    "desc": desc,
                    "elapsed": float("nan"),
                    "n_rows": 0,
                    "peak_mem_mb": float("nan"),
                    "out_dir": "",
                }
            )
            print(f"TIMEOUT (>{CONV_TIMEOUT_S}s)")
        except RuntimeError as e:
            conv_results.append(
                {
                    "label": label,
                    "desc": desc,
                    "elapsed": float("nan"),
                    "n_rows": 0,
                    "peak_mem_mb": float("nan"),
                    "out_dir": "",
                }
            )
            print(f"FAILED: {e}")

    # ── ACAT benchmark ───────────────────────────────────────────
    acat_results = []
    if meta_output_paths:
        # Two pipeline variants:
        #   1. rust_sweep → use rust_v1 conversion output (schema_dir) from conv benchmark
        #   2. rust_tsv   → pass TSVs directly; ACAT converts internally
        rust_v1_schema = conv_schema_dirs.get("rust_v1", "")
        runs = [
            ("rust_sweep", "sweep", "parquet", rust_v1_schema),
            ("rust_arrow", "arrow", "parquet", rust_v1_schema),
            ("rust_tsv", "sweep", "tsv", ""),
        ]

        acat_paths = dict(paths)
        acat_paths.update(meta_output_paths)

        for label, backend, fmt, schema_dir in runs:
            if label == "rust_sweep" and not schema_dir:
                print(f"\n  ACAT {label} SKIPPED (rust_v1 conversion did not produce output)")
                continue
            print(f"\n  Running ACAT: {label} ...")

            metrics, peak_mb, _delta_mb = _fork_run(_acat_case, acat_paths, backend, label, schema_dir=schema_dir)

            acat_results.append(
                {
                    "backend": label,
                    "format": fmt,
                    "total_s": metrics["elapsed"],
                    "load_s": metrics["load_s"],
                    "index_s": metrics["index_s"],
                    "compute_s": metrics["compute_s"],
                    "peak_mem_mb": peak_mb,
                    "n_results": metrics["n_results"],
                    "internal": metrics["internal"],
                }
            )

            print(
                f"    {metrics['elapsed']:.2f}s total | {metrics['load_s']:.2f}s load | "
                f"{metrics['index_s']:.2f}s index | {metrics['compute_s']:.2f}s compute | "
                f"{peak_mb:.0f}MB peak | {metrics['n_results']:,} rows"
            )

    # ── Summary ──────────────────────────────────────────────────
    _print_summary(meta_runs, conv_results, acat_results)


def _print_summary(meta_runs: list[tuple], conv_results: list[dict], acat_results: list[dict]) -> None:
    """Print the detailed benchmark summary table."""

    print(f"\n\n{'=' * 80}")
    print("  BENCHMARK SUMMARY")
    print(
        f"  {N_SINGLE_VARIANTS:,} singles | {N_AGGREGATES:,} aggregates | "
        f"{N_COHORTS} cohorts | formats: {', '.join(f for f, _ in INPUT_FORMATS)}"
    )
    print(f"{'=' * 80}")

    # Meta-analysis section
    print("\n  META-ANALYSIS  (per-type subprocess isolation, kernel-tracked RSS)")
    print(f"  {'Backend':<18} {'Type':<22} {'Time (s)':<12} {'Peak RSS (MB)':<14}")
    print(f"  {'-' * 66}")
    for label, timings, mem in meta_runs:
        for rt, _ in _META_TYPES:
            t = timings.get(rt, float("nan"))
            m = mem.get(rt, float("nan"))
            t_s = f"{t:.2f}" if not math.isnan(t) else "TIMEOUT"
            m_s = f"{m:.0f}" if not math.isnan(m) else "—"
            print(f"  {label:<18} {rt:<22} {t_s:<12} {m_s:<14}")
        total = sum(v for v in timings.values() if not math.isnan(v))
        max_mem = max((v for v in mem.values() if not math.isnan(v)), default=0)
        print(f"  {label:<18} {'TOTAL':<22} {total:<12.2f} {max_mem:<14.0f}")
        print(f"  {'-' * 66}")

    # Meta speedup (compare streaming vs merge on gz format)
    meta_by_label = {label: (timings, mem) for label, timings, mem in meta_runs}
    if "streaming/gz" in meta_by_label and "merge/gz" in meta_by_label:
        st_t = meta_by_label["streaming/gz"][0]
        mg_t = meta_by_label["merge/gz"][0]
        st_total = sum(v for v in st_t.values() if not math.isnan(v))
        mg_total = sum(v for v in mg_t.values() if not math.isnan(v))
        ratio = st_total / mg_total if mg_total > 0 else float("inf")
        print(f"\n  Meta ratio (gz): streaming/merge = {ratio:.2f}x")

    # Conversion section
    if conv_results:
        print("\n  TSV → PARQUET CONVERSION")
        print(f"  {'Method':<28} {'Time (s)':<12} {'Mem (MB)':<12} {'Rows':<12}")
        print(f"  {'-' * 64}")
        for r in conv_results:
            t_s = f"{r['elapsed']:.2f}" if r["elapsed"] == r["elapsed"] else "FAILED"
            m_s = f"{r['peak_mem_mb']:.0f}" if r["peak_mem_mb"] == r["peak_mem_mb"] else "—"
            print(f"  {r['desc']:<28} {t_s:<12} {m_s:<12} {r['n_rows']:<12,}")
        # Speedup vs Python Polars baseline
        py_r = next((r for r in conv_results if r["label"] == "python_polars"), None)
        if py_r and py_r["elapsed"] == py_r["elapsed"]:
            print()
            for r in conv_results:
                if r["label"] == "python_polars" or r["elapsed"] != r["elapsed"]:
                    continue
                speedup = py_r["elapsed"] / r["elapsed"] if r["elapsed"] > 0 else float("inf")
                print(f"  {r['desc']}: {speedup:.2f}x vs Python Polars")

    # ACAT section — two paths: rust_sweep (parquet via convert) vs rust_tsv (TSV direct)
    print("\n  ACAT PIPELINE")
    print(
        f"  {'Path':<16} {'Total (s)':<10} {'Load (s)':<10} "
        f"{'Index (s)':<10} {'Compute (s)':<12} {'Mem (MB)':<10} {'Rows':<10}"
    )
    print(f"  {'-' * 78}")
    for r in acat_results:
        print(
            f"  {r['backend']:<16} {r['total_s']:<10.2f} "
            f"{r['load_s']:<10.2f} {r['index_s']:<10.2f} {r['compute_s']:<12.2f} "
            f"{r['peak_mem_mb']:<10.0f} {r['n_results']:<10,}"
        )

    # Side-by-side speedup: rust_sweep (parquet path) vs rust_tsv (TSV direct)
    by_label = {r["backend"]: r for r in acat_results}
    tsv_r = by_label.get("rust_tsv")
    swp_r = by_label.get("rust_sweep")
    if tsv_r and swp_r:
        speedup = tsv_r["total_s"] / swp_r["total_s"] if swp_r["total_s"] > 0 else float("inf")
        mem_red = (1 - swp_r["peak_mem_mb"] / tsv_r["peak_mem_mb"]) * 100 if tsv_r["peak_mem_mb"] > 0 else 0
        print("\n  ACAT SPEEDUP  rust_sweep (parquet) vs rust_tsv (TSV direct)")
        print(f"  {'Metric':<20} {'rust_tsv':<14} {'rust_sweep':<14} {'Ratio':<10}")
        print(f"  {'-' * 58}")
        print(f"  {'Total (s)':<20} {tsv_r['total_s']:<14.2f} {swp_r['total_s']:<14.2f} {speedup:.2f}x")
        print(
            f"  {'Load (s)':<20} {tsv_r['load_s']:<14.2f} {swp_r['load_s']:<14.2f} "
            f"{tsv_r['load_s'] / swp_r['load_s'] if swp_r['load_s'] > 0 else float('inf'):.2f}x"
        )
        print(
            f"  {'Compute (s)':<20} {tsv_r['compute_s']:<14.2f} {swp_r['compute_s']:<14.2f} "
            f"{tsv_r['compute_s'] / swp_r['compute_s'] if swp_r['compute_s'] > 0 else float('inf'):.2f}x"
        )
        print(
            f"  {'Peak mem (MB)':<20} {tsv_r['peak_mem_mb']:<14.0f} {swp_r['peak_mem_mb']:<14.0f} {mem_red:.1f}% less"
        )

    # Detailed ACAT internal breakdown per path
    print("\n  ACAT INTERNAL BREAKDOWN (per-function)")
    for r in acat_results:
        if not r["internal"]:
            continue
        print(f"\n  [{r['backend']}]")
        print(f"  {'Function':<40} {'Total (s)':<12} {'Calls':<8} {'Mean (ms)':<10}")
        print(f"  {'-' * 70}")
        for fn in sorted(r["internal"], key=lambda x: x["total_s"], reverse=True):
            print(f"  {fn['function']:<40} {fn['total_s']:<12.3f} {fn['calls']:<8} {fn['mean_ms']:<10.2f}")

    print(f"\n{'=' * 80}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    """Run the benchmark."""
    parser = argparse.ArgumentParser(
        description="Benchmark meta-analysis -> ACAT pipeline",
    )
    parser.add_argument(
        "--regenerate",
        action="store_true",
        help="Force regeneration of benchmark data",
    )
    parser.add_argument(
        "--skip-datagen",
        action="store_true",
        help="Skip data generation (assume files exist)",
    )
    args = parser.parse_args()

    if args.skip_datagen:
        if not BENCH_DIR.exists():
            print(f"ERROR: {BENCH_DIR} does not exist. Run without --skip-datagen first.")
            sys.exit(1)
        paths = _collect_paths(BENCH_DIR)
    else:
        paths = generate_data(force=args.regenerate)

    run_benchmark(paths)


if __name__ == "__main__":
    main()
