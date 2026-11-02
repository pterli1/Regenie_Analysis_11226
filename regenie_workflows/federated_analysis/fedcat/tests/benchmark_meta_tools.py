#!/usr/bin/env python3
"""Benchmark fedcat meta-analysis against GWAMA and METAL.

This script generates synthetic per-cohort summary-stat files, then runs and
times three meta-analysis engines on the same data:

1. fedcat (native Rust-backed meta-analysis)
2. GWAMA (-qt mode)
3. METAL (SCHEME STDERR)

It also records peak RSS from ``/usr/bin/time -v``. The benchmark is intended
for Linux environments where GWAMA and METAL binaries are available.

Usage:

    python fedcat/tests/benchmark_meta_tools.py
    python fedcat/tests/benchmark_meta_tools.py --regenerate
    python fedcat/tests/benchmark_meta_tools.py --n-variants 500000 --n-cohorts 12

Notes:
- fedcat is benchmarked through a worker mode of this same script so all tools
  are measured with the same timing/memory method.
- GWAMA and METAL executable names can be overridden with CLI flags.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
_REGENIE_WORKFLOWS = Path(__file__).resolve().parents[2]
if str(_REGENIE_WORKFLOWS) not in sys.path:
    sys.path.insert(0, str(_REGENIE_WORKFLOWS))

os.environ.setdefault("WORKSPACE_BUCKET", "gs://dummy-bucket")

from fedcat.stats.meta_analysis import meta_analyse_regenie, z_to_log10p

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_BENCH_DIR = Path("./bench_outputs")
DEFAULT_N_VARIANTS = 200_000
DEFAULT_N_COHORTS = 10
DEFAULT_SEED = 42

FEDCAT_SHARED_COLUMN_CONFIG = {
    "chrom": "CHROM",
    "id": "SNP",
    "allele0": "NEA",
    "allele1": "EA",
    "a1freq": "EAF",
    "n": "N",
    "test": "TEST",
    "beta": "BETA",
    "se": "SE",
    "log10p": "LOG10P",
    "separator": " ",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@dataclass
class ToolRun:
    """Outcome and resource use of one tool's run."""

    tool: str
    ok: bool
    elapsed_s: float
    peak_rss_mb: float
    cpu_user_s: float
    cpu_system_s: float
    cpu_total_s: float
    output_rows: int
    output_path: str
    note: str = ""
    stdout_log: str = ""
    stderr_log: str = ""
    command_log: str = ""


def _format_elapsed_hms(seconds: float) -> str:
    """Format seconds as `h:mm:ss`, or `NA` if not finite."""
    if not math.isfinite(seconds):
        return "NA"
    m, s = divmod(round(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h:d}:{m:02d}:{s:02d}"


def _parse_time_v_stderr(stderr_text: str) -> tuple[float, float, float, float]:
    """Parse GNU ``/usr/bin/time -v`` output.

    Returns ``(elapsed_seconds, peak_rss_mb, user_cpu_s, system_cpu_s)``.
    """
    elapsed_s = float("nan")
    peak_mb = float("nan")
    user_s = float("nan")
    system_s = float("nan")

    rss_match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr_text)
    if rss_match:
        peak_mb = int(rss_match.group(1)) / 1024.0

    user_match = re.search(r"User time \(seconds\):\s*([0-9.]+)", stderr_text)
    if user_match:
        user_s = float(user_match.group(1))

    system_match = re.search(r"System time \(seconds\):\s*([0-9.]+)", stderr_text)
    if system_match:
        system_s = float(system_match.group(1))

    elapsed_match = re.search(r"Elapsed \(wall clock\) time .*?:\s*([0-9:.]+)", stderr_text)
    if elapsed_match:
        raw = elapsed_match.group(1).strip()
        parts = raw.split(":")
        try:
            if len(parts) == 3:
                h = float(parts[0])
                m = float(parts[1])
                s = float(parts[2])
                elapsed_s = h * 3600 + m * 60 + s
            elif len(parts) == 2:
                m = float(parts[0])
                s = float(parts[1])
                elapsed_s = m * 60 + s
            elif len(parts) == 1:
                elapsed_s = float(parts[0])
        except ValueError:
            pass

    return elapsed_s, peak_mb, user_s, system_s


def _run_with_time_v(
    command: list[str], cwd: Path | None = None
) -> tuple[int, float, float, float, float, str, str]:
    """Run command with GNU ``time -v`` and parse wall-time/RSS.

    Returns ``(exit_code, elapsed_s, peak_rss_mb, user_cpu_s, system_cpu_s, stdout, stderr)``.
    """
    time_bin = shutil.which("time") or "/usr/bin/time"
    wrapped = [time_bin, "-v", *command]
    t0 = time.perf_counter()

    proc = subprocess.Popen(
        wrapped,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    assert proc.stdout is not None
    assert proc.stderr is not None
    stdout_chunks: list[str] = []
    stderr_chunks: list[str] = []

    def _pump(stream, sink, chunks: list[str]) -> None:
        for line in iter(stream.readline, ""):
            sink.write(line)
            sink.flush()
            chunks.append(line)
        stream.close()

    t_out = threading.Thread(target=_pump, args=(proc.stdout, sys.stdout, stdout_chunks), daemon=True)
    t_err = threading.Thread(target=_pump, args=(proc.stderr, sys.stderr, stderr_chunks), daemon=True)
    t_out.start()
    t_err.start()
    rc = proc.wait()
    t_out.join()
    t_err.join()

    stdout_text = "".join(stdout_chunks)
    stderr_text = "".join(stderr_chunks)

    elapsed_s, peak_mb, user_s, system_s = _parse_time_v_stderr(stderr_text)
    if not math.isfinite(elapsed_s):
        elapsed_s = time.perf_counter() - t0
    return rc, elapsed_s, peak_mb, user_s, system_s, stdout_text, stderr_text


def _tail_lines(text: str, n: int = 12) -> str:
    """Return the last `n` non-empty lines of `text`."""
    lines = [ln for ln in text.strip().splitlines() if ln.strip()]
    if not lines:
        return ""
    return "\n".join(lines[-n:])


def _write_tool_logs(
    logs_dir: Path,
    tag: str,
    command: list[str],
    stdout: str,
    stderr: str,
) -> tuple[Path, Path, Path]:
    """Write a tool run's command, stdout and stderr to `logs_dir`, returning the three paths."""
    logs_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = logs_dir / f"{tag}.stdout.log"
    stderr_path = logs_dir / f"{tag}.stderr.log"
    command_path = logs_dir / f"{tag}.command.txt"

    with open(stdout_path, "w") as fh:
        fh.write(stdout)
    with open(stderr_path, "w") as fh:
        fh.write(stderr)
    with open(command_path, "w") as fh:
        fh.write(" ".join(command) + "\n")

    return stdout_path, stderr_path, command_path


def _abspath_from(base_dir: Path, p: str) -> str:
    """Return `p` as an absolute path, relative to `base_dir` if not already absolute."""
    path = Path(p)
    if path.is_absolute():
        return str(path)
    return str((base_dir / path).resolve())


def _normalize_manifest_paths(manifest: dict, base_dir: Path) -> dict:
    """Return a copy of manifest with key file paths converted to absolute.

    Relative paths are interpreted relative to *base_dir* (the launch CWD).
    """
    out = dict(manifest)

    if "bench_dir" in out:
        out["bench_dir"] = _abspath_from(base_dir, out["bench_dir"])
    if "gwama_in" in out:
        out["gwama_in"] = _abspath_from(base_dir, out["gwama_in"])
    if "metal_script" in out:
        out["metal_script"] = _abspath_from(base_dir, out["metal_script"])
    if "fedcat_output" in out:
        out["fedcat_output"] = _abspath_from(base_dir, out["fedcat_output"])
    if "gwama_output_root" in out:
        out["gwama_output_root"] = _abspath_from(base_dir, out["gwama_output_root"])
    if "metal_output_root" in out:
        out["metal_output_root"] = _abspath_from(base_dir, out["metal_output_root"])

    if "fedcat_files" in out:
        fixed = {}
        for cohort, paths in out["fedcat_files"].items():
            fixed[cohort] = [_abspath_from(base_dir, p) for p in paths]
        out["fedcat_files"] = fixed

    return out


def _refresh_external_tool_inputs(manifest: dict) -> None:
    """Rewrite GWAMA and METAL control files using absolute cohort file paths.

    This prevents failures when cached benchmark folders were generated with
    older relative-path manifests.
    """
    cohorts = manifest.get("cohorts", [])
    fedcat_files = manifest.get("fedcat_files", {})

    shared_files: list[str] = []
    for cohort in cohorts:
        cohort_paths = fedcat_files.get(cohort, [])
        if not cohort_paths:
            raise FileNotFoundError(f"Missing cohort path entry in manifest for: {cohort}")
        p = Path(cohort_paths[0])
        if not p.exists():
            raise FileNotFoundError(f"Missing shared cohort file: {p}")
        shared_files.append(str(p.resolve()))

    gwama_in = Path(manifest["gwama_in"])
    gwama_in.parent.mkdir(parents=True, exist_ok=True)
    with open(gwama_in, "w") as fh:
        for p in shared_files:
            fh.write(f"{p}\n")

    metal_script = Path(manifest["metal_script"])
    metal_script.parent.mkdir(parents=True, exist_ok=True)
    metal_out_root = Path(manifest["metal_output_root"])
    with open(metal_script, "w") as fh:
        fh.write("SCHEME STDERR\n")
        fh.write("SEPARATOR WHITESPACE\n")
        fh.write("MARKER SNP\n")
        fh.write("ALLELE EA NEA\n")
        fh.write("EFFECT BETA\n")
        fh.write("STDERR SE\n")
        fh.write("WEIGHT N\n")
        for p in shared_files:
            fh.write(f"PROCESS {p}\n")
        fh.write(f"OUTFILE {metal_out_root} .tbl\n")
        fh.write("ANALYZE\n")
        fh.write("QUIT\n")


def _chrom_rank_py(chrom: str) -> int:
    """Sort rank of a chromosome: 1-22, X, Y, MT, then others."""
    c = chrom.upper()
    if c == "X" or c == "23":
        return 23
    if c == "Y":
        return 24
    if c == "MT":
        return 25
    try:
        return int(c)
    except ValueError:
        return 26


def _count_data_rows(path: Path) -> int:
    """Return the rows of a file after its header, or 0 if it is missing."""
    if not path.exists():
        return 0
    with open(path, "r") as fh:
        n = sum(1 for _ in fh)
    return max(0, n - 1)


def _make_variant_ids(chrom: np.ndarray, pos: np.ndarray, a0: np.ndarray, a1: np.ndarray) -> np.ndarray:
    """Return `chr:pos:ref:alt` IDs."""
    return np.array([f"{c}:{p}:{r}:{a}" for c, p, r, a in zip(chrom, pos, a0, a1)], dtype=object)


def _write_shared_txt(path: Path, rows: dict[str, np.ndarray]) -> None:
    """Write a cohort's results in the space-separated input format GWAMA and METAL read."""
    with open(path, "w") as fh:
        fh.write("SNP CHROM GENPOS EA NEA EAF N BETA SE LOG10P TEST\n")
        fh.writelines(
            f"{rows['ID'][i]} {rows['CHROM'][i]} {int(rows['GENPOS'][i])} "
            f"{rows['ALLELE1'][i]} {rows['ALLELE0'][i]} {rows['A1FREQ'][i]:.6f} "
            f"{int(rows['N'][i])} {rows['BETA'][i]:.6f} {rows['SE'][i]:.6f} "
            f"{rows['LOG10P'][i]:.6f} ADD\n"
            for i in range(len(rows["ID"]))
        )


def generate_data(bench_dir: Path, n_variants: int, n_cohorts: int, seed: int, force: bool = False) -> dict:
    """Generate synthetic inputs shared by fedcat, GWAMA, and METAL."""
    bench_dir = bench_dir.resolve()
    marker = bench_dir / "generated_manifest.json"
    if marker.exists() and not force:
        with open(marker) as fh:
            return json.load(fh)

    rng = np.random.default_rng(seed)

    (bench_dir / "inputs" / "shared").mkdir(parents=True, exist_ok=True)
    (bench_dir / "inputs" / "gwama").mkdir(parents=True, exist_ok=True)
    (bench_dir / "inputs" / "metal").mkdir(parents=True, exist_ok=True)
    (bench_dir / "outputs").mkdir(parents=True, exist_ok=True)

    # Shared base marker universe.
    chrom = rng.integers(1, 23, size=n_variants).astype(str)
    pos = rng.integers(1, 250_000_000, size=n_variants)
    bases = np.array(["A", "C", "G", "T"], dtype=object)
    ref_i = rng.integers(0, 4, size=n_variants)
    alt_i = (ref_i + rng.integers(1, 4, size=n_variants)) % 4
    a0 = bases[ref_i]
    a1 = bases[alt_i]
    marker_id = _make_variant_ids(chrom, pos, a0, a1)

    # Force unique IDs to avoid collisions.
    _, uniq_idx = np.unique(marker_id, return_index=True)
    keep = np.zeros(n_variants, dtype=bool)
    keep[uniq_idx] = True
    chrom = chrom[keep]
    pos = pos[keep]
    a0 = a0[keep]
    a1 = a1[keep]
    marker_id = marker_id[keep]

    # Sort once globally so each cohort subset preserves merge-compatible order.
    chrom_rank = np.fromiter((_chrom_rank_py(c) for c in chrom), dtype=np.int32, count=chrom.size)
    sort_idx = np.lexsort((marker_id, pos, chrom_rank))
    chrom = chrom[sort_idx]
    pos = pos[sort_idx]
    a0 = a0[sort_idx]
    a1 = a1[sort_idx]
    marker_id = marker_id[sort_idx]

    n_effective = marker_id.size

    true_beta = rng.normal(0.0, 0.03, size=n_effective)
    base_se = rng.uniform(0.01, 0.06, size=n_effective)
    base_af = rng.uniform(0.0005, 0.5, size=n_effective)

    cohort_names = [f"cohort_{i + 1:02d}" for i in range(n_cohorts)]
    fedcat_files: dict[str, list[str]] = {}
    shared_files: list[str] = []

    print(f"Generating synthetic data: {n_effective:,} unique markers x {n_cohorts} cohorts")
    t0 = time.perf_counter()

    for ci, cohort in enumerate(cohort_names):
        cohort_seed = seed + 10_000 + ci
        crng = np.random.default_rng(cohort_seed)

        # Cohort-level missingness and perturbation around shared truth.
        miss = crng.uniform(0.05, 0.15)
        keep_mask = crng.random(n_effective) >= miss

        c_beta = true_beta[keep_mask] + crng.normal(0.0, 0.02, size=keep_mask.sum())
        c_se = np.clip(base_se[keep_mask] * crng.uniform(0.9, 1.1, size=keep_mask.sum()), 0.001, None)
        z = c_beta / c_se
        c_log10p = z_to_log10p(z)
        c_chisq = z * z
        c_af = np.clip(base_af[keep_mask] + crng.normal(0.0, 0.02, size=keep_mask.sum()), 1e-6, 0.999999)
        c_n = crng.integers(35_000, 120_000, size=keep_mask.sum())

        rows = {
            "CHROM": chrom[keep_mask],
            "GENPOS": pos[keep_mask],
            "ID": marker_id[keep_mask],
            "ALLELE0": a0[keep_mask],
            "ALLELE1": a1[keep_mask],
            "A1FREQ": c_af,
            "N": c_n,
            "BETA": c_beta,
            "SE": c_se,
            "CHISQ": c_chisq,
            "LOG10P": c_log10p,
        }

        shared_path = bench_dir / "inputs" / "shared" / f"{cohort}.txt"
        _write_shared_txt(shared_path, rows)

        fedcat_files[cohort] = [str(shared_path)]
        shared_files.append(str(shared_path))

        if (ci + 1) % 5 == 0 or (ci + 1) == n_cohorts:
            print(f"  wrote {ci + 1}/{n_cohorts} cohorts")

    gwama_in = bench_dir / "inputs" / "gwama" / "gwama.in"
    with open(gwama_in, "w") as fh:
        fh.writelines(f"{p}\n" for p in shared_files)

    # METAL script (STDERR scheme using BETA + SE).
    metal_script = bench_dir / "inputs" / "metal" / "metal_script.txt"
    metal_out_root = bench_dir / "outputs" / "metal_meta_"
    with open(metal_script, "w") as fh:
        fh.write("SCHEME STDERR\n")
        fh.write("SEPARATOR WHITESPACE\n")
        fh.write("MARKER SNP\n")
        fh.write("ALLELE EA NEA\n")
        fh.write("EFFECT BETA\n")
        fh.write("STDERR SE\n")
        fh.write("WEIGHT N\n")
        fh.writelines(f"PROCESS {p}\n" for p in shared_files)
        fh.write(f"OUTFILE {metal_out_root} .tbl\n")
        fh.write("ANALYZE\n")
        fh.write("QUIT\n")

    manifest = {
        "bench_dir": str(bench_dir),
        "seed": seed,
        "n_variants_requested": int(n_variants),
        "n_variants_effective": int(n_effective),
        "n_cohorts": int(n_cohorts),
        "cohorts": cohort_names,
        "fedcat_files": fedcat_files,
        "fedcat_column_config": FEDCAT_SHARED_COLUMN_CONFIG,
        "gwama_in": str(gwama_in),
        "metal_script": str(metal_script),
        "fedcat_output": str(bench_dir / "outputs" / "fedcat_META_single_variants.tsv"),
        "gwama_output_root": str(bench_dir / "outputs" / "gwama_meta"),
        "metal_output_root": str(metal_out_root),
    }

    with open(marker, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"Data generation finished in {time.perf_counter() - t0:.1f}s")
    return manifest


def _worker_fedcat(manifest_path: Path, backend: str, output_override: Path | None = None) -> int:
    """Run fedcat's meta-analysis on a manifest's inputs, as a worker process to measure."""
    with open(manifest_path) as fh:
        manifest = json.load(fh)
    per_cohort_files = {k: [Path(p) for p in v] for k, v in manifest["fedcat_files"].items()}
    output_file = output_override if output_override is not None else Path(manifest["fedcat_output"])
    output_file.parent.mkdir(parents=True, exist_ok=True)

    meta_analyse_regenie(
        per_cohort_files=per_cohort_files,
        output_file=output_file,
        res_type="single_variants",
        column_config=manifest.get("fedcat_column_config", None),
        backend=backend,
    )
    return 0


def _find_latest_metal_table(out_root: str) -> Path | None:
    # METAL writes OUTFILE root + index + suffix, e.g. metal_meta_1.tbl
    """Return METAL's newest output table for an output root, or None."""
    root = Path(out_root)
    parent = root.parent
    prefix = root.name
    candidates = sorted(parent.glob(f"{prefix}*.tbl"), key=lambda p: p.stat().st_mtime)
    return candidates[-1] if candidates else None


def run_benchmark(
    manifest_path: Path,
    gwama_bin: str,
    metal_bin: str,
    fedcat_backends: list[str],
    require_all: bool,
    logs_dir: Path | None = None,
) -> list[ToolRun]:
    """Run fedcat (each backend), GWAMA and METAL on the manifest's inputs, measuring each."""
    launch_cwd = Path.cwd().resolve()
    with open(manifest_path) as fh:
        manifest = _normalize_manifest_paths(json.load(fh), launch_cwd)

    # Enforce known-good fedcat mapping even with cached old manifests.
    manifest["fedcat_column_config"] = dict(FEDCAT_SHARED_COLUMN_CONFIG)

    bench_dir = Path(manifest["bench_dir"])
    out_dir = bench_dir / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = logs_dir or (out_dir / "logs")
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Persist normalized manifest so worker subprocess can read stable absolute paths.
    manifest_abs_path = out_dir / "generated_manifest.abs.json"
    with open(manifest_abs_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    # Refresh GWAMA and METAL control files from normalized absolute paths.
    _refresh_external_tool_inputs(manifest)

    results: list[ToolRun] = []
    py = sys.executable
    this_script = Path(__file__).resolve()

    for fedcat_backend in fedcat_backends:
        print(f"\nRunning fedcat benchmark ({fedcat_backend})...")
        fedcat_out = out_dir / f"fedcat_META_single_variants_{fedcat_backend}.tsv"
        fedcat_cmd = [
            py,
            str(this_script),
            "--worker",
            "fedcat",
            "--manifest",
            str(manifest_abs_path),
            "--fedcat-backend",
            fedcat_backend,
            "--fedcat-output",
            str(fedcat_out),
        ]
        rc, elapsed, peak_mb, user_s, system_s, stdout, stderr = _run_with_time_v(fedcat_cmd, cwd=bench_dir)
        f_stdout, f_stderr, f_cmd = _write_tool_logs(
            logs_dir,
            f"fedcat_{fedcat_backend}",
            fedcat_cmd,
            stdout,
            stderr,
        )
        fedcat_rows = _count_data_rows(fedcat_out)
        results.append(
            ToolRun(
                tool=f"fedcat/{fedcat_backend}",
                ok=(rc == 0 and fedcat_out.exists()),
                elapsed_s=elapsed,
                peak_rss_mb=peak_mb,
                cpu_user_s=user_s,
                cpu_system_s=system_s,
                cpu_total_s=(
                    (user_s + system_s) if math.isfinite(user_s) and math.isfinite(system_s) else float("nan")
                ),
                output_rows=fedcat_rows,
                output_path=str(fedcat_out),
                note="" if rc == 0 else _tail_lines(stderr),
                stdout_log=str(f_stdout),
                stderr_log=str(f_stderr),
                command_log=str(f_cmd),
            )
        )

    gwama_path = shutil.which(gwama_bin)
    if gwama_path:
        print("Running GWAMA benchmark...")
        gwama_path = str(Path(gwama_path).resolve())
        gwama_out_root = Path(manifest["gwama_output_root"])
        gwama_cmd = [
            gwama_path,
            "-qt",
            "-i",
            manifest["gwama_in"],
            "-o",
            str(gwama_out_root),
            "--name_marker",
            "SNP",
            "--name_ea",
            "EA",
            "--name_nea",
            "NEA",
            "--name_beta",
            "BETA",
            "--name_se",
            "SE",
            "--name_n",
            "N",
        ]
        rc, elapsed, peak_mb, user_s, system_s, stdout, stderr = _run_with_time_v(gwama_cmd, cwd=bench_dir)
        g_stdout, g_stderr, g_cmd = _write_tool_logs(logs_dir, "gwama", gwama_cmd, stdout, stderr)
        gwama_out = Path(str(gwama_out_root) + ".out")
        rows = _count_data_rows(gwama_out)
        results.append(
            ToolRun(
                tool="GWAMA",
                ok=(rc == 0 and gwama_out.exists()),
                elapsed_s=elapsed,
                peak_rss_mb=peak_mb,
                cpu_user_s=user_s,
                cpu_system_s=system_s,
                cpu_total_s=(
                    (user_s + system_s) if math.isfinite(user_s) and math.isfinite(system_s) else float("nan")
                ),
                output_rows=rows,
                output_path=str(gwama_out),
                note="" if rc == 0 else _tail_lines(stderr),
                stdout_log=str(g_stdout),
                stderr_log=str(g_stderr),
                command_log=str(g_cmd),
            )
        )
    else:
        msg = f"GWAMA binary not found: {gwama_bin}"
        if require_all:
            raise RuntimeError(msg)
        results.append(
            ToolRun(
                tool="GWAMA",
                ok=False,
                elapsed_s=float("nan"),
                peak_rss_mb=float("nan"),
                cpu_user_s=float("nan"),
                cpu_system_s=float("nan"),
                cpu_total_s=float("nan"),
                output_rows=0,
                output_path="",
                note=msg,
                stdout_log="",
                stderr_log="",
                command_log="",
            )
        )

    metal_path = shutil.which(metal_bin)
    if metal_path:
        print("Running METAL benchmark...")
        metal_path = str(Path(metal_path).resolve())
        metal_cmd = [metal_path, manifest["metal_script"]]
        rc, elapsed, peak_mb, user_s, system_s, stdout, stderr = _run_with_time_v(metal_cmd, cwd=bench_dir)
        m_stdout, m_stderr, m_cmd = _write_tool_logs(logs_dir, "metal", metal_cmd, stdout, stderr)
        metal_out = _find_latest_metal_table(manifest["metal_output_root"])
        rows = _count_data_rows(metal_out) if metal_out else 0
        results.append(
            ToolRun(
                tool="METAL",
                ok=(rc == 0 and metal_out is not None and metal_out.exists()),
                elapsed_s=elapsed,
                peak_rss_mb=peak_mb,
                cpu_user_s=user_s,
                cpu_system_s=system_s,
                cpu_total_s=(
                    (user_s + system_s) if math.isfinite(user_s) and math.isfinite(system_s) else float("nan")
                ),
                output_rows=rows,
                output_path=str(metal_out) if metal_out else "",
                note="" if rc == 0 else _tail_lines(stderr),
                stdout_log=str(m_stdout),
                stderr_log=str(m_stderr),
                command_log=str(m_cmd),
            )
        )
    else:
        msg = f"METAL binary not found: {metal_bin}"
        if require_all:
            raise RuntimeError(msg)
        results.append(
            ToolRun(
                tool="METAL",
                ok=False,
                elapsed_s=float("nan"),
                peak_rss_mb=float("nan"),
                cpu_user_s=float("nan"),
                cpu_system_s=float("nan"),
                cpu_total_s=float("nan"),
                output_rows=0,
                output_path="",
                note=msg,
                stdout_log="",
                stderr_log="",
                command_log="",
            )
        )

    summary_path = out_dir / "benchmark_summary.json"
    with open(summary_path, "w") as fh:
        json.dump([r.__dict__ for r in results], fh, indent=2)
    print(f"\nSaved machine-readable summary: {summary_path}")
    print(f"Saved run logs: {logs_dir}")

    return results


def print_summary(results: list[ToolRun]) -> None:
    """Print a table of the tools' runs."""
    print("\n" + "=" * 82)
    print("META-ANALYSIS BENCHMARK SUMMARY")
    print("=" * 82)
    print(
        f"{'Tool':<20} {'Status':<8} {'Wall (s)':<10} {'CPU (s)':<10} {'Peak RSS (MB)':<14} {'Rows':<12} {'Output'}"
    )
    print("-" * 82)
    for r in results:
        status = "OK" if r.ok else "FAIL"
        t_str = f"{r.elapsed_s:.2f}" if math.isfinite(r.elapsed_s) else "NA"
        c_str = f"{r.cpu_total_s:.2f}" if math.isfinite(r.cpu_total_s) else "NA"
        m_str = f"{r.peak_rss_mb:.0f}" if math.isfinite(r.peak_rss_mb) else "NA"
        out = r.output_path if r.output_path else "-"
        print(f"{r.tool:<20} {status:<8} {t_str:<10} {c_str:<10} {m_str:<14} {r.output_rows:<12,} {out}")
        if r.note and not r.ok:
            print(f"  note: {r.note.strip().splitlines()[-1]}")
            if r.stderr_log:
                print(f"  stderr log: {r.stderr_log}")
            if r.stdout_log:
                print(f"  stdout log: {r.stdout_log}")

    ok_runs = [r for r in results if r.ok and math.isfinite(r.elapsed_s)]
    if len(ok_runs) >= 2:
        baseline = next((r for r in ok_runs if r.tool.startswith("fedcat/")), ok_runs[0])
        print("\nRelative to baseline:")
        for r in ok_runs:
            ratio = r.elapsed_s / baseline.elapsed_s if baseline.elapsed_s > 0 else float("inf")
            print(f"  {r.tool:<20} {ratio:>8.2f}x")


def parse_args() -> argparse.Namespace:
    """Parse the benchmark's arguments."""
    parser = argparse.ArgumentParser(description="Benchmark fedcat vs GWAMA vs METAL meta-analysis")
    parser.add_argument("--bench-dir", type=Path, default=DEFAULT_BENCH_DIR)
    parser.add_argument("--n-variants", "--n-var", type=int, default=DEFAULT_N_VARIANTS)
    parser.add_argument("--n-cohorts", type=int, default=DEFAULT_N_COHORTS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--regenerate", action="store_true", help="Force re-generation of synthetic inputs")
    parser.add_argument(
        "--skip-datagen", action="store_true", help="Skip data generation and use existing manifest"
    )
    parser.add_argument("--gwama-bin", default="./GWAMA", help="GWAMA executable name/path")
    parser.add_argument("--metal-bin", default="./metal", help="METAL executable name/path")
    parser.add_argument(
        "--logs-dir",
        type=Path,
        default=None,
        help="Directory for per-tool stdout/stderr/command logs (default: <bench-dir>/outputs/logs)",
    )
    parser.add_argument(
        "--fedcat-backend",
        choices=["streaming", "merge", "both"],
        default="both",
        help="fedcat backend mode: one backend or both",
    )
    parser.add_argument(
        "--require-all",
        action="store_true",
        help="Fail immediately if GWAMA or METAL executables are not found",
    )

    # Internal worker mode used so fedcat can be measured by /usr/bin/time too.
    parser.add_argument("--worker", choices=["fedcat"], default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument(
        "--fedcat-output",
        type=Path,
        default=None,
        help="Worker mode only: override fedcat output path for this backend run",
    )

    return parser.parse_args()


def main() -> int:
    """Run the benchmark, or one tool as a worker process."""
    args = parse_args()

    if args.worker:
        if args.manifest is None:
            raise ValueError("--manifest is required in --worker mode")
        match args.worker:
            case "fedcat":
                return _worker_fedcat(
                    args.manifest, backend=args.fedcat_backend, output_override=args.fedcat_output
                )
        return 1

    bench_dir = args.bench_dir
    manifest_path = bench_dir / "generated_manifest.json"

    if args.skip_datagen:
        if not manifest_path.exists():
            raise FileNotFoundError(f"Manifest not found: {manifest_path}")
    else:
        generate_data(
            bench_dir=bench_dir,
            n_variants=args.n_variants,
            n_cohorts=args.n_cohorts,
            seed=args.seed,
            force=args.regenerate,
        )

    fedcat_backends = ["streaming", "merge"] if args.fedcat_backend == "both" else [args.fedcat_backend]

    results = run_benchmark(
        manifest_path=manifest_path,
        gwama_bin=args.gwama_bin,
        metal_bin=args.metal_bin,
        fedcat_backends=fedcat_backends,
        require_all=args.require_all,
        logs_dir=args.logs_dir,
    )
    print_summary(results)

    return 0 if any(r.ok for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
