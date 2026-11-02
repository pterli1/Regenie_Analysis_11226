#!/usr/bin/env python3
"""Compare manc-cojo's back-calculated sample sizes across the SNPs of a conditioning step.

GCTA effect mode replaces each SNP's N with N_eff = (Vp/V - b^2)/se^2 + 1, where Vp is the
median of V*N*(se^2 + b^2/(N-1)) over the .ma, and weights the adjustment of a tested SNP t
for a conditioning SNP c by min(N_c, N_t)/N_t. Every SNP in a job comes from the same
samples, so imputeNA mode (which uses no N) gives the same conditional effects when N_eff is
flat, and GCTA mode weakens the conditioning wherever conditioning SNPs have the lower N_eff.

Reads <job>_{mac10,burden}_mancojo_input.ma.gz and <job>.conditional_list from a step folder.
V is taken as 2pq; with --var-from-ld manc-cojo uses the genotype variance, which is close.

  python3 mancojo_effective_n.py <step folder> [--match eur_chr4] [--tsv rows.tsv]

<step folder> is gs://..., a local folder, or a path under $WORKSPACE_BUCKET, e.g.
WGS_Analysis/<pheno>_binary_wgs/meta_analysis/<pheno>_aggregate_conditioning_<suffix>/agg_cond_step11
"""

from __future__ import annotations

import argparse
import bisect
import contextlib
import gzip
import os
import re
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

SUFFIX = re.compile(r"_(mac10|burden)_mancojo_input\.ma\.gz$")
VARIANT_ID = re.compile(r"^[0-9XYMT]+:\d+:[^:]+:[^:]+$")  # everything else is a mask
QUANTILES = (0.05, 0.25, 0.5, 0.75, 0.95)


def fetch(folder: str, match: str, workdir: Path) -> list[Path]:
    """Return local paths of a step's `.ma.gz` files and their conditional lists, downloading from GCS if needed.

    Args:
        folder: Step folder: `gs://...`, local, or relative to `$WORKSPACE_BUCKET`.
        match: Substring the file names must contain.
        workdir: Local directory to download to.
    """
    if not folder.startswith("gs://") and not Path(folder).is_dir():
        bucket = os.environ.get("WORKSPACE_BUCKET")
        if not bucket:
            sys.exit(f"{folder} is not a local folder and WORKSPACE_BUCKET is not set")
        folder = f"{bucket.rstrip('/')}/{folder.lstrip('/')}"
    if not folder.startswith("gs://"):
        return sorted(p for p in Path(folder).glob("*_mancojo_input.ma.gz") if match in p.name)

    listing = subprocess.run(
        ["gcloud", "storage", "ls", f"{folder.rstrip('/')}/*_mancojo_input.ma.gz"],
        capture_output=True,
        text=True,
        check=False,  # no match is reported below
    )
    remote = [u for u in listing.stdout.split() if match in u.rsplit("/", 1)[-1]]
    if not remote:
        sys.exit(f"No *_mancojo_input.ma.gz matching '{match}' in {folder}\n{listing.stderr}")
    subprocess.run(["gcloud", "storage", "cp", *remote, str(workdir)], check=True)
    # A job without a conditional list is reported with 0 conditioning SNPs
    lists = sorted({f"{SUFFIX.sub('', u)}.conditional_list" for u in remote})
    subprocess.run(["gcloud", "storage", "cp", *lists, str(workdir)], check=False)
    return sorted(workdir / u.rsplit("/", 1)[-1] for u in remote)


def quantiles(values: list[float]) -> str:
    """Format the `QUANTILES` of values, or `-` if there are none."""
    if not values:
        return "-"
    s = sorted(values)
    return "  ".join(f"{s[min(len(s) - 1, int(q * len(s)))]:.3f}" for q in QUANTILES)


def analyse(ma: Path, rows_out) -> None:
    """Print N_eff/N quantiles of one job's `.ma.gz` by SNP group, and the GCTA weights they imply.

    Args:
        ma: manc-cojo input `.ma.gz`.
        rows_out: Open file to write per-SNP rows to, or None.
    """
    cond_path = Path(SUFFIX.sub("", str(ma)) + ".conditional_list")
    cond = set(cond_path.read_text().split()) if cond_path.exists() else set()

    rows, vp_terms = [], []
    with gzip.open(ma, "rt") as fh:
        next(fh)
        for line in fh:
            snp, _a1, _a2, freq, b, se, _p, n = line.split()
            freq, b, se, n = float(freq), float(b), float(se), float(n)
            v = 2 * freq * (1 - freq)
            vp_terms.append(v * n * (se * se + b * b / (n - 1)))
            if v > 0 and se > 0:
                rows.append((snp, freq, b, se, n, v))
    vp = statistics.median(vp_terms)

    groups: dict[str, list[float]] = {}
    neff_by_group: dict[str, list[float]] = {}
    for snp, freq, b, se, n, v in rows:
        neff = (vp / v - b * b) / (se * se) + 1
        if snp in cond:
            group = "conditioning variants"
        elif not VARIANT_ID.match(snp):
            group = "tested masks"
        else:
            group = "other variants, MAC " + ("<= 10" if 2 * min(freq, 1 - freq) * n <= 10.5 else ">= 11")
        groups.setdefault(group, []).append(neff / n)
        neff_by_group.setdefault(group, []).append(neff)
        if rows_out:
            rows_out.write(f"{ma.name}\t{snp}\t{group}\t{freq:.6g}\t{n:.0f}\t{neff:.6g}\t{neff / n:.4f}\n")

    print(f"\n{ma.name}\n  Vp {vp:.4g}; {len(rows)} SNPs; {len(cond)} in conditional list")
    print(f"  N_eff / N quantiles ({', '.join(f'{q:.0%}' for q in QUANTILES)}):")
    for group in sorted(groups):
        print(f"    {group:<30} n={len(groups[group]):<6} {quantiles(groups[group])}")

    # GCTA weights the adjustment for conditioning SNP c by min(N_c, N_t) / N_t; imputeNA by 1
    cond_neff = sorted(neff_by_group.get("conditioning variants", []))
    masks = neff_by_group.get("tested masks", [])
    if cond_neff and masks:
        median_cond = statistics.median(cond_neff)
        typical = [min(1.0, median_cond / nt) for nt in masks]
        weak = [bisect.bisect_left(cond_neff, 0.9 * nt) / len(cond_neff) for nt in masks]
        print("  GCTA weighting of conditioning adjustments per tested mask (imputeNA = 1):")
        print(f"    for the median conditioning variant      {quantiles(typical)}")
        print(f"    share of conditioning variants below 0.9 {quantiles(weak)}")


def main() -> None:
    """Parse the arguments and analyse each matching file."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", help="conditioning step folder: gs://..., local, or relative to $WORKSPACE_BUCKET")
    ap.add_argument("--match", default="", help="only files whose name contains this, e.g. eur_chr4")
    ap.add_argument("--tsv", help="also write per-SNP values to this file")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        files = fetch(args.folder, args.match, Path(tmp))
        if not files:
            sys.exit("No *_mancojo_input.ma.gz files found")
        with open(args.tsv, "w") if args.tsv else contextlib.nullcontext() as rows_out:
            if rows_out:
                rows_out.write("file\tSNP\tgroup\tfreq\tN\tN_eff\tN_eff_over_N\n")
            for ma in files:
                analyse(ma, rows_out)
    print(
        "\nIf conditioning variants and tested masks have similar N_eff/N, and the weights are near 1,\n"
        "GCTA and imputeNA give nearly the same conditional effects. Weights well below 1 mean GCTA\n"
        "applies only that fraction of the conditioning, leaving residual signal in those masks."
    )


if __name__ == "__main__":
    main()
