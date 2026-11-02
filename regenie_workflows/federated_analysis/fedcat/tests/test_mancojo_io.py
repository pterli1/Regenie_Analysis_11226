"""Tests for the Manc-COJO ``.cma.cojo`` reformatter.

The fixtures below are real ``manc_cojo --cojo-cond`` output, not hand-written
approximations of it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from fedcat.core.io import (
    MANCOJO_OUTPUT_SCHEMA,
    reformat_mancojo_output,
)

HEADER = "Chr\tSNP\tbp\tA1\tA2\tfreq\tb\tse\tp\tn\tbC\tbC_se\tpC"


def _write(tmp_path: Path, rows: list[str], name: str = "out.cma.cojo") -> Path:
    """Write a `.cma.cojo` file with the given rows."""
    f = tmp_path / name
    f.write_text(HEADER + "\n" + "\n".join(rows) + "\n")
    return f


def test_reformats_real_mancojo_output(tmp_path):
    """Column mapping against genuine manc_cojo output for a chrX locus."""
    src = _write(
        tmp_path,
        [
            (
                "23\tX:200:C:T\t200\tT\tC\t0.3\t0.05\t0.01\t1e-08\t8836.124108\t"
                "0.03690692659\t0.01001357141\t0.0002280832572"
            )
        ],
    )
    out = tmp_path / "reformatted.tsv"
    reformat_mancojo_output([src], out)
    df = pl.read_csv(out, separator=" ", schema_overrides={"CHROM": pl.Utf8})

    assert df.columns == list(MANCOJO_OUTPUT_SCHEMA)
    row = df.row(0, named=True)
    assert row["CHROM"] == "X"  # 23 -> X, as for REGENIE output
    assert row["ID"] == "X:200:C:T"
    assert row["GENPOS"] == 200
    assert row["ALLELE1"] == "T"  # COJO A1 is the effect allele
    assert row["ALLELE0"] == "C"
    assert row["A1FREQ"] == pytest.approx(0.3)
    assert row["BETA"] == pytest.approx(0.03690692659)  # bC, not the marginal b
    assert row["SE"] == pytest.approx(0.01001357141)
    assert row["N"] == 8836  # fractional effective N rounded
    assert row["TEST"] == "ADD"


def test_log10p_recomputed_not_read_from_pC(tmp_path):
    """LOG10P comes from bC/bC_se, so it survives where pC has underflowed to 0."""
    src = _write(
        tmp_path,
        # |z| ~ 40 -> p far below the ~1e-308 double floor; manc-cojo writes 0
        ["1\t1:100:A:G\t100\tG\tA\t0.25\t0.4\t0.01\t0\t50000\t0.4\t0.01\t0"],
    )
    out = tmp_path / "reformatted.tsv"
    reformat_mancojo_output([src], out)
    df = pl.read_csv(out, separator=" ", schema_overrides={"CHROM": pl.Utf8})

    log10p = df["LOG10P"][0]
    assert np.isfinite(log10p)
    assert log10p > 308, "underflowed pC must not cap the reported significance"


def test_autosome_chrom_untouched_and_multiple_files_concatenated(tmp_path):
    a = _write(tmp_path, ["1\t1:10:A:G\t10\tG\tA\t0.1\t0.2\t0.02\t1e-20\t1000\t0.15\t0.02\t1e-12"], "a.cma.cojo")
    b = _write(tmp_path, ["7\t7:20:C:T\t20\tT\tC\t0.4\t0.3\t0.03\t1e-22\t1000\t0.25\t0.03\t1e-15"], "b.cma.cojo")
    out = tmp_path / "reformatted.tsv"
    reformat_mancojo_output([a, b], out)
    df = pl.read_csv(out, separator=" ", schema_overrides={"CHROM": pl.Utf8})

    assert df.height == 2
    assert set(df["CHROM"]) == {"1", "7"}


def test_unfittable_rows_dropped(tmp_path):
    """A zero SE would give an infinite Z rather than an obviously-missing result."""
    src = _write(
        tmp_path,
        [
            "1\t1:10:A:G\t10\tG\tA\t0.1\t0.2\t0.02\t1e-20\t1000\t0.15\t0.02\t1e-12",
            "1\t1:11:A:G\t11\tG\tA\t0.1\t0.2\t0.02\t1e-20\t1000\t0.15\t0\t1",
            "1\t1:12:A:G\t12\tG\tA\t0.1\t0.2\t0.02\t1e-20\t1000\tNA\tNA\tNA",
        ],
    )
    out = tmp_path / "reformatted.tsv"
    reformat_mancojo_output([src], out)
    df = pl.read_csv(out, separator=" ", schema_overrides={"CHROM": pl.Utf8})

    assert df["ID"].to_list() == ["1:10:A:G"]


def test_empty_input_writes_header_only(tmp_path):
    src = tmp_path / "empty.cma.cojo"
    src.write_text(HEADER + "\n")
    out = tmp_path / "reformatted.tsv"
    reformat_mancojo_output([src], out)
    df = pl.read_csv(out, separator=" ", schema_overrides={"CHROM": pl.Utf8})

    assert df.height == 0
    assert df.columns == list(MANCOJO_OUTPUT_SCHEMA)


def test_test_label_is_configurable(tmp_path):
    src = _write(tmp_path, ["1\tPCSK9.M1.0.001\t55039974\tA\tR\t0.001\t0.5\t0.1\t1e-06\t1000\t0.4\t0.1\t1e-04"])
    out = tmp_path / "reformatted.tsv"
    reformat_mancojo_output([src], out, test_label="BURDEN")
    df = pl.read_csv(out, separator=" ", schema_overrides={"CHROM": pl.Utf8})

    assert df["TEST"][0] == "BURDEN"
    assert df["ID"][0] == "PCSK9.M1.0.001"
