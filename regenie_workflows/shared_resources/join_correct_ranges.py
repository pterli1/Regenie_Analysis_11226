#!/usr/bin/env python3
import os

os.environ["POLARS_MAX_THREADS"] = "16"
import polars as pl
import argparse


def _norm_chr() -> pl.Expr:
    return pl.col("CHROM").cast(pl.Utf8).str.replace_all(r"^chr", "").alias("CHROM")


def read_positions(paths: list[str]) -> pl.LazyFrame:
    return (
        pl.scan_csv(
            paths,
            separator="\t",
            has_header=True,
            try_parse_dates=False,
            quote_char=None,
            low_memory=True,
            schema={
                "UNIT": pl.Utf8,
                "MASK": pl.Utf8,
                "CHROM": pl.Utf8,
                "START": pl.Int32,
                "END": pl.Int32,
                "N_VAR": pl.Int32,
            },
        )
        .with_columns(_norm_chr())
        .with_columns((pl.col("UNIT") + "." + pl.col("MASK")).alias("UNIT_MASK"))
    )


def main():
    # Read pos_lf and summary_lf, join on CHROM and ID and replace START/END in summary_lf
    parser = argparse.ArgumentParser(description="Join correct ranges to summary_lf")
    parser.add_argument(
        "-p", "--pos-file", type=str, nargs="+", required=True, help="Path(s) to unit-mask position files"
    )
    parser.add_argument("-s", "--summary-file", type=str, required=True, help="Path to summary file")
    parser.add_argument("-o", "--out-file", type=str, required=True, help="Path to output file")
    args = parser.parse_args()
    pos_lf = read_positions(args.pos_file)
    summary_lf = pl.scan_csv(
        args.summary_file,
        separator="\t",
        has_header=True,
        try_parse_dates=False,
        quote_char=None,
        low_memory=True,
        schema_overrides={
            "CHROM": pl.Utf8,
        },
    )
    summary_lf = summary_lf.with_columns(
        [
            _norm_chr(),
            pl.col("ID").str.replace_all(r"(\.singleton|\.0\.\d+)$", "").alias("UNIT_MASK"),
        ]
    ).drop(["START", "END"])
    merged = summary_lf.join(
        pos_lf.select(["CHROM", "UNIT_MASK", "START", "END"]), on=["CHROM", "UNIT_MASK"], validate="m:1"
    ).drop("UNIT_MASK")
    merged.sink_csv(args.out_file, separator="\t", null_value="NA")


if __name__ == "__main__":
    main()
