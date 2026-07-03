"""Command-line entry point for Stage S3.

Usage::

    python -m stride_s3 --input-dir outputs --output-dir outputs_s3

By default the S0 STRIDE table is read from ``--input-dir`` by its standard
filename; it can be overridden with an explicit ``--stride-table`` path.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import S3Error
from .models.schema import IN_STRIDE_TABLE, OUT_HIERARCHY_SUMMARY
from .s3 import run_s3


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_s3",
        description="Stage S3: build the hierarchy reduction layer.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("outputs"),
        help="directory holding the S0 STRIDE table (default: outputs)",
    )
    parser.add_argument(
        "--stride-table",
        type=Path,
        default=None,
        help=(
            f"path to {IN_STRIDE_TABLE} "
            f"(default: <input-dir>/{IN_STRIDE_TABLE})"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_s3"),
        help="directory to write S3 artifacts (default: outputs_s3)",
    )
    args = parser.parse_args(argv)

    stride = args.stride_table or (args.input_dir / IN_STRIDE_TABLE)

    try:
        _tables, report = run_s3(stride, args.output_dir)
    except S3Error as exc:
        print(f"S3 FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S3 OK: {report.n_scale_curve} scale-curve, "
        f"{report.n_resolution_gap} gap, "
        f"{report.n_monotonicity_audit} monotonicity, "
        f"{report.n_chain_contrast} chain-contrast rows, "
        f"{len(report.checks)} checks passed."
    )
    print(f"Artifacts written to {args.output_dir}/ (see {OUT_HIERARCHY_SUMMARY})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
