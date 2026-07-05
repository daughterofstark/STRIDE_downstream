"""Command-line entry point for Stage S7 (figures & tables).

Usage::

    python -m stride_s7 --s2-input-dir outputs_s2 \\
                        --s3-input-dir outputs_s3 \\
                        --s4-input-dir outputs_s4 \\
                        --s5-input-dir outputs_s5 \\
                        --s6-input-dir outputs_s6 \\
                        --output-dir outputs_s7

Each ``--sN-input-dir`` points at the output directory of the corresponding prior
stage; S7 reads the reduction tables it needs from those directories (never the raw
STRIDE files). All required inputs must be present — a missing one is a hard
failure.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import S7Error
from .models.schema import OUT_SUMMARY
from .s7 import DEFAULT_STAGE_DIRS, run_s7


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_s7",
        description=(
            "Stage S7: assemble publication figures (F1-F8) and manuscript "
            "tables (T1-T5) from the S2-S6 outputs. Reporting only."
        ),
    )
    parser.add_argument(
        "--s2-input-dir",
        type=Path,
        default=Path(DEFAULT_STAGE_DIRS["s2"]),
        help=f"S2 output directory (default: {DEFAULT_STAGE_DIRS['s2']})",
    )
    parser.add_argument(
        "--s3-input-dir",
        type=Path,
        default=Path(DEFAULT_STAGE_DIRS["s3"]),
        help=f"S3 output directory (default: {DEFAULT_STAGE_DIRS['s3']})",
    )
    parser.add_argument(
        "--s4-input-dir",
        type=Path,
        default=Path(DEFAULT_STAGE_DIRS["s4"]),
        help=f"S4 output directory (default: {DEFAULT_STAGE_DIRS['s4']})",
    )
    parser.add_argument(
        "--s5-input-dir",
        type=Path,
        default=Path(DEFAULT_STAGE_DIRS["s5"]),
        help=f"S5 output directory (default: {DEFAULT_STAGE_DIRS['s5']})",
    )
    parser.add_argument(
        "--s6-input-dir",
        type=Path,
        default=Path(DEFAULT_STAGE_DIRS["s6"]),
        help=f"S6 output directory (default: {DEFAULT_STAGE_DIRS['s6']})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_s7"),
        help="directory to write S7 artifacts (default: outputs_s7)",
    )
    args = parser.parse_args(argv)

    stage_dirs = {
        "s2": args.s2_input_dir,
        "s3": args.s3_input_dir,
        "s4": args.s4_input_dir,
        "s5": args.s5_input_dir,
        "s6": args.s6_input_dir,
    }

    try:
        _artifacts, report = run_s7(stage_dirs, args.output_dir)
    except S7Error as exc:
        print(f"S7 FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S7 OK: {len(report.figures)} figures, {len(report.tables)} tables, "
        f"{len(report.limitations)} limitation rows, "
        f"{len(report.checks)} checks passed."
    )
    print(
        f"artifacts written to {args.output_dir}/ (see {OUT_SUMMARY})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
