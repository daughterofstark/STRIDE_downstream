"""Command-line entry point for Stage S1A.

Usage::

    python -m stride_s1a --stride-table outputs/stride_table.parquet \\
                         --replicate-table outputs/replicate_table.parquet \\
                         --output-dir outputs_s1a

If ``--stride-table`` / ``--replicate-table`` are omitted, they default to
``stride_table.parquet`` / ``replicate_table.parquet`` inside ``--input-dir``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import S1AError
from .models.schema import OUT_DATASET_SUMMARY
from .s1a import run_s1a


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_s1a",
        description="Stage S1A: build the reusable biological data layer.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("outputs"),
        help="directory holding the S0 canonical tables (default: outputs)",
    )
    parser.add_argument(
        "--stride-table",
        type=Path,
        default=None,
        help="path to stride_table.parquet (default: <input-dir>/stride_table.parquet)",
    )
    parser.add_argument(
        "--replicate-table",
        type=Path,
        default=None,
        help="path to replicate_table.parquet (default: <input-dir>/replicate_table.parquet)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_s1a"),
        help="directory to write S1A artifacts (default: outputs_s1a)",
    )
    args = parser.parse_args(argv)

    stride_path = args.stride_table or (args.input_dir / "stride_table.parquet")
    replicate_path = args.replicate_table or (
        args.input_dir / "replicate_table.parquet"
    )

    try:
        tables, report = run_s1a(stride_path, replicate_path, args.output_dir)
    except S1AError as exc:
        print(f"S1A FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S1A OK: {report.n_canonical_residues} canonical residues, "
        f"{report.n_domains} domains, "
        f"{report.n_conserved_all}/{report.n_union} conserved in all serotypes, "
        f"{len(report.checks)} checks passed."
    )
    print(f"Artifacts written to {args.output_dir}/ (see {OUT_DATASET_SUMMARY})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
