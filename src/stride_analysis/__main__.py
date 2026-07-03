"""Command-line entry point for Stage S0.

Usage::

    python -m stride_analysis --data-root path/to/data --output-dir outputs
    python -m stride_analysis --data-root examples/small_synthetic_dataset --output-dir outputs
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import StrideAnalysisError
from .s0 import run_s0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_analysis",
        description="Stage S0: ingest STRIDE outputs into canonical tables.",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("data"),
        help="dataset root (default: data). Never hardcoded; anything works.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="directory to write S0 artifacts (default: outputs)",
    )
    parser.add_argument(
        "--no-require-replicates",
        action="store_true",
        help="allow serotypes without Level-1 replicate data",
    )
    parser.add_argument(
        "--no-require-summaries",
        action="store_true",
        help="allow serotypes without Level-2 summaries",
    )
    parser.add_argument(
        "--allow-unequal-replicates",
        action="store_true",
        help="do not enforce equal replicate counts across serotypes",
    )
    parser.add_argument(
        "--strict-cross-level",
        action="store_true",
        help="require replicate residue labels to be a subset of profile labels",
    )
    args = parser.parse_args(argv)

    try:
        rep, stride, report = run_s0(
            args.data_root,
            args.output_dir,
            require_replicates=not args.no_require_replicates,
            require_summaries=not args.no_require_summaries,
            enforce_equal_replicate_counts=not args.allow_unequal_replicates,
            strict_cross_level=args.strict_cross_level,
        )
    except StrideAnalysisError as exc:
        print(f"S0 FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S0 OK: {len(report.serotype_facts)} serotype(s), "
        f"{report.replicate_rows} replicate rows, "
        f"{report.stride_rows} stride rows, "
        f"{len(report.checks)} checks passed."
    )
    print(f"Artifacts written to {args.output_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
