"""Command-line entry point for Stage S4.

Usage::

    python -m stride_s4 --input-dir outputs --output-dir outputs_s4

By default the S0 STRIDE table is read from ``--input-dir`` by its standard
filename; it can be overridden with an explicit ``--stride-table`` path.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import S4Error
from .models.schema import IN_STRIDE_TABLE, OUT_UNCERTAINTY_SUMMARY
from .s4 import run_s4


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_s4",
        description="Stage S4: build the uncertainty layer.",
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
        default=Path("outputs_s4"),
        help="directory to write S4 artifacts (default: outputs_s4)",
    )
    args = parser.parse_args(argv)

    stride = args.stride_table or (args.input_dir / IN_STRIDE_TABLE)

    try:
        _tables, report = run_s4(stride, args.output_dir)
    except S4Error as exc:
        print(f"S4 FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S4 OK: {report.n_variance_budget} variance-budget, "
        f"{report.n_residue_variance} residue-variance, "
        f"{report.n_significance_screen} significance-screen, "
        f"{report.n_domain_effect_summary} domain-effect rows, "
        f"{len(report.checks)} checks passed."
    )
    print(
        f"Artifacts written to {args.output_dir}/ "
        f"(see {OUT_UNCERTAINTY_SUMMARY})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
