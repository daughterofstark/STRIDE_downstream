"""Command-line entry point for Stage S5.

Usage::

    python -m stride_s5 --input-dir outputs \\
                        --conservation-input-dir outputs_s1a \\
                        --output-dir outputs_s5

By default the S0 STRIDE table is read from ``--input-dir`` and the S1A
conservation table from ``--conservation-input-dir`` by their standard filenames;
each can be overridden with an explicit path. The provisional ρ* defaults to the
frozen 0.5 and can be overridden with ``--rho-star``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import S5Error
from .models.schema import (
    IN_CONSERVATION_TABLE,
    IN_STRIDE_TABLE,
    OUT_CONSERVATION_SUMMARY,
    PROVISIONAL_RHO_STAR,
)
from .s5 import run_s5


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_s5",
        description="Stage S5: build the cross-serotype layer (n=4).",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("outputs"),
        help="directory holding the S0 STRIDE table (default: outputs)",
    )
    parser.add_argument(
        "--conservation-input-dir",
        type=Path,
        default=Path("outputs_s1a"),
        help=(
            "directory holding the S1A conservation table (default: outputs_s1a)"
        ),
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
        "--conservation-table",
        type=Path,
        default=None,
        help=(
            f"path to {IN_CONSERVATION_TABLE} "
            f"(default: <conservation-input-dir>/{IN_CONSERVATION_TABLE})"
        ),
    )
    parser.add_argument(
        "--rho-star",
        type=float,
        default=PROVISIONAL_RHO_STAR,
        metavar="RHO",
        help=(
            "provisional gate threshold for the reproducibility determination "
            f"(default: {PROVISIONAL_RHO_STAR})"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_s5"),
        help="directory to write S5 artifacts (default: outputs_s5)",
    )
    args = parser.parse_args(argv)

    stride = args.stride_table or (args.input_dir / IN_STRIDE_TABLE)
    conservation = args.conservation_table or (
        args.conservation_input_dir / IN_CONSERVATION_TABLE
    )

    try:
        _tables, report = run_s5(
            stride,
            conservation,
            args.output_dir,
            rho_star=args.rho_star,
        )
    except S5Error as exc:
        print(f"S5 FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S5 OK: {report.n_position_conservation} position-conservation, "
        f"{report.n_direction_concordance} direction-concordance, "
        f"{report.n_domain_serotype_matrix} domain-serotype-matrix, "
        f"{report.n_cross_serotype_scorecard} scorecard rows, "
        f"{len(report.checks)} checks passed."
    )
    print(
        f"Artifacts written to {args.output_dir}/ "
        f"(see {OUT_CONSERVATION_SUMMARY})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
