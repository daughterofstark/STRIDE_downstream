"""Command-line entry point for Stage S2.

Usage::

    python -m stride_s2 --stride-input-dir outputs \\
                        --annotation-input-dir outputs_s1b \\
                        --output-dir outputs_s2

By default the S0 STRIDE table is read from ``--stride-input-dir`` and the two
S1B annotation tables from ``--annotation-input-dir`` by their standard
filenames; each can be overridden with an explicit path. The ρ* band defaults to
the frozen 0.50–0.90 sweep and can be overridden with ``--rho-star``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import S2Error
from .models.schema import (
    IN_DOMAIN_ANNOTATION,
    IN_RESIDUE_ANNOTATION,
    IN_STRIDE_TABLE,
    OUT_REDUCTION_SUMMARY,
)
from .s2 import run_s2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_s2",
        description="Stage S2: build the per-serotype reduction layer.",
    )
    parser.add_argument(
        "--stride-input-dir",
        type=Path,
        default=Path("outputs"),
        help="directory holding the S0 STRIDE table (default: outputs)",
    )
    parser.add_argument(
        "--annotation-input-dir",
        type=Path,
        default=Path("outputs_s1b"),
        help="directory holding the S1B annotation tables (default: outputs_s1b)",
    )
    parser.add_argument(
        "--stride-table",
        type=Path,
        default=None,
        help=f"path to {IN_STRIDE_TABLE} (default: <stride-input-dir>/{IN_STRIDE_TABLE})",
    )
    parser.add_argument(
        "--residue-annotation",
        type=Path,
        default=None,
        help=(
            f"path to {IN_RESIDUE_ANNOTATION} "
            f"(default: <annotation-input-dir>/{IN_RESIDUE_ANNOTATION})"
        ),
    )
    parser.add_argument(
        "--domain-annotation",
        type=Path,
        default=None,
        help=(
            f"path to {IN_DOMAIN_ANNOTATION} "
            f"(default: <annotation-input-dir>/{IN_DOMAIN_ANNOTATION})"
        ),
    )
    parser.add_argument(
        "--rho-star",
        type=float,
        nargs="+",
        default=None,
        metavar="RHO",
        help=(
            "ρ* band to sweep (space-separated). Default: the frozen "
            "0.50–0.90 sweep in steps of 0.05."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_s2"),
        help="directory to write S2 artifacts (default: outputs_s2)",
    )
    args = parser.parse_args(argv)

    stride = args.stride_table or (args.stride_input_dir / IN_STRIDE_TABLE)
    residue = args.residue_annotation or (
        args.annotation_input_dir / IN_RESIDUE_ANNOTATION
    )
    domain = args.domain_annotation or (
        args.annotation_input_dir / IN_DOMAIN_ANNOTATION
    )

    try:
        _tables, report = run_s2(
            stride,
            residue,
            domain,
            args.output_dir,
            rho_star_band=args.rho_star,
        )
    except S2Error as exc:
        print(f"S2 FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S2 OK: {report.n_resolution_census} census, "
        f"{report.n_residue_landscape} residue-landscape, "
        f"{report.n_domain_reproducibility} domain, "
        f"{report.n_signed_screen} signed-screen, "
        f"{report.n_serotype_summary} serotype rows, "
        f"{len(report.checks)} checks passed."
    )
    print(f"Artifacts written to {args.output_dir}/ (see {OUT_REDUCTION_SUMMARY})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
