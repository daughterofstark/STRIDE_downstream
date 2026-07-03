"""Command-line entry point for Stage S1B.

Usage::

    python -m stride_s1b --input-dir outputs_s1a --output-dir outputs_s1b

By default the four S1A tables are read from ``--input-dir`` by their standard
filenames; each can be overridden with an explicit path.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import S1BError
from .models.schema import (
    IN_CANONICAL_RESIDUES,
    IN_CONSERVATION_TABLE,
    IN_DOMAIN_TABLE,
    IN_REPLICATE_INVENTORY,
    OUT_ANNOTATION_SUMMARY,
)
from .s1b import run_s1b


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_s1b",
        description="Stage S1B: build the biological annotation layer.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("outputs_s1a"),
        help="directory holding the S1A tables (default: outputs_s1a)",
    )
    parser.add_argument(
        "--canonical-residues",
        type=Path,
        default=None,
        help=f"path to {IN_CANONICAL_RESIDUES} (default: <input-dir>/{IN_CANONICAL_RESIDUES})",
    )
    parser.add_argument(
        "--domain-table",
        type=Path,
        default=None,
        help=f"path to {IN_DOMAIN_TABLE} (default: <input-dir>/{IN_DOMAIN_TABLE})",
    )
    parser.add_argument(
        "--replicate-inventory",
        type=Path,
        default=None,
        help=f"path to {IN_REPLICATE_INVENTORY} (default: <input-dir>/{IN_REPLICATE_INVENTORY})",
    )
    parser.add_argument(
        "--conservation-table",
        type=Path,
        default=None,
        help=f"path to {IN_CONSERVATION_TABLE} (default: <input-dir>/{IN_CONSERVATION_TABLE})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_s1b"),
        help="directory to write S1B artifacts (default: outputs_s1b)",
    )
    args = parser.parse_args(argv)

    canonical = args.canonical_residues or (args.input_dir / IN_CANONICAL_RESIDUES)
    domain = args.domain_table or (args.input_dir / IN_DOMAIN_TABLE)
    replicate = args.replicate_inventory or (args.input_dir / IN_REPLICATE_INVENTORY)
    conservation = args.conservation_table or (
        args.input_dir / IN_CONSERVATION_TABLE
    )

    try:
        tables, report = run_s1b(
            canonical, domain, replicate, conservation, args.output_dir
        )
    except S1BError as exc:
        print(f"S1B FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S1B OK: {report.n_residue_annotations} residue, "
        f"{report.n_domain_annotations} domain, "
        f"{report.n_hierarchy_annotations} hierarchy, "
        f"{report.n_serotype_annotations} serotype annotations, "
        f"{len(report.checks)} checks passed."
    )
    print(f"Artifacts written to {args.output_dir}/ (see {OUT_ANNOTATION_SUMMARY})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
