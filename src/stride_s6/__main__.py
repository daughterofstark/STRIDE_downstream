"""Command-line entry point for Stage S6.

Usage::

    python -m stride_s6 --input-dir outputs \\
                        --inventory-input-dir outputs_s1a \\
                        --output-dir outputs_s6

By default the S1A replicate inventory is read from ``--inventory-input-dir`` and
the optional S0 replicate table from ``--input-dir`` by their standard filenames;
each can be overridden with an explicit path. The replicate table is optional: if
it is absent (the design's "per-run correlation CSVs unavailable" state, §4.1) the
per-run analyses are recorded as blocked rather than approximated. Pass
``--no-replicate-table`` to force the blocked path even when the file exists.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models.errors import S6Error
from .models.schema import (
    IN_REPLICATE_INVENTORY,
    IN_REPLICATE_TABLE,
    OUT_REPLICATE_SUMMARY,
)
from .s6 import run_s6


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="stride_s6",
        description="Stage S6: build the replicate layer (per-replicate axis).",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("outputs"),
        help=(
            "directory holding the optional S0 replicate table (default: outputs)"
        ),
    )
    parser.add_argument(
        "--inventory-input-dir",
        type=Path,
        default=Path("outputs_s1a"),
        help="directory holding the S1A replicate inventory (default: outputs_s1a)",
    )
    parser.add_argument(
        "--replicate-table",
        type=Path,
        default=None,
        help=(
            f"path to {IN_REPLICATE_TABLE} "
            f"(default: <input-dir>/{IN_REPLICATE_TABLE}; optional)"
        ),
    )
    parser.add_argument(
        "--replicate-inventory",
        type=Path,
        default=None,
        help=(
            f"path to {IN_REPLICATE_INVENTORY} "
            f"(default: <inventory-input-dir>/{IN_REPLICATE_INVENTORY})"
        ),
    )
    parser.add_argument(
        "--no-replicate-table",
        action="store_true",
        help=(
            "ignore the S0 replicate table even if present, forcing the "
            "per-run analyses onto the documented blocked path"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_s6"),
        help="directory to write S6 artifacts (default: outputs_s6)",
    )
    args = parser.parse_args(argv)

    inventory = args.replicate_inventory or (
        args.inventory_input_dir / IN_REPLICATE_INVENTORY
    )
    if args.no_replicate_table:
        replicate_table: Path | None = None
    else:
        replicate_table = args.replicate_table or (
            args.input_dir / IN_REPLICATE_TABLE
        )

    try:
        _tables, report = run_s6(
            inventory,
            args.output_dir,
            replicate_table_path=replicate_table,
        )
    except S6Error as exc:
        print(f"S6 FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print(
        f"S6 OK: {report.n_replicate_regime} regime, "
        f"{report.n_replicate_effect_spread} effect-spread, "
        f"{report.n_replicate_concordance} concordance, "
        f"{report.n_replicate_blocked_analyses} blocked-analysis rows, "
        f"{len(report.checks)} checks passed."
    )
    print(
        f"per-run effects available: {report.per_replicate_effects_available}; "
        f"artifacts written to {args.output_dir}/ (see {OUT_REPLICATE_SUMMARY})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
