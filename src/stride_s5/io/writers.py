"""Writers for the S5 output artifacts.

Writes the four cross-serotype parquet tables plus ``conservation_summary.json``
(facts, provenance header, and validation outcomes).
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..models import S5Report
from ..models.schema import (
    OUT_CONSERVATION_SUMMARY,
    OUT_CROSS_SEROTYPE_SCORECARD,
    OUT_DIRECTION_CONCORDANCE,
    OUT_DOMAIN_SEROTYPE_MATRIX,
    OUT_POSITION_CONSERVATION,
)


def write_tables(
    position_conservation: pd.DataFrame,
    direction_concordance: pd.DataFrame,
    domain_serotype_matrix: pd.DataFrame,
    cross_serotype_scorecard: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write the four S5 tables to ``output_dir``; return the paths."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        OUT_POSITION_CONSERVATION: out / OUT_POSITION_CONSERVATION,
        OUT_DIRECTION_CONCORDANCE: out / OUT_DIRECTION_CONCORDANCE,
        OUT_DOMAIN_SEROTYPE_MATRIX: out / OUT_DOMAIN_SEROTYPE_MATRIX,
        OUT_CROSS_SEROTYPE_SCORECARD: out / OUT_CROSS_SEROTYPE_SCORECARD,
    }
    position_conservation.to_parquet(
        paths[OUT_POSITION_CONSERVATION], index=False
    )
    direction_concordance.to_parquet(
        paths[OUT_DIRECTION_CONCORDANCE], index=False
    )
    domain_serotype_matrix.to_parquet(
        paths[OUT_DOMAIN_SEROTYPE_MATRIX], index=False
    )
    cross_serotype_scorecard.to_parquet(
        paths[OUT_CROSS_SEROTYPE_SCORECARD], index=False
    )
    return paths


def write_conservation_summary(
    report: S5Report, output_dir: str | Path
) -> Path:
    """Write the machine-readable ``conservation_summary.json``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "S5",
        "all_checks_passed": report.all_passed,
        "serotypes": report.serotypes,
        "provenance": report.provenance,
        "n_position_conservation": report.n_position_conservation,
        "n_direction_concordance": report.n_direction_concordance,
        "n_domain_serotype_matrix": report.n_domain_serotype_matrix,
        "n_cross_serotype_scorecard": report.n_cross_serotype_scorecard,
        "facts": report.facts,
        "checks": [asdict(c) for c in report.checks],
    }
    path = out / OUT_CONSERVATION_SUMMARY
    path.write_text(json.dumps(payload, indent=2))
    return path
