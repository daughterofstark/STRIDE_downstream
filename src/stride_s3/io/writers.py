"""Writers for the S3 output artifacts.

Writes the four hierarchy-reduction parquet tables plus
``hierarchy_summary.json`` (facts, provenance header, and validation outcomes).
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..models import S3Report
from ..models.schema import (
    OUT_CHAIN_CONTRAST,
    OUT_HIERARCHY_SUMMARY,
    OUT_MONOTONICITY_AUDIT,
    OUT_RESOLUTION_GAP,
    OUT_SCALE_CURVE,
)


def write_tables(
    scale_curve: pd.DataFrame,
    resolution_gap: pd.DataFrame,
    monotonicity_audit: pd.DataFrame,
    chain_contrast: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write the four S3 reduction tables to ``output_dir``; return the paths."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        OUT_SCALE_CURVE: out / OUT_SCALE_CURVE,
        OUT_RESOLUTION_GAP: out / OUT_RESOLUTION_GAP,
        OUT_MONOTONICITY_AUDIT: out / OUT_MONOTONICITY_AUDIT,
        OUT_CHAIN_CONTRAST: out / OUT_CHAIN_CONTRAST,
    }
    scale_curve.to_parquet(paths[OUT_SCALE_CURVE], index=False)
    resolution_gap.to_parquet(paths[OUT_RESOLUTION_GAP], index=False)
    monotonicity_audit.to_parquet(paths[OUT_MONOTONICITY_AUDIT], index=False)
    chain_contrast.to_parquet(paths[OUT_CHAIN_CONTRAST], index=False)
    return paths


def write_hierarchy_summary(report: S3Report, output_dir: str | Path) -> Path:
    """Write the machine-readable ``hierarchy_summary.json``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "S3",
        "all_checks_passed": report.all_passed,
        "serotypes": report.serotypes,
        "provenance": report.provenance,
        "n_scale_curve": report.n_scale_curve,
        "n_resolution_gap": report.n_resolution_gap,
        "n_monotonicity_audit": report.n_monotonicity_audit,
        "n_chain_contrast": report.n_chain_contrast,
        "facts": report.facts,
        "checks": [asdict(c) for c in report.checks],
    }
    path = out / OUT_HIERARCHY_SUMMARY
    path.write_text(json.dumps(payload, indent=2))
    return path
