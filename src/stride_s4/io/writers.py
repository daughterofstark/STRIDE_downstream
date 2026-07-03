"""Writers for the S4 output artifacts.

Writes the four uncertainty-layer parquet tables plus
``uncertainty_summary.json`` (facts, provenance header, and validation outcomes).
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..models import S4Report
from ..models.schema import (
    OUT_DOMAIN_EFFECT_SUMMARY,
    OUT_RESIDUE_VARIANCE,
    OUT_SIGNIFICANCE_SCREEN,
    OUT_UNCERTAINTY_SUMMARY,
    OUT_VARIANCE_BUDGET,
)


def write_tables(
    variance_budget: pd.DataFrame,
    residue_variance: pd.DataFrame,
    significance_screen: pd.DataFrame,
    domain_effect_summary: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write the four S4 tables to ``output_dir``; return the paths."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        OUT_VARIANCE_BUDGET: out / OUT_VARIANCE_BUDGET,
        OUT_RESIDUE_VARIANCE: out / OUT_RESIDUE_VARIANCE,
        OUT_SIGNIFICANCE_SCREEN: out / OUT_SIGNIFICANCE_SCREEN,
        OUT_DOMAIN_EFFECT_SUMMARY: out / OUT_DOMAIN_EFFECT_SUMMARY,
    }
    variance_budget.to_parquet(paths[OUT_VARIANCE_BUDGET], index=False)
    residue_variance.to_parquet(paths[OUT_RESIDUE_VARIANCE], index=False)
    significance_screen.to_parquet(paths[OUT_SIGNIFICANCE_SCREEN], index=False)
    domain_effect_summary.to_parquet(
        paths[OUT_DOMAIN_EFFECT_SUMMARY], index=False
    )
    return paths


def write_uncertainty_summary(
    report: S4Report, output_dir: str | Path
) -> Path:
    """Write the machine-readable ``uncertainty_summary.json``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "S4",
        "all_checks_passed": report.all_passed,
        "serotypes": report.serotypes,
        "provenance": report.provenance,
        "n_variance_budget": report.n_variance_budget,
        "n_residue_variance": report.n_residue_variance,
        "n_significance_screen": report.n_significance_screen,
        "n_domain_effect_summary": report.n_domain_effect_summary,
        "facts": report.facts,
        "checks": [asdict(c) for c in report.checks],
    }
    path = out / OUT_UNCERTAINTY_SUMMARY
    path.write_text(json.dumps(payload, indent=2))
    return path
