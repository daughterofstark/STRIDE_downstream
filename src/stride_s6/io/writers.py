"""Writers for the S6 output artifacts.

Writes the four replicate-layer parquet tables plus ``replicate_summary.json``
(facts, provenance header, the blocked-analysis ledger, and validation outcomes).
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..models import S6Report
from ..models.schema import (
    OUT_REPLICATE_BLOCKED_ANALYSES,
    OUT_REPLICATE_CONCORDANCE,
    OUT_REPLICATE_EFFECT_SPREAD,
    OUT_REPLICATE_REGIME,
    OUT_REPLICATE_SUMMARY,
)


def write_tables(
    replicate_regime: pd.DataFrame,
    replicate_effect_spread: pd.DataFrame,
    replicate_concordance: pd.DataFrame,
    replicate_blocked_analyses: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write the four S6 tables to ``output_dir``; return the paths."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        OUT_REPLICATE_REGIME: out / OUT_REPLICATE_REGIME,
        OUT_REPLICATE_EFFECT_SPREAD: out / OUT_REPLICATE_EFFECT_SPREAD,
        OUT_REPLICATE_CONCORDANCE: out / OUT_REPLICATE_CONCORDANCE,
        OUT_REPLICATE_BLOCKED_ANALYSES: out / OUT_REPLICATE_BLOCKED_ANALYSES,
    }
    replicate_regime.to_parquet(paths[OUT_REPLICATE_REGIME], index=False)
    replicate_effect_spread.to_parquet(
        paths[OUT_REPLICATE_EFFECT_SPREAD], index=False
    )
    replicate_concordance.to_parquet(
        paths[OUT_REPLICATE_CONCORDANCE], index=False
    )
    replicate_blocked_analyses.to_parquet(
        paths[OUT_REPLICATE_BLOCKED_ANALYSES], index=False
    )
    return paths


def write_replicate_summary(report: S6Report, output_dir: str | Path) -> Path:
    """Write the machine-readable ``replicate_summary.json``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "S6",
        "all_checks_passed": report.all_passed,
        "serotypes": report.serotypes,
        "per_replicate_effects_available": report.per_replicate_effects_available,
        "provenance": report.provenance,
        "n_replicate_regime": report.n_replicate_regime,
        "n_replicate_effect_spread": report.n_replicate_effect_spread,
        "n_replicate_concordance": report.n_replicate_concordance,
        "n_replicate_blocked_analyses": report.n_replicate_blocked_analyses,
        "facts": report.facts,
        "blocked_analyses": report.blocked_analyses,
        "checks": [asdict(c) for c in report.checks],
    }
    path = out / OUT_REPLICATE_SUMMARY
    path.write_text(json.dumps(payload, indent=2))
    return path
