"""Writers for the S2 output artifacts.

Writes the five reduction parquet tables plus ``reduction_summary.json`` (facts,
provenance header, and validation outcomes).
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..models import S2Report
from ..models.schema import (
    OUT_DOMAIN_REPRODUCIBILITY,
    OUT_REDUCTION_SUMMARY,
    OUT_RESIDUE_LANDSCAPE,
    OUT_RESOLUTION_CENSUS,
    OUT_SEROTYPE_SUMMARY,
    OUT_SIGNED_SCREEN,
)


def write_tables(
    resolution_census: pd.DataFrame,
    residue_landscape: pd.DataFrame,
    domain_reproducibility: pd.DataFrame,
    signed_screen: pd.DataFrame,
    serotype_summary: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write the five S2 reduction tables to ``output_dir``; return the paths."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        OUT_RESOLUTION_CENSUS: out / OUT_RESOLUTION_CENSUS,
        OUT_RESIDUE_LANDSCAPE: out / OUT_RESIDUE_LANDSCAPE,
        OUT_DOMAIN_REPRODUCIBILITY: out / OUT_DOMAIN_REPRODUCIBILITY,
        OUT_SIGNED_SCREEN: out / OUT_SIGNED_SCREEN,
        OUT_SEROTYPE_SUMMARY: out / OUT_SEROTYPE_SUMMARY,
    }
    resolution_census.to_parquet(paths[OUT_RESOLUTION_CENSUS], index=False)
    residue_landscape.to_parquet(paths[OUT_RESIDUE_LANDSCAPE], index=False)
    domain_reproducibility.to_parquet(
        paths[OUT_DOMAIN_REPRODUCIBILITY], index=False
    )
    signed_screen.to_parquet(paths[OUT_SIGNED_SCREEN], index=False)
    serotype_summary.to_parquet(paths[OUT_SEROTYPE_SUMMARY], index=False)
    return paths


def write_reduction_summary(report: S2Report, output_dir: str | Path) -> Path:
    """Write the machine-readable ``reduction_summary.json``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "S2",
        "all_checks_passed": report.all_passed,
        "serotypes": report.serotypes,
        "rho_star_band": report.rho_star_band,
        "provenance": report.provenance,
        "n_resolution_census": report.n_resolution_census,
        "n_residue_landscape": report.n_residue_landscape,
        "n_domain_reproducibility": report.n_domain_reproducibility,
        "n_signed_screen": report.n_signed_screen,
        "n_serotype_summary": report.n_serotype_summary,
        "facts": report.facts,
        "checks": [asdict(c) for c in report.checks],
    }
    path = out / OUT_REDUCTION_SUMMARY
    path.write_text(json.dumps(payload, indent=2))
    return path
