"""Writers for the S1B output artifacts.

Writes the four annotation parquet tables plus ``annotation_summary.json``.
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..models import S1BReport
from ..models.schema import (
    OUT_ANNOTATION_SUMMARY,
    OUT_DOMAIN_ANNOTATION,
    OUT_HIERARCHY_ANNOTATION,
    OUT_RESIDUE_ANNOTATION,
    OUT_SEROTYPE_ANNOTATION,
)


def write_tables(
    residue_annotation: pd.DataFrame,
    domain_annotation: pd.DataFrame,
    hierarchy_annotation: pd.DataFrame,
    serotype_annotation: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write the four S1B annotation tables to ``output_dir``; return the paths."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        OUT_RESIDUE_ANNOTATION: out / OUT_RESIDUE_ANNOTATION,
        OUT_DOMAIN_ANNOTATION: out / OUT_DOMAIN_ANNOTATION,
        OUT_HIERARCHY_ANNOTATION: out / OUT_HIERARCHY_ANNOTATION,
        OUT_SEROTYPE_ANNOTATION: out / OUT_SEROTYPE_ANNOTATION,
    }
    residue_annotation.to_parquet(paths[OUT_RESIDUE_ANNOTATION], index=False)
    domain_annotation.to_parquet(paths[OUT_DOMAIN_ANNOTATION], index=False)
    hierarchy_annotation.to_parquet(paths[OUT_HIERARCHY_ANNOTATION], index=False)
    serotype_annotation.to_parquet(paths[OUT_SEROTYPE_ANNOTATION], index=False)
    return paths


def write_annotation_summary(report: S1BReport, output_dir: str | Path) -> Path:
    """Write the machine-readable ``annotation_summary.json``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "S1B",
        "all_checks_passed": report.all_passed,
        "serotypes": report.serotypes,
        "n_residue_annotations": report.n_residue_annotations,
        "n_domain_annotations": report.n_domain_annotations,
        "n_hierarchy_annotations": report.n_hierarchy_annotations,
        "n_serotype_annotations": report.n_serotype_annotations,
        "facts": report.facts,
        "checks": [asdict(c) for c in report.checks],
    }
    path = out / OUT_ANNOTATION_SUMMARY
    path.write_text(json.dumps(payload, indent=2))
    return path
