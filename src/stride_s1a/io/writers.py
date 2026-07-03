"""Writers for the S1A output artifacts.

Writes the four parquet tables plus ``dataset_summary.json``. Object-dtype list
columns (e.g. member identifiers) are serialised natively by pyarrow.
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..models import S1AReport
from ..models.schema import (
    OUT_CANONICAL_RESIDUES,
    OUT_CONSERVATION_TABLE,
    OUT_DATASET_SUMMARY,
    OUT_DOMAIN_TABLE,
    OUT_REPLICATE_INVENTORY,
)


def write_tables(
    canonical_residues: pd.DataFrame,
    domain_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
    conservation_table: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write the four S1A tables to ``output_dir``; return the written paths."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        OUT_CANONICAL_RESIDUES: out / OUT_CANONICAL_RESIDUES,
        OUT_DOMAIN_TABLE: out / OUT_DOMAIN_TABLE,
        OUT_REPLICATE_INVENTORY: out / OUT_REPLICATE_INVENTORY,
        OUT_CONSERVATION_TABLE: out / OUT_CONSERVATION_TABLE,
    }
    canonical_residues.to_parquet(paths[OUT_CANONICAL_RESIDUES], index=False)
    domain_table.to_parquet(paths[OUT_DOMAIN_TABLE], index=False)
    replicate_inventory.to_parquet(paths[OUT_REPLICATE_INVENTORY], index=False)
    conservation_table.to_parquet(paths[OUT_CONSERVATION_TABLE], index=False)
    return paths


def write_dataset_summary(report: S1AReport, output_dir: str | Path) -> Path:
    """Write the machine-readable ``dataset_summary.json``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "S1A",
        "all_checks_passed": report.all_passed,
        "serotypes": report.serotypes,
        "n_canonical_residues": report.n_canonical_residues,
        "n_domains": report.n_domains,
        "n_conserved_all_serotypes": report.n_conserved_all,
        "n_union_residues": report.n_union,
        "facts": report.facts,
        "checks": [asdict(c) for c in report.checks],
    }
    path = out / OUT_DATASET_SUMMARY
    path.write_text(json.dumps(payload, indent=2))
    return path
