"""Stage-S1B orchestration.

Thin composition of the reusable subpackages: load the four S1A tables → build
the four annotation tables → validate the annotation layer → (optionally) write
artifacts. No statistics, no ranking, no figures.

Public entry points:

- :func:`build_s1b` — load + build + validate; returns the tables + report, no
  file writes.
- :func:`run_s1b`   — additionally writes the artifacts to ``output_dir``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .build import (
    build_domain_annotation,
    build_hierarchy_annotation,
    build_residue_annotation,
    build_serotype_annotation,
)
from .io import (
    load_canonical_residues,
    load_conservation_table,
    load_domain_table,
    load_replicate_inventory,
    write_annotation_summary,
    write_tables,
)
from .models import S1BReport
from .models.schema import (
    AVAIL_ALL,
    AVAIL_NONE,
    AVAIL_SOME,
    CONSERVATION_PAN,
    CONSERVATION_PARTIAL,
    CONSERVATION_UNIQUE,
    DOMAIN_STATUS_ASSIGNED,
    DOMAIN_STATUS_UNASSIGNED,
)
from .validation import (
    validate_domain_membership,
    validate_one_annotation_per_residue,
    validate_referential_integrity,
    validate_serotype_references,
    validate_unique_hierarchy_paths,
)


@dataclass
class S1BTables:
    """The four S1B annotation tables, returned together."""

    residue_annotation: pd.DataFrame
    domain_annotation: pd.DataFrame
    hierarchy_annotation: pd.DataFrame
    serotype_annotation: pd.DataFrame


def _column_counts(series: pd.Series, categories: tuple[str, ...]) -> dict[str, int]:
    """Count occurrences of each category in ``series`` (missing → 0)."""
    vc = series.value_counts().to_dict()
    return {cat: int(vc.get(cat, 0)) for cat in categories}


def build_s1b(
    canonical_residues_path: str | Path,
    domain_table_path: str | Path,
    replicate_inventory_path: str | Path,
    conservation_table_path: str | Path,
) -> tuple[S1BTables, S1BReport]:
    """Load the four S1A tables, build the four S1B tables, and validate them.

    Parameters
    ----------
    canonical_residues_path, domain_table_path, replicate_inventory_path, conservation_table_path
        Paths to the four S1A parquet outputs.

    Returns
    -------
    (tables, report)
        The four annotation tables and the run report.

    Raises
    ------
    S1BError
        Any subclass on the first input or consistency problem encountered.
    """
    report = S1BReport()

    canonical_residues = load_canonical_residues(canonical_residues_path)
    domain_table = load_domain_table(domain_table_path)
    replicate_inventory = load_replicate_inventory(replicate_inventory_path)
    conservation_table = load_conservation_table(conservation_table_path)

    report.serotypes = sorted(canonical_residues["serotype"].unique().tolist())
    n_serotypes_total = len(report.serotypes)
    report.add(
        "inputs loaded",
        "global",
        True,
        f"canonical_residues={len(canonical_residues)}, "
        f"domain_table={len(domain_table)}, "
        f"replicate_inventory={len(replicate_inventory)}, "
        f"conservation_table={len(conservation_table)}",
    )

    # -- build --------------------------------------------------------------
    residue_annotation = build_residue_annotation(
        canonical_residues,
        conservation_table,
        replicate_inventory,
        n_serotypes_total,
    )
    domain_annotation = build_domain_annotation(domain_table, residue_annotation)
    hierarchy_annotation = build_hierarchy_annotation(canonical_residues)
    serotype_annotation = build_serotype_annotation(
        residue_annotation, domain_annotation
    )

    # -- validate (structural only) -----------------------------------------
    validate_one_annotation_per_residue(
        residue_annotation, canonical_residues, report
    )
    validate_unique_hierarchy_paths(hierarchy_annotation, report)
    validate_domain_membership(domain_annotation, residue_annotation, report)
    validate_serotype_references(
        serotype_annotation, canonical_residues, residue_annotation, report
    )
    validate_referential_integrity(
        residue_annotation,
        domain_annotation,
        hierarchy_annotation,
        serotype_annotation,
        report,
    )

    # -- report facts -------------------------------------------------------
    report.n_residue_annotations = int(len(residue_annotation))
    report.n_domain_annotations = int(len(domain_annotation))
    report.n_hierarchy_annotations = int(len(hierarchy_annotation))
    report.n_serotype_annotations = int(len(serotype_annotation))
    if not residue_annotation.empty:
        report.facts = {
            "conservation_class_counts": _column_counts(
                residue_annotation["conservation_class"],
                (CONSERVATION_PAN, CONSERVATION_PARTIAL, CONSERVATION_UNIQUE),
            ),
            "availability_class_counts": _column_counts(
                residue_annotation["availability_class"],
                (AVAIL_ALL, AVAIL_SOME, AVAIL_NONE),
            ),
            "domain_status_counts": _column_counts(
                residue_annotation["domain_status"],
                (DOMAIN_STATUS_ASSIGNED, DOMAIN_STATUS_UNASSIGNED),
            ),
            "n_fully_resolved_hierarchy_paths": (
                int(hierarchy_annotation["fully_resolved"].sum())
                if not hierarchy_annotation.empty
                else 0
            ),
        }
    report.add(
        "s1b annotation tables built",
        "global",
        True,
        f"{report.n_residue_annotations} residue, "
        f"{report.n_domain_annotations} domain, "
        f"{report.n_hierarchy_annotations} hierarchy, "
        f"{report.n_serotype_annotations} serotype annotations",
    )

    tables = S1BTables(
        residue_annotation=residue_annotation,
        domain_annotation=domain_annotation,
        hierarchy_annotation=hierarchy_annotation,
        serotype_annotation=serotype_annotation,
    )
    return tables, report


def run_s1b(
    canonical_residues_path: str | Path,
    domain_table_path: str | Path,
    replicate_inventory_path: str | Path,
    conservation_table_path: str | Path,
    output_dir: str | Path,
) -> tuple[S1BTables, S1BReport]:
    """Full S1B: build + validate, then write the five artifacts to ``output_dir``.

    Writes ``residue_annotation.parquet``, ``domain_annotation.parquet``,
    ``hierarchy_annotation.parquet``, ``serotype_annotation.parquet`` and
    ``annotation_summary.json``. Returns ``(tables, report)``.
    """
    tables, report = build_s1b(
        canonical_residues_path,
        domain_table_path,
        replicate_inventory_path,
        conservation_table_path,
    )
    write_tables(
        tables.residue_annotation,
        tables.domain_annotation,
        tables.hierarchy_annotation,
        tables.serotype_annotation,
        output_dir,
    )
    write_annotation_summary(report, output_dir)
    return tables, report
