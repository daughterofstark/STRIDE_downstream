"""Stage-S1A orchestration.

Thin composition of the reusable subpackages: load the S0 canonical tables →
build the four biological tables → validate the derived layer → (optionally)
write artifacts. No statistics, no interpretation.

Public entry points:

- :func:`build_s1a` — load + build + validate; returns the tables + report, no
  file writes.
- :func:`run_s1a`   — additionally writes the artifacts to ``output_dir``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .build import (
    build_canonical_residues,
    build_conservation_table,
    build_domain_table,
    build_replicate_inventory,
)
from .io import (
    load_replicate_table,
    load_stride_table,
    write_dataset_summary,
    write_tables,
)
from .models import S1AReport
from .validation import (
    validate_annotation_consistency,
    validate_locus_mapping,
    validate_replicate_mapping,
    validate_single_hierarchy_path,
)


@dataclass
class S1ATables:
    """The four S1A biological-data-layer tables, returned together."""

    canonical_residues: pd.DataFrame
    domain_table: pd.DataFrame
    replicate_inventory: pd.DataFrame
    conservation_table: pd.DataFrame


def build_s1a(
    stride_table_path: str | Path,
    replicate_table_path: str | Path,
) -> tuple[S1ATables, S1AReport]:
    """Load the S0 tables, build the four S1A tables, and validate them.

    Returns ``(tables, report)``. Raises an :class:`~stride_s1a.S1AError`
    subclass on any input or consistency problem. Does not write to disk.
    """
    report = S1AReport()

    stride = load_stride_table(stride_table_path)
    replicate = load_replicate_table(replicate_table_path)
    report.serotypes = sorted(stride["serotype"].unique().tolist())
    report.add(
        "inputs loaded",
        "global",
        True,
        f"stride_table rows={len(stride)}, replicate_table rows={len(replicate)}",
    )

    # -- build (Tasks 1, 3, 4, 2) -------------------------------------------
    canonical_residues = build_canonical_residues(stride)
    domain_table = build_domain_table(canonical_residues)
    replicate_inventory = build_replicate_inventory(replicate, canonical_residues)
    conservation_table = build_conservation_table(canonical_residues)

    # -- validate (Task 5) --------------------------------------------------
    validate_locus_mapping(stride, canonical_residues, report)
    validate_single_hierarchy_path(canonical_residues, report)
    validate_replicate_mapping(replicate, canonical_residues, report)
    validate_annotation_consistency(canonical_residues, report)

    # -- report facts -------------------------------------------------------
    report.n_canonical_residues = int(len(canonical_residues))
    report.n_domains = int(len(domain_table))
    report.n_union = int(len(conservation_table))
    report.n_conserved_all = (
        int(conservation_table["in_all_serotypes"].sum())
        if not conservation_table.empty
        else 0
    )
    report.facts = {
        "residues_per_serotype": (
            canonical_residues.groupby("serotype").size().to_dict()
            if not canonical_residues.empty
            else {}
        ),
        "domains_per_serotype": (
            domain_table.groupby("serotype").size().to_dict()
            if not domain_table.empty
            else {}
        ),
        "n_residues_with_replicates": (
            int((replicate_inventory["n_replicates"] > 0).sum())
            if not replicate_inventory.empty
            else 0
        ),
    }
    report.add(
        "s1a tables built",
        "global",
        True,
        f"{report.n_canonical_residues} residues, {report.n_domains} domains, "
        f"{report.n_conserved_all}/{report.n_union} conserved in all serotypes",
    )

    tables = S1ATables(
        canonical_residues=canonical_residues,
        domain_table=domain_table,
        replicate_inventory=replicate_inventory,
        conservation_table=conservation_table,
    )
    return tables, report


def run_s1a(
    stride_table_path: str | Path,
    replicate_table_path: str | Path,
    output_dir: str | Path,
) -> tuple[S1ATables, S1AReport]:
    """Full S1A: build + validate, then write the five artifacts to ``output_dir``.

    Writes ``canonical_residues.parquet``, ``domain_table.parquet``,
    ``replicate_inventory.parquet``, ``conservation_table.parquet`` and
    ``dataset_summary.json``. Returns ``(tables, report)``.
    """
    tables, report = build_s1a(stride_table_path, replicate_table_path)
    write_tables(
        tables.canonical_residues,
        tables.domain_table,
        tables.replicate_inventory,
        tables.conservation_table,
        output_dir,
    )
    write_dataset_summary(report, output_dir)
    return tables, report
