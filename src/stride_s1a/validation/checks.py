"""Task 5 — dataset validation.

Verifies the derived biological data layer is internally coherent and correctly
mapped back to the S0 inputs:

- every STRIDE locus maps to exactly one canonical residue;
- every canonical residue has exactly one hierarchy path;
- every replicate row maps to a canonical residue (no orphan replicate rows);
- no orphan canonical residues (every residue is backed by a STRIDE locus);
- a shared ``canon_label`` carries consistent structural annotation across
  serotypes.

Each check appends a :class:`ValidationCheck` to the report and raises
:class:`ConsistencyError` on failure (the layer must be trustworthy before later
stages consume it).
"""
from __future__ import annotations

import pandas as pd

from ..models import S1AReport
from ..models.errors import ConsistencyError
from ..models.schema import RESIDUE_SCALE_LEVEL


def validate_locus_mapping(
    stride_table: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    report: S1AReport,
) -> None:
    """Every STRIDE locus maps to exactly one canonical residue, and vice versa.

    The mapping key is (serotype, canon_label) at residue scale. Raises on any
    STRIDE residue-locus with no canonical residue, or any canonical residue not
    backed by a STRIDE locus (orphan residue).
    """
    res = stride_table[stride_table["scale_level"] == RESIDUE_SCALE_LEVEL]
    stride_keys = set(
        zip(res["serotype"], res["canon_label"], strict=True)
    )
    canon_keys = set(
        zip(
            canonical_residues["serotype"],
            canonical_residues["canon_label"],
            strict=True,
        )
    )

    unmapped = stride_keys - canon_keys  # STRIDE locus with no canonical residue
    orphan = canon_keys - stride_keys    # canonical residue with no STRIDE locus
    if unmapped:
        raise ConsistencyError(
            f"{len(unmapped)} STRIDE residue-locus/loci map to no canonical "
            f"residue, e.g. {sorted(unmapped)[:3]}"
        )
    if orphan:
        raise ConsistencyError(
            f"{len(orphan)} canonical residue(s) are orphaned (no STRIDE "
            f"locus), e.g. {sorted(orphan)[:3]}"
        )
    report.add(
        "every STRIDE locus maps to exactly one canonical residue",
        "global",
        True,
        f"{len(stride_keys)} residue loci mapped, no orphans",
    )


def validate_single_hierarchy_path(
    canonical_residues: pd.DataFrame, report: S1AReport
) -> None:
    """Every canonical residue belongs to exactly one hierarchy path."""
    if canonical_residues.empty:
        report.add("every residue has one hierarchy path", "global", True, "empty")
        return
    counts = canonical_residues.groupby(["serotype", "canon_label"])[
        "hierarchy_path"
    ].nunique()
    multi = counts[counts > 1]
    if len(multi):
        raise ConsistencyError(
            f"{len(multi)} residue(s) have more than one hierarchy path, e.g. "
            f"{list(multi.index[:3])}"
        )
    report.add(
        "every canonical residue has exactly one hierarchy path",
        "global",
        True,
        f"{len(counts)} residues",
    )


def validate_replicate_mapping(
    replicate_table: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    report: S1AReport,
) -> None:
    """Every replicate row maps to a canonical residue (no orphan replicate rows).

    A replicate observation of a residue absent from the canonical set would
    mean the two S0 tables disagree; that is a hard error. (The converse — a
    canonical residue with no replicate — is allowed and recorded as
    unavailable by the inventory.)
    """
    if replicate_table.empty:
        report.add(
            "every replicate maps to a canonical residue",
            "global",
            True,
            "no replicate rows (summaries-only dataset)",
        )
        return
    rep_keys = set(
        zip(
            replicate_table["serotype"],
            replicate_table["canon_label"],
            strict=True,
        )
    )
    canon_keys = set(
        zip(
            canonical_residues["serotype"],
            canonical_residues["canon_label"],
            strict=True,
        )
    )
    orphan_rep = rep_keys - canon_keys
    if orphan_rep:
        raise ConsistencyError(
            f"{len(orphan_rep)} replicate residue(s) map to no canonical "
            f"residue, e.g. {sorted(orphan_rep)[:3]}"
        )
    report.add(
        "every replicate maps to a canonical residue",
        "global",
        True,
        f"{len(rep_keys)} (serotype, residue) replicate keys mapped",
    )


def validate_annotation_consistency(
    canonical_residues: pd.DataFrame, report: S1AReport
) -> None:
    """A shared canon_label carries consistent structural annotation.

    For each ``canon_label`` present in multiple serotypes, the (chain, domain,
    motif, secondary_structure) tuple must be identical across serotypes — else
    the cross-serotype conservation mapping would be conflating distinct
    structural positions.
    """
    if canonical_residues.empty:
        report.add("shared canon_label annotation consistent", "global", True, "empty")
        return
    annot_cols = ["chain", "domain", "motif", "secondary_structure"]
    per_label = canonical_residues.groupby("canon_label")[annot_cols].nunique()
    inconsistent = per_label[(per_label > 1).any(axis=1)]
    if len(inconsistent):
        raise ConsistencyError(
            f"{len(inconsistent)} canon_label(s) have inconsistent structural "
            f"annotation across serotypes, e.g. {list(inconsistent.index[:3])}"
        )
    report.add(
        "shared canon_label has consistent annotation across serotypes",
        "global",
        True,
        f"{len(per_label)} distinct canon_labels",
    )
