"""Task 5 — dataset validation tests."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s1a.build import build_canonical_residues
from stride_s1a.models import S1AReport
from stride_s1a.models.errors import ConsistencyError
from stride_s1a.validation import (
    validate_annotation_consistency,
    validate_locus_mapping,
    validate_replicate_mapping,
    validate_single_hierarchy_path,
)


def test_locus_mapping_passes(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    report = S1AReport()
    validate_locus_mapping(stride_table, cr, report)
    assert report.all_passed


def test_orphan_residue_raises(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    # inject a canonical residue with no STRIDE locus
    extra = cr.iloc[[0]].copy()
    extra["canon_label"] = "NS3:999"
    cr2 = pd.concat([cr, extra], ignore_index=True)
    report = S1AReport()
    with pytest.raises(ConsistencyError, match="orphan"):
        validate_locus_mapping(stride_table, cr2, report)


def test_unmapped_locus_raises(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    # drop a residue from the canonical set that still exists in stride_table
    cr2 = cr[~((cr.serotype == "DENVA") & (cr.canon_label == "NS3:51"))]
    report = S1AReport()
    with pytest.raises(ConsistencyError, match="map to no canonical residue"):
        validate_locus_mapping(stride_table, cr2, report)


def test_single_hierarchy_path_passes(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    report = S1AReport()
    validate_single_hierarchy_path(cr, report)
    assert report.all_passed


def test_multiple_hierarchy_paths_raises(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    # corrupt one residue to have two different paths
    dupe = cr[(cr.serotype == "DENVA") & (cr.canon_label == "NS3:51")].copy()
    dupe["hierarchy_path"] = "DIFFERENT/PATH"
    cr2 = pd.concat([cr, dupe], ignore_index=True)
    report = S1AReport()
    with pytest.raises(ConsistencyError, match="more than one hierarchy path"):
        validate_single_hierarchy_path(cr2, report)


def test_replicate_mapping_passes(stride_table: pd.DataFrame, replicate_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    report = S1AReport()
    validate_replicate_mapping(replicate_table, cr, report)
    assert report.all_passed


def test_orphan_replicate_raises(stride_table: pd.DataFrame, replicate_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    # a replicate row for a residue not in the canonical set
    bad = replicate_table.iloc[[0]].copy()
    bad["canon_label"] = "NS3:404"
    rep2 = pd.concat([replicate_table, bad], ignore_index=True)
    report = S1AReport()
    with pytest.raises(ConsistencyError, match="replicate residue"):
        validate_replicate_mapping(rep2, cr, report)


def test_replicate_mapping_empty_ok(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    report = S1AReport()
    validate_replicate_mapping(replicate_table_empty(), cr, report)
    assert report.all_passed


def replicate_table_empty() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["serotype", "replicate", "replicate_index", "canon_label"]
    )


def test_annotation_consistency_passes(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    report = S1AReport()
    validate_annotation_consistency(cr, report)
    assert report.all_passed


def test_inconsistent_annotation_raises(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    # make NS3:51 have a different domain in DENVB than DENVA
    mask = (cr.serotype == "DENVB") & (cr.canon_label == "NS3:51")
    cr.loc[mask, "domain"] = "Totally Different Domain"
    report = S1AReport()
    with pytest.raises(ConsistencyError, match="inconsistent structural annotation"):
        validate_annotation_consistency(cr, report)


def test_empty_canonical_residues_validation_passes() -> None:
    # all checks must handle an empty canonical set without error
    empty = pd.DataFrame(columns=list(
        ["serotype", "canon_label", "chain", "domain", "motif",
         "secondary_structure", "hierarchy_path", "complex", "protein", "residue"]
    ))
    report = S1AReport()
    validate_single_hierarchy_path(empty, report)
    validate_annotation_consistency(empty, report)
    assert report.all_passed
