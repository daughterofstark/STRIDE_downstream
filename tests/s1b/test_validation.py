"""Tests for S1B structural validation checks."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s1b.build import (
    build_domain_annotation,
    build_hierarchy_annotation,
    build_residue_annotation,
    build_serotype_annotation,
)
from stride_s1b.models import S1BReport
from stride_s1b.models.errors import ConsistencyError
from stride_s1b.validation import (
    validate_domain_membership,
    validate_one_annotation_per_residue,
    validate_referential_integrity,
    validate_serotype_references,
    validate_unique_hierarchy_paths,
)
from tests.s1b.fixtures import N_SEROTYPES_TOTAL

_Built = tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]


@pytest.fixture
def built(
    canonical_residues: pd.DataFrame,
    domain_table: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> _Built:
    """All four annotation tables built from the standard fixtures."""
    ra = build_residue_annotation(
        canonical_residues, conservation_table, replicate_inventory,
        N_SEROTYPES_TOTAL,
    )
    da = build_domain_annotation(domain_table, ra)
    ha = build_hierarchy_annotation(canonical_residues)
    sa = build_serotype_annotation(ra, da)
    return ra, da, ha, sa


# ---------------------------------------------------------------------------
# one annotation per residue / no orphans
# ---------------------------------------------------------------------------
def test_one_annotation_per_residue_passes(built: _Built, canonical_residues: pd.DataFrame) -> None:
    ra, _, _, _ = built
    report = S1BReport()
    validate_one_annotation_per_residue(ra, canonical_residues, report)
    assert report.all_passed


def test_missing_annotation_raises(built: _Built, canonical_residues: pd.DataFrame) -> None:
    ra, _, _, _ = built
    # drop an annotation that still exists in canonical residues
    ra2 = ra[~((ra.serotype == "DENVA") & (ra.canon_label == "NS3:51"))]
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="no annotation"):
        validate_one_annotation_per_residue(ra2, canonical_residues, report)


def test_orphan_annotation_raises(built: _Built, canonical_residues: pd.DataFrame) -> None:
    ra, _, _, _ = built
    extra = ra.iloc[[0]].copy()
    extra["canon_label"] = "NS3:999"
    ra2 = pd.concat([ra, extra], ignore_index=True)
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="orphan"):
        validate_one_annotation_per_residue(ra2, canonical_residues, report)


def test_duplicate_annotation_raises(built: _Built, canonical_residues: pd.DataFrame) -> None:
    ra, _, _, _ = built
    ra2 = pd.concat([ra, ra.iloc[[0]]], ignore_index=True)
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="duplicate"):
        validate_one_annotation_per_residue(ra2, canonical_residues, report)


# ---------------------------------------------------------------------------
# unique hierarchy paths
# ---------------------------------------------------------------------------
def test_unique_hierarchy_paths_passes(built: _Built) -> None:
    _, _, ha, _ = built
    report = S1BReport()
    validate_unique_hierarchy_paths(ha, report)
    assert report.all_passed


def test_duplicate_hierarchy_path_raises(built: _Built) -> None:
    _, _, ha, _ = built
    ha2 = pd.concat([ha, ha.iloc[[0]]], ignore_index=True)
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="not unique"):
        validate_unique_hierarchy_paths(ha2, report)


# ---------------------------------------------------------------------------
# domain membership
# ---------------------------------------------------------------------------
def test_domain_membership_passes(built: _Built) -> None:
    ra, da, _, _ = built
    report = S1BReport()
    validate_domain_membership(da, ra, report)
    assert report.all_passed


def test_domain_membership_unknown_domain_raises(built: _Built) -> None:
    ra, da, _, _ = built
    # corrupt a residue to reference a domain not in the domain annotation
    ra2 = ra.copy()
    ra2.loc[ra2.index[0], "domain"] = "Ghost Domain"
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="absent from the domain annotation"):
        validate_domain_membership(da, ra2, report)


def test_domain_membership_count_mismatch_raises(built: _Built) -> None:
    ra, da, _, _ = built
    da2 = da.copy()
    da2.loc[da2.index[0], "n_residues"] = 999
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="n_residues"):
        validate_domain_membership(da2, ra, report)


# ---------------------------------------------------------------------------
# serotype references
# ---------------------------------------------------------------------------
def test_serotype_references_passes(built: _Built, canonical_residues: pd.DataFrame) -> None:
    ra, _, _, sa = built
    report = S1BReport()
    validate_serotype_references(sa, canonical_residues, ra, report)
    assert report.all_passed


def test_serotype_unknown_reference_raises(built: _Built, canonical_residues: pd.DataFrame) -> None:
    ra, _, _, sa = built
    extra = sa.iloc[[0]].copy()
    extra["serotype"] = "DENVZ"
    sa2 = pd.concat([sa, extra], ignore_index=True)
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="no canonical residues"):
        validate_serotype_references(sa2, canonical_residues, ra, report)


def test_serotype_count_mismatch_raises(built: _Built, canonical_residues: pd.DataFrame) -> None:
    ra, _, _, sa = built
    sa2 = sa.copy()
    sa2.loc[sa2.index[0], "n_residues"] = 999
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="n_residues"):
        validate_serotype_references(sa2, canonical_residues, ra, report)


# ---------------------------------------------------------------------------
# referential integrity
# ---------------------------------------------------------------------------
def test_referential_integrity_passes(built: _Built) -> None:
    ra, da, ha, sa = built
    report = S1BReport()
    validate_referential_integrity(ra, da, ha, sa, report)
    assert report.all_passed


def test_referential_integrity_hierarchy_mismatch_raises(built: _Built) -> None:
    ra, da, ha, sa = built
    # drop a hierarchy row that the residue annotation still references
    ha2 = ha.iloc[1:]
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="hierarchy_path"):
        validate_referential_integrity(ra, da, ha2, sa, report)


def test_referential_integrity_serotype_coverage_raises(built: _Built) -> None:
    ra, da, ha, sa = built
    # remove a serotype from the serotype annotation
    sa2 = sa[sa.serotype != "DENVA"]
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="serotype annotation"):
        validate_referential_integrity(ra, da, ha, sa2, report)


# ---------------------------------------------------------------------------
# empty-input and remaining error branches
# ---------------------------------------------------------------------------
def test_unique_hierarchy_paths_empty_ok() -> None:
    report = S1BReport()
    validate_unique_hierarchy_paths(pd.DataFrame(), report)
    assert report.all_passed


def test_referential_integrity_empty_ok() -> None:
    report = S1BReport()
    validate_referential_integrity(
        pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), report
    )
    assert report.all_passed


def test_duplicate_domain_key_raises(built: _Built) -> None:
    ra, da, _, _ = built
    da2 = pd.concat([da, da.iloc[[0]]], ignore_index=True)
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="domain annotation key"):
        validate_domain_membership(da2, ra, report)


def test_duplicate_serotype_row_raises(built: _Built, canonical_residues: pd.DataFrame) -> None:
    ra, _, _, sa = built
    sa2 = pd.concat([sa, sa.iloc[[0]]], ignore_index=True)
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="duplicate serotype"):
        validate_serotype_references(sa2, canonical_residues, ra, report)


def test_referential_integrity_domain_subset_raises(built: _Built) -> None:
    ra, da, ha, sa = built
    # remove a domain row that residues still reference
    da2 = da[da.domain != "Catalytic Triad"]
    report = S1BReport()
    with pytest.raises(ConsistencyError, match="absent from the domain"):
        validate_referential_integrity(ra, da2, ha, sa, report)
