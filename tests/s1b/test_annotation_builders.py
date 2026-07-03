"""Tests for the domain, hierarchy, and serotype annotation builders."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s1b.build import (
    build_domain_annotation,
    build_hierarchy_annotation,
    build_residue_annotation,
    build_serotype_annotation,
)
from stride_s1b.models.errors import ConsistencyError
from stride_s1b.models.schema import (
    DOMAIN_ANNOTATION_COLUMNS,
    DOMAIN_STATUS_ASSIGNED,
    DOMAIN_STATUS_UNASSIGNED,
    HIERARCHY_ANNOTATION_COLUMNS,
    N_HIERARCHY_LEVELS,
    SEROTYPE_ANNOTATION_COLUMNS,
)
from tests.s1b.fixtures import N_SEROTYPES_TOTAL


def _residue_annotation(
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> pd.DataFrame:
    return build_residue_annotation(
        canonical_residues,
        conservation_table,
        replicate_inventory,
        N_SEROTYPES_TOTAL,
    )


# ---------------------------------------------------------------------------
# domain annotation
# ---------------------------------------------------------------------------
def test_domain_annotation_columns_and_key(
    domain_table: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _residue_annotation(canonical_residues, conservation_table, replicate_inventory)
    da = build_domain_annotation(domain_table, ra)
    assert list(da.columns) == list(DOMAIN_ANNOTATION_COLUMNS)
    assert not da.duplicated(["serotype", "chain", "domain"]).any()


def test_domain_status_and_counts(
    domain_table: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _residue_annotation(canonical_residues, conservation_table, replicate_inventory)
    da = build_domain_annotation(domain_table, ra)
    # DENVA Catalytic Triad: NS3:51 (pan) + NS3:72 (partial) = 2 residues, 1 pan
    row = da[
        (da.serotype == "DENVA") & (da.domain == "Catalytic Triad")
    ].iloc[0]
    assert row["domain_status"] == DOMAIN_STATUS_ASSIGNED
    assert row["n_residues"] == 2
    assert row["n_pan_serotype_residues"] == 1
    assert row["fully_conserved"] == False  # noqa: E712
    # DENVA unassigned domain: NS3:99 only, unassigned status
    row99 = da[(da.serotype == "DENVA") & (da.domain == "unassigned")].iloc[0]
    assert row99["domain_status"] == DOMAIN_STATUS_UNASSIGNED


def test_domain_fully_conserved(
    domain_table: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _residue_annotation(canonical_residues, conservation_table, replicate_inventory)
    da = build_domain_annotation(domain_table, ra)
    # DENVC Catalytic Triad contains only NS3:51 (pan) -> fully conserved
    row = da[(da.serotype == "DENVC") & (da.domain == "Catalytic Triad")].iloc[0]
    assert row["n_residues"] == 1
    assert row["n_pan_serotype_residues"] == 1
    assert row["fully_conserved"] == True  # noqa: E712


def test_domain_annotation_empty() -> None:
    da = build_domain_annotation(pd.DataFrame(), pd.DataFrame())
    assert da.empty
    assert list(da.columns) == list(DOMAIN_ANNOTATION_COLUMNS)


def test_domain_annotation_empty_residue_annotation(domain_table: pd.DataFrame) -> None:
    # domain table present but no residue annotation -> pan counts all zero
    da = build_domain_annotation(domain_table, pd.DataFrame())
    assert len(da) == len(domain_table.drop_duplicates(["serotype", "chain", "domain"]))
    assert (da["n_pan_serotype_residues"] == 0).all()


# ---------------------------------------------------------------------------
# hierarchy annotation
# ---------------------------------------------------------------------------
def test_hierarchy_annotation_columns_and_key(canonical_residues: pd.DataFrame) -> None:
    ha = build_hierarchy_annotation(canonical_residues)
    assert list(ha.columns) == list(HIERARCHY_ANNOTATION_COLUMNS)
    assert not ha.duplicated(["serotype", "hierarchy_path"]).any()
    assert len(ha) == len(canonical_residues)


def test_hierarchy_resolution_counts(canonical_residues: pd.DataFrame) -> None:
    ha = build_hierarchy_annotation(canonical_residues)
    # NS3:51: complex/protein/chain/domain resolved, motif "none" unresolved,
    # ss "helix" resolved, residue resolved -> 6 resolved of 7
    r51 = ha[(ha.serotype == "DENVA") & (ha.canon_label == "NS3:51")].iloc[0]
    assert r51["n_levels_total"] == N_HIERARCHY_LEVELS
    assert r51["n_levels_resolved"] == 6
    assert r51["fully_resolved"] == False  # noqa: E712
    # NS3:99: domain "unassigned" + motif "none" unresolved, ss "sheet" resolved
    # -> complex,protein,chain,ss,residue resolved = 5
    r99 = ha[(ha.serotype == "DENVA") & (ha.canon_label == "NS3:99")].iloc[0]
    assert r99["n_levels_resolved"] == 5


def test_hierarchy_duplicate_path_raises(canonical_residues: pd.DataFrame) -> None:
    dup = pd.concat(
        [canonical_residues, canonical_residues.iloc[[0]]], ignore_index=True
    )
    with pytest.raises(ConsistencyError, match="not unique"):
        build_hierarchy_annotation(dup)


def test_hierarchy_annotation_empty() -> None:
    empty = pd.DataFrame(
        columns=[
            "serotype", "canon_label", "chain", "domain", "motif",
            "secondary_structure", "hierarchy_path", "complex", "protein",
            "residue",
        ]
    )
    ha = build_hierarchy_annotation(empty)
    assert ha.empty
    assert list(ha.columns) == list(HIERARCHY_ANNOTATION_COLUMNS)


# ---------------------------------------------------------------------------
# serotype annotation
# ---------------------------------------------------------------------------
def test_serotype_annotation_columns_and_key(
    domain_table: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _residue_annotation(canonical_residues, conservation_table, replicate_inventory)
    da = build_domain_annotation(domain_table, ra)
    sa = build_serotype_annotation(ra, da)
    assert list(sa.columns) == list(SEROTYPE_ANNOTATION_COLUMNS)
    assert not sa.duplicated(["serotype"]).any()


def test_serotype_counts(
    domain_table: pd.DataFrame,
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _residue_annotation(canonical_residues, conservation_table, replicate_inventory)
    da = build_domain_annotation(domain_table, ra)
    sa = build_serotype_annotation(ra, da)
    # DENVA has 3 residues: NS3:51(pan), NS3:72(partial), NS3:99(unique)
    a = sa[sa.serotype == "DENVA"].iloc[0]
    assert a["n_residues"] == 3
    assert a["n_pan_serotype_residues"] == 1
    assert a["n_partial_residues"] == 1
    assert a["n_serotype_unique_residues"] == 1
    # domain-status tally: NS3:51/NS3:72 assigned, NS3:99 unassigned
    assert a["n_assigned_domain_residues"] == 2
    assert a["n_unassigned_domain_residues"] == 1
    # counts sum to n_residues
    assert (
        a["n_pan_serotype_residues"]
        + a["n_partial_residues"]
        + a["n_serotype_unique_residues"]
        == a["n_residues"]
    )


def test_serotype_annotation_empty() -> None:
    sa = build_serotype_annotation(pd.DataFrame(), pd.DataFrame())
    assert sa.empty
    assert list(sa.columns) == list(SEROTYPE_ANNOTATION_COLUMNS)


def test_serotype_annotation_empty_domain_annotation(
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    # residues present but no domain annotation -> n_domains all zero
    ra = _residue_annotation(canonical_residues, conservation_table, replicate_inventory)
    sa = build_serotype_annotation(ra, pd.DataFrame())
    assert (sa["n_domains"] == 0).all()
    assert len(sa) == ra["serotype"].nunique()
