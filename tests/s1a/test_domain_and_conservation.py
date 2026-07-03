"""Task 2 (conservation) and Task 3 (domain) table tests."""
from __future__ import annotations

import pandas as pd

from stride_s1a.build import (
    build_canonical_residues,
    build_conservation_table,
    build_domain_table,
)
from stride_s1a.models.schema import (
    CONSERVATION_TABLE_COLUMNS,
    DOMAIN_TABLE_COLUMNS,
)


# ---------------------------------------------------------------------------
# Task 3 — domain table
# ---------------------------------------------------------------------------
def test_domain_table_columns_and_key(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    dt = build_domain_table(cr)
    assert list(dt.columns) == list(DOMAIN_TABLE_COLUMNS)
    assert not dt.duplicated(["serotype", "chain", "domain"]).any()


def test_domain_residue_counts(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    dt = build_domain_table(cr)
    # DENVA Catalytic Triad has NS3:51 + NS3:72 = 2 residues
    row = dt[
        (dt.serotype == "DENVA") & (dt.domain == "Catalytic Triad")
    ].iloc[0]
    assert row["n_residues"] == 2
    assert sorted(row["canon_labels"]) == ["NS3:51", "NS3:72"]
    assert row["hierarchy_path"] == "CPLX/protease/NS3/Catalytic Triad"


def test_domain_count_matches_residue_sum(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    dt = build_domain_table(cr)
    # total residues across domains equals number of canonical residues
    assert int(dt["n_residues"].sum()) == len(cr)


def test_domain_table_empty_input() -> None:
    empty = pd.DataFrame(columns=list(
        ["serotype", "canon_label", "chain", "domain", "motif",
         "secondary_structure", "hierarchy_path", "complex", "protein", "residue"]
    ))
    dt = build_domain_table(empty)
    assert dt.empty
    assert list(dt.columns) == list(DOMAIN_TABLE_COLUMNS)


# ---------------------------------------------------------------------------
# Task 2 — conservation table
# ---------------------------------------------------------------------------
def test_conservation_columns_and_key(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    ct = build_conservation_table(cr)
    assert list(ct.columns) == list(CONSERVATION_TABLE_COLUMNS)
    assert not ct.duplicated(["canon_label"]).any()


def test_conservation_union_and_intersection(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    ct = build_conservation_table(cr)
    # union = 3 residues
    assert len(ct) == 3
    # intersection (in all 3 serotypes) = NS3:51 only
    in_all = ct[ct.in_all_serotypes]
    assert list(in_all["canon_label"]) == ["NS3:51"]
    assert int(in_all.iloc[0]["n_serotypes"]) == 3


def test_conservation_presence_absence(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    ct = build_conservation_table(cr)
    # NS3:72 present in DENVA, DENVB; absent from DENVC
    row = ct[ct.canon_label == "NS3:72"].iloc[0]
    assert sorted(row["serotypes_present"]) == ["DENVA", "DENVB"]
    assert sorted(row["serotypes_absent"]) == ["DENVC"]
    assert row["in_all_serotypes"] == False  # noqa: E712
    assert row["in_any_serotype"] == True     # noqa: E712
    # NS3:99 unique to DENVA
    row99 = ct[ct.canon_label == "NS3:99"].iloc[0]
    assert list(row99["serotypes_present"]) == ["DENVA"]
    assert sorted(row99["serotypes_absent"]) == ["DENVB", "DENVC"]
    assert int(row99["n_serotypes"]) == 1


def test_conservation_carries_annotation(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    ct = build_conservation_table(cr)
    row = ct[ct.canon_label == "NS3:51"].iloc[0]
    assert row["chain"] == "NS3"
    assert row["domain"] == "Catalytic Triad"


def test_conservation_empty_input() -> None:
    empty = pd.DataFrame(columns=list(
        ["serotype", "canon_label", "chain", "domain", "motif",
         "secondary_structure", "hierarchy_path", "complex", "protein", "residue"]
    ))
    ct = build_conservation_table(empty)
    assert ct.empty
    assert list(ct.columns) == list(CONSERVATION_TABLE_COLUMNS)
