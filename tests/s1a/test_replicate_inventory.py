"""Task 4 — replicate inventory tests (availability only)."""
from __future__ import annotations

import pandas as pd

from stride_s1a.build import build_canonical_residues, build_replicate_inventory
from stride_s1a.models.schema import REPLICATE_INVENTORY_COLUMNS
from tests.s1a.fixtures import make_replicate_table


def test_inventory_columns_and_key(stride_table: pd.DataFrame, replicate_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    ri = build_replicate_inventory(replicate_table, cr)
    assert list(ri.columns) == list(REPLICATE_INVENTORY_COLUMNS)
    assert not ri.duplicated(["serotype", "canon_label"]).any()
    # one row per canonical residue
    assert len(ri) == len(cr)


def test_full_availability(stride_table: pd.DataFrame, replicate_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    ri = build_replicate_inventory(replicate_table, cr)
    row = ri[(ri.serotype == "DENVA") & (ri.canon_label == "NS3:51")].iloc[0]
    assert row["n_replicates"] == 3
    assert sorted(row["replicates"]) == ["1st_run", "2nd_run", "3rd_run"]
    assert sorted(row["replicate_indices"]) == [1, 2, 3]
    assert row["available"] == True          # noqa: E712
    assert row["in_all_replicates"] == True  # noqa: E712


def test_partial_availability(stride_table: pd.DataFrame, replicate_table: pd.DataFrame) -> None:
    # NS3:72 in DENVB observed in only the first replicate
    cr = build_canonical_residues(stride_table)
    ri = build_replicate_inventory(replicate_table, cr)
    row = ri[(ri.serotype == "DENVB") & (ri.canon_label == "NS3:72")].iloc[0]
    assert row["n_replicates"] == 1
    assert list(row["replicates"]) == ["1st_run"]
    assert row["available"] == True           # noqa: E712
    assert row["in_all_replicates"] == False  # noqa: E712 (DENVB has 3 replicates)


def test_no_averaging_columns_present() -> None:
    # the inventory must not carry any averaged/statistical column
    forbidden = {"mean", "avg", "median", "std", "sem", "pvalue", "p_value",
                 "effect", "r", "abs_r", "rho", "beta"}
    assert forbidden.isdisjoint(set(REPLICATE_INVENTORY_COLUMNS))


def test_summaries_only_all_unavailable(stride_table: pd.DataFrame) -> None:
    # empty replicate table -> every residue recorded as unavailable, not dropped
    cr = build_canonical_residues(stride_table)
    empty_rep = make_replicate_table().iloc[0:0]
    ri = build_replicate_inventory(empty_rep, cr)
    assert len(ri) == len(cr)
    assert (ri["n_replicates"] == 0).all()
    assert (~ri["available"]).all()
    assert (~ri["in_all_replicates"]).all()


def test_residue_absent_from_replicates_recorded(stride_table: pd.DataFrame) -> None:
    # NS3:99 exists in canonical (DENVA) but drop it from replicates entirely
    cr = build_canonical_residues(stride_table)
    rep = make_replicate_table()
    rep = rep[~((rep.serotype == "DENVA") & (rep.canon_label == "NS3:99"))]
    ri = build_replicate_inventory(rep, cr)
    row = ri[(ri.serotype == "DENVA") & (ri.canon_label == "NS3:99")].iloc[0]
    assert row["n_replicates"] == 0
    assert row["available"] == False  # noqa: E712
    # still present as a row (not dropped)
    assert len(ri) == len(cr)
