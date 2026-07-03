"""Task 1 — canonical residue identifier tests."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s1a.build import build_canonical_residues
from stride_s1a.models.errors import ConsistencyError
from stride_s1a.models.schema import CANONICAL_RESIDUES_COLUMNS


def test_builds_one_row_per_residue(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    # 3 residues in DENVA, 2 in DENVB, 1 in DENVC = 6
    assert len(cr) == 6
    assert list(cr.columns) == list(CANONICAL_RESIDUES_COLUMNS)


def test_key_is_unique(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    assert not cr.duplicated(["serotype", "canon_label"]).any()


def test_hierarchy_fields_populated(stride_table: pd.DataFrame) -> None:
    cr = build_canonical_residues(stride_table)
    row = cr[(cr.serotype == "DENVA") & (cr.canon_label == "NS3:51")].iloc[0]
    assert row["chain"] == "NS3"
    assert row["domain"] == "Catalytic Triad"
    assert row["secondary_structure"] == "helix"
    assert row["complex"] == "CPLX"
    assert row["protein"] == "protease"
    assert row["hierarchy_path"].endswith("NS3:51")


def test_empty_input_yields_empty_table() -> None:
    empty = pd.DataFrame(
        columns=["serotype", "canon_label", "scale_level", "locus", "region_id",
                 "h_complex", "h_protein", "h_chain", "h_domain", "h_motif",
                 "h_secondary_structure", "h_residue"]
    )
    cr = build_canonical_residues(empty)
    assert cr.empty
    assert list(cr.columns) == list(CANONICAL_RESIDUES_COLUMNS)


def test_duplicate_residue_key_raises(stride_table: pd.DataFrame) -> None:
    # duplicate a residue-scale row for the same (serotype, canon_label)
    res_row = stride_table[
        (stride_table.scale_level == "residue")
        & (stride_table.serotype == "DENVA")
        & (stride_table.canon_label == "NS3:51")
    ]
    dup = pd.concat([stride_table, res_row], ignore_index=True)
    with pytest.raises(ConsistencyError, match="not unique"):
        build_canonical_residues(dup)


def test_sentinels_preserved_verbatim(stride_table: pd.DataFrame) -> None:
    # NS3:99 has motif 'none' and secondary_structure 'unknown' — must be kept
    cr = build_canonical_residues(stride_table)
    row = cr[(cr.serotype == "DENVA") & (cr.canon_label == "NS3:99")].iloc[0]
    assert row["motif"] == "none"
    assert row["secondary_structure"] == "unknown"
