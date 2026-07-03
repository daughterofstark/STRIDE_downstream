"""Tests for the residue annotation builder."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s1b.build import build_residue_annotation
from stride_s1b.models.errors import ConsistencyError
from stride_s1b.models.schema import (
    AVAIL_ALL,
    AVAIL_NONE,
    AVAIL_SOME,
    CONSERVATION_PAN,
    CONSERVATION_PARTIAL,
    CONSERVATION_UNIQUE,
    DOMAIN_STATUS_ASSIGNED,
    DOMAIN_STATUS_UNASSIGNED,
    RESIDUE_ANNOTATION_COLUMNS,
    SS_STATUS_RESOLVED,
    SS_STATUS_UNRESOLVED,
)
from tests.s1b.fixtures import (
    N_SEROTYPES_TOTAL,
    make_empty_replicate_inventory,
)


def _build(
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


def test_columns_and_key(
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _build(canonical_residues, conservation_table, replicate_inventory)
    assert list(ra.columns) == list(RESIDUE_ANNOTATION_COLUMNS)
    assert not ra.duplicated(["serotype", "canon_label"]).any()
    assert len(ra) == len(canonical_residues)


def test_conservation_classes(
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _build(canonical_residues, conservation_table, replicate_inventory)

    def cls(serotype: str, label: str) -> str:
        row = ra[(ra.serotype == serotype) & (ra.canon_label == label)].iloc[0]
        return row["conservation_class"]

    assert cls("DENVA", "NS3:51") == CONSERVATION_PAN      # in all 3
    assert cls("DENVA", "NS3:72") == CONSERVATION_PARTIAL  # in 2 of 3
    assert cls("DENVA", "NS3:99") == CONSERVATION_UNIQUE   # in 1 of 3


def test_domain_and_ss_status(
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _build(canonical_residues, conservation_table, replicate_inventory)
    r51 = ra[(ra.serotype == "DENVA") & (ra.canon_label == "NS3:51")].iloc[0]
    assert r51["domain_status"] == DOMAIN_STATUS_ASSIGNED
    assert r51["secondary_structure_status"] == SS_STATUS_RESOLVED  # helix
    r99 = ra[(ra.serotype == "DENVA") & (ra.canon_label == "NS3:99")].iloc[0]
    assert r99["domain_status"] == DOMAIN_STATUS_UNASSIGNED  # "unassigned"
    r72 = ra[(ra.serotype == "DENVA") & (ra.canon_label == "NS3:72")].iloc[0]
    assert r72["secondary_structure_status"] == SS_STATUS_UNRESOLVED  # unknown


def test_availability_classes(
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    ra = _build(canonical_residues, conservation_table, replicate_inventory)
    # NS3:51 DENVA observed in all replicates
    r51 = ra[(ra.serotype == "DENVA") & (ra.canon_label == "NS3:51")].iloc[0]
    assert r51["availability_class"] == AVAIL_ALL
    assert r51["n_replicates"] == 3
    # NS3:72 in DENVB observed in only one replicate
    r72b = ra[(ra.serotype == "DENVB") & (ra.canon_label == "NS3:72")].iloc[0]
    assert r72b["availability_class"] == AVAIL_SOME
    assert r72b["n_replicates"] == 1


def test_no_replicates_when_inventory_empty(
    canonical_residues: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    ri = make_empty_replicate_inventory()
    ra = _build(canonical_residues, conservation_table, ri)
    assert (ra["availability_class"] == AVAIL_NONE).all()
    assert (ra["n_replicates"] == 0).all()


def test_empty_canonical_residues() -> None:
    ra = build_residue_annotation(
        pd.DataFrame(
            columns=[
                "serotype", "canon_label", "chain", "domain", "motif",
                "secondary_structure", "hierarchy_path", "complex", "protein",
                "residue",
            ]
        ),
        pd.DataFrame(),
        pd.DataFrame(),
        0,
    )
    assert ra.empty
    assert list(ra.columns) == list(RESIDUE_ANNOTATION_COLUMNS)


def test_duplicate_residue_key_raises(
    canonical_residues: pd.DataFrame,
    conservation_table: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
) -> None:
    dup = pd.concat(
        [canonical_residues, canonical_residues.iloc[[0]]], ignore_index=True
    )
    with pytest.raises(ConsistencyError, match="not unique"):
        _build(dup, conservation_table, replicate_inventory)


def test_no_statistical_columns() -> None:
    forbidden = {"mean", "median", "std", "sem", "pvalue", "p_value", "score",
                 "rank", "effect", "rho", "beta"}
    assert forbidden.isdisjoint(set(RESIDUE_ANNOTATION_COLUMNS))


def test_empty_conservation_table_defaults(
    canonical_residues: pd.DataFrame, replicate_inventory: pd.DataFrame
) -> None:
    # conservation table empty -> every residue defaults to 1 serotype present
    ra = build_residue_annotation(
        canonical_residues, pd.DataFrame(), replicate_inventory,
        N_SEROTYPES_TOTAL,
    )
    # with n_present defaulting to 1 of 3, all become serotype_unique
    assert (ra["n_serotypes_present"] == 1).all()
    assert (ra["conservation_class"] == CONSERVATION_UNIQUE).all()


def test_empty_replicate_inventory_defaults(
    canonical_residues: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    # replicate inventory empty -> every residue defaults to no replicates
    ra = build_residue_annotation(
        canonical_residues, conservation_table, pd.DataFrame(),
        N_SEROTYPES_TOTAL,
    )
    assert (ra["n_replicates"] == 0).all()
    assert (ra["availability_class"] == AVAIL_NONE).all()
