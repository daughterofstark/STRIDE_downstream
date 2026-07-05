"""Unit tests for the deterministic S5 frame and classification helpers."""
from __future__ import annotations

import math

import pandas as pd
import pytest

from stride_s5.build._classify import (
    concordance_class,
    conservation_class,
    is_catalytic_domain,
    is_catalytic_triad,
    majority_direction,
)
from stride_s5.build._frames import (
    domain_regions,
    gated_directions,
    position_frame,
    residue_slice,
    scale_tier,
)
from stride_s5.models.errors import ConsistencyError
from stride_s5.models.schema import (
    CONCORDANCE_AGREE,
    CONCORDANCE_CONFLICT,
    CONCORDANCE_MAJORITY,
    CONSERVATION_ALL,
    CONSERVATION_MAJORITY,
    CONSERVATION_NONE,
    CONSERVATION_SOME,
    DIRECTION_DECREASE,
    DIRECTION_INCREASE,
    DIRECTION_NONE,
    PROVISIONAL_RHO_STAR,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)


# ---------------------------------------------------------------------------
# catalytic membership
# ---------------------------------------------------------------------------
def test_is_catalytic_triad() -> None:
    assert is_catalytic_triad("NS3:51")
    assert is_catalytic_triad("NS3:75")
    assert is_catalytic_triad("NS3:135")
    assert not is_catalytic_triad("NS3:200")
    assert not is_catalytic_triad("NS2B:-1")


def test_is_catalytic_domain() -> None:
    assert is_catalytic_domain("Catalytic Triad")
    assert is_catalytic_domain("Oxyanion Loop")
    assert not is_catalytic_domain("C-Terminal Tail")
    assert not is_catalytic_domain("Cofactor Interface")


# ---------------------------------------------------------------------------
# conservation_class
# ---------------------------------------------------------------------------
def test_conservation_class_all() -> None:
    assert conservation_class(4, 4) == CONSERVATION_ALL
    assert conservation_class(3, 3) == CONSERVATION_ALL


def test_conservation_class_none() -> None:
    assert conservation_class(0, 4) == CONSERVATION_NONE
    assert conservation_class(0, 0) == CONSERVATION_NONE  # absent everywhere


def test_conservation_class_majority() -> None:
    assert conservation_class(3, 4) == CONSERVATION_MAJORITY
    assert conservation_class(2, 3) == CONSERVATION_MAJORITY


def test_conservation_class_some() -> None:
    # exactly half is "some", not "majority"
    assert conservation_class(2, 4) == CONSERVATION_SOME
    assert conservation_class(1, 4) == CONSERVATION_SOME


# ---------------------------------------------------------------------------
# concordance_class / majority_direction
# ---------------------------------------------------------------------------
def test_concordance_class_agree() -> None:
    assert concordance_class(4, 0) == CONCORDANCE_AGREE
    assert concordance_class(0, 3) == CONCORDANCE_AGREE


def test_concordance_class_conflict() -> None:
    assert concordance_class(2, 2) == CONCORDANCE_CONFLICT
    assert concordance_class(1, 1) == CONCORDANCE_CONFLICT


def test_concordance_class_majority() -> None:
    assert concordance_class(3, 1) == CONCORDANCE_MAJORITY
    assert concordance_class(1, 2) == CONCORDANCE_MAJORITY


def test_majority_direction() -> None:
    assert majority_direction(3, 1) == DIRECTION_INCREASE
    assert majority_direction(1, 2) == DIRECTION_DECREASE
    assert majority_direction(2, 2) == DIRECTION_NONE
    assert majority_direction(0, 0) == DIRECTION_NONE


# ---------------------------------------------------------------------------
# scale_tier
# ---------------------------------------------------------------------------
def test_scale_tier() -> None:
    assert scale_tier("residue") == TIER_EXPLORATORY
    assert scale_tier("secondary_structure") == TIER_EXPLORATORY
    assert scale_tier("domain") == TIER_LICENSED
    assert scale_tier("chain") == TIER_LICENSED
    assert scale_tier("complex") == TIER_LICENSED
    assert scale_tier("nonsense") == TIER_EXPLORATORY


# ---------------------------------------------------------------------------
# residue_slice
# ---------------------------------------------------------------------------
def test_residue_slice_one_row_per_locus(stride_table: pd.DataFrame) -> None:
    res = residue_slice(stride_table)
    assert (res["scale_level"] == "residue").all()
    # unique per (serotype, canon_label)
    assert not res.duplicated(["serotype", "canon_label"]).any()


def test_residue_slice_rejects_duplicates(stride_table: pd.DataFrame) -> None:
    dup = pd.concat(
        [stride_table, stride_table[stride_table["scale_level"] == "residue"]],
        ignore_index=True,
    )
    with pytest.raises(ConsistencyError, match="duplicate residue-scale"):
        residue_slice(dup)


# ---------------------------------------------------------------------------
# gated_directions
# ---------------------------------------------------------------------------
def test_gated_directions_one_per_locus(stride_table: pd.DataFrame) -> None:
    g = gated_directions(stride_table)
    assert not g.duplicated(["serotype", "canon_label"]).any()
    # NS3:51 gates at residue with an increase in every serotype
    ns351 = g[g["canon_label"] == "NS3:51"]
    assert (ns351["gated_scale_level"] == "residue").all()
    assert (ns351["direction"] == DIRECTION_INCREASE).all()
    # NS3:200 is not reproducible → gates coarser than residue
    ns3200 = g[g["canon_label"] == "NS3:200"]
    assert (ns3200["gated_scale_level"] != "residue").all()


def test_gated_directions_rejects_multiple(stride_table: pd.DataFrame) -> None:
    extra = stride_table[stride_table["is_gated_scale"]].copy()
    dup = pd.concat([stride_table, extra], ignore_index=True)
    with pytest.raises(ConsistencyError, match="multiple gated rows"):
        gated_directions(dup)


# ---------------------------------------------------------------------------
# position_frame
# ---------------------------------------------------------------------------
def test_position_frame_reproducible_and_signed(
    stride_table: pd.DataFrame,
) -> None:
    pf = position_frame(stride_table, PROVISIONAL_RHO_STAR)
    # NS3:51 reproducible + signed increase everywhere
    ns351 = pf[pf["canon_label"] == "NS3:51"]
    assert ns351["reproducible"].all()
    assert ns351["is_signed"].all()
    assert (ns351["direction"] == DIRECTION_INCREASE).all()
    # NS3:200 reproducible in nobody
    ns3200 = pf[pf["canon_label"] == "NS3:200"]
    assert not ns3200["reproducible"].any()
    assert not ns3200["is_signed"].any()
    assert (ns3200["direction"] == DIRECTION_NONE).all()


def test_position_frame_reproducible_but_mixed(
    stride_table: pd.DataFrame,
) -> None:
    pf = position_frame(stride_table, PROVISIONAL_RHO_STAR)
    # NS2B:-1 in DENV2 is reproducible (ρ=0.6) but mixed → not signed
    row = pf[(pf["canon_label"] == "NS2B:-1") & (pf["serotype"] == "DENV2")]
    assert bool(row["reproducible"].iloc[0]) is True
    assert bool(row["is_signed"].iloc[0]) is False


def test_position_frame_threshold_moves_with_rho_star(
    stride_table: pd.DataFrame,
) -> None:
    # at ρ*=0.65, NS3:135 in DENV3 (ρ=0.62) is no longer reproducible
    pf = position_frame(stride_table, 0.65)
    row = pf[(pf["canon_label"] == "NS3:135") & (pf["serotype"] == "DENV3")]
    assert bool(row["reproducible"].iloc[0]) is False


# ---------------------------------------------------------------------------
# domain_regions
# ---------------------------------------------------------------------------
def test_domain_regions_region_constant(stride_table: pd.DataFrame) -> None:
    reg = domain_regions(stride_table)
    assert not reg.duplicated(["serotype", "chain", "domain"]).any()
    for val in reg["rho"]:
        assert 0.0 <= float(val) <= 1.0


def test_domain_regions_rejects_nonconstant(
    stride_table: pd.DataFrame,
) -> None:
    tampered = stride_table.copy()
    mask = (
        (tampered["serotype"] == "DENV1")
        & (tampered["h_domain"] == "Catalytic Triad")
        & (tampered["scale_level"] == "domain")
    )
    idx = tampered[mask].index[0]
    tampered.loc[idx, "rho"] = 0.123456  # break region-constancy
    with pytest.raises(ConsistencyError, match="non-constant rho"):
        domain_regions(tampered)


def test_helpers_do_not_mutate_input(stride_table: pd.DataFrame) -> None:
    before = stride_table.copy()
    residue_slice(stride_table)
    gated_directions(stride_table)
    position_frame(stride_table, PROVISIONAL_RHO_STAR)
    domain_regions(stride_table)
    pd.testing.assert_frame_equal(stride_table, before)


def test_median_even_and_odd() -> None:
    from stride_s5.build.position_conservation import _median

    assert _median([0.2, 0.8]) == 0.5
    assert _median([0.1, 0.5, 0.9]) == 0.5
    assert math.isnan(_median([]))
