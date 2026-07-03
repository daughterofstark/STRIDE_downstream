"""Tests for the four S3 hierarchy-reduction builders against hand-computed values."""
from __future__ import annotations

import math

import pandas as pd

from stride_s3.build import (
    build_chain_contrast,
    build_monotonicity_audit,
    build_resolution_gap,
    build_scale_curve,
)
from stride_s3.models.schema import (
    CHAIN_CONTRAST_COLUMNS,
    MONOTONICITY_AUDIT_COLUMNS,
    N_SCALES,
    RESOLUTION_GAP_COLUMNS,
    SCALE_CURVE_COLUMNS,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)


# ---------------------------------------------------------------------------
# scale curve
# ---------------------------------------------------------------------------
def test_scale_curve_schema_and_row_count(stride_table: pd.DataFrame) -> None:
    sc = build_scale_curve(stride_table)
    assert list(sc.columns) == list(SCALE_CURVE_COLUMNS)
    # 2 serotypes × 4 loci × 7 scales = 56 rows
    assert len(sc) == 2 * 4 * N_SCALES


def test_scale_curve_gains_and_gating(stride_table: pd.DataFrame) -> None:
    sc = build_scale_curve(stride_table)
    locus = sc[
        (sc.serotype == "DENVA") & (sc.canon_label == "NS3:200")
    ].sort_values("scale_index")
    rhos = locus["rho"].tolist()
    assert rhos == [0.30, 0.35, 0.40, 0.70, 0.75, 0.80, 0.85]
    # residue row: rho_prev / step_gain are NaN
    res = locus.iloc[0]
    assert math.isnan(res.rho_prev)
    assert math.isnan(res.rho_step_gain)
    assert res.rho_cumulative_gain == 0.0
    # domain row: big step gain, gated here
    dom = locus[locus.scale_level == "domain"].iloc[0]
    assert round(dom.rho_step_gain, 6) == round(0.30, 6)
    assert round(dom.rho_cumulative_gain, 6) == round(0.40, 6)
    assert dom.is_gated_scale


def test_scale_curve_tiers(stride_table: pd.DataFrame) -> None:
    sc = build_scale_curve(stride_table)
    assert (sc[sc.scale_index <= 2].tier == TIER_EXPLORATORY).all()
    assert (sc[sc.scale_index >= 3].tier == TIER_LICENSED).all()


def test_scale_curve_exactly_one_gated_scale_per_locus(
    stride_table: pd.DataFrame,
) -> None:
    sc = build_scale_curve(stride_table)
    per_locus = sc.groupby(["serotype", "canon_label"])["is_gated_scale"].sum()
    assert (per_locus == 1).all()


# ---------------------------------------------------------------------------
# resolution gap
# ---------------------------------------------------------------------------
def test_resolution_gap_schema_and_key(stride_table: pd.DataFrame) -> None:
    rg = build_resolution_gap(stride_table)
    assert list(rg.columns) == list(RESOLUTION_GAP_COLUMNS)
    assert len(rg) == 8  # 2 serotypes × 4 loci
    assert not rg.duplicated(["serotype", "canon_label"]).any()


def test_resolution_gap_distributed_flag(stride_table: pd.DataFrame) -> None:
    rg = build_resolution_gap(stride_table).set_index(["serotype", "canon_label"])
    # NS3:200 residue 0.30 < 0.5 <= domain 0.70 → distributed
    assert rg.loc[("DENVA", "NS3:200"), "is_distributed"]
    assert round(rg.loc[("DENVA", "NS3:200"), "delta_rho_domain_residue"], 6) == round(0.40, 6)
    # NS3:51 residue 0.80 already reproducible → not distributed
    assert not rg.loc[("DENVA", "NS3:51"), "is_distributed"]


def test_resolution_gap_min_max_and_gated(stride_table: pd.DataFrame) -> None:
    rg = build_resolution_gap(stride_table).set_index(["serotype", "canon_label"])
    row = rg.loc[("DENVA", "NS3:200")]
    assert row.rho_min == 0.30
    assert row.rho_max == 0.85
    assert row.gated_scale_level == "domain"
    assert row.gated_scale_index == 3
    assert row.rho_at_gated == 0.70


# ---------------------------------------------------------------------------
# monotonicity audit
# ---------------------------------------------------------------------------
def test_monotonicity_audit_schema(stride_table: pd.DataFrame) -> None:
    ma = build_monotonicity_audit(stride_table)
    assert list(ma.columns) == list(MONOTONICITY_AUDIT_COLUMNS)
    assert len(ma) == 8


def test_monotonicity_flags_non_monotone_locus(
    stride_table: pd.DataFrame,
) -> None:
    ma = build_monotonicity_audit(stride_table).set_index(
        ["serotype", "canon_label"]
    )
    # NS3:99 dips at motif → non-monotone with one violation
    ns99 = ma.loc[("DENVA", "NS3:99")]
    assert not ns99.is_monotone
    assert ns99.n_violations == 1
    assert round(ns99.max_decrease, 6) == round(0.05, 6)
    assert ns99.first_violation_scale_index == 1
    # the other loci are monotone
    assert ma.loc[("DENVA", "NS3:51"), "is_monotone"]
    assert ma.loc[("DENVA", "NS3:200"), "is_monotone"]
    assert ma.loc[("DENVA", "NS2B:-1"), "is_monotone"]


def test_monotonicity_context_rho(stride_table: pd.DataFrame) -> None:
    ma = build_monotonicity_audit(stride_table).set_index(
        ["serotype", "canon_label"]
    )
    row = ma.loc[("DENVA", "NS3:51")]
    assert row.rho_residue == 0.80
    assert row.rho_complex == 0.95
    assert row.n_scales == N_SCALES


# ---------------------------------------------------------------------------
# chain contrast
# ---------------------------------------------------------------------------
def test_chain_contrast_schema_and_two_chains(
    stride_table: pd.DataFrame,
) -> None:
    cc = build_chain_contrast(stride_table)
    assert list(cc.columns) == list(CHAIN_CONTRAST_COLUMNS)
    # 2 serotypes × 2 chains (NS2B, NS3) = 4 rows
    assert len(cc) == 4
    assert set(cc.chain.unique()) == {"NS2B", "NS3"}
    assert (cc.tier == TIER_LICENSED).all()


def test_chain_contrast_ns3_aggregates(stride_table: pd.DataFrame) -> None:
    cc = build_chain_contrast(stride_table).set_index(["serotype", "chain"])
    ns3 = cc.loc[("DENVA", "NS3")]
    # NS3 has 3 loci: residue ρ = 0.80, 0.30, 0.42
    assert ns3.n_loci == 3
    assert round(ns3.rho_residue_min, 6) == 0.30
    assert round(ns3.rho_residue_max, 6) == 0.80
    # direction: NS3:51 increase, NS3:200 mixed, NS3:99 decrease
    assert ns3.n_increase == 1
    assert ns3.n_decrease == 1
    assert ns3.n_mixed == 1
    assert ns3.n_mechanisms == 3
    assert ns3.n_signed == 2


def test_chain_contrast_ns2b_present(stride_table: pd.DataFrame) -> None:
    cc = build_chain_contrast(stride_table).set_index(["serotype", "chain"])
    ns2b = cc.loc[("DENVA", "NS2B")]
    assert ns2b.n_loci == 1
    assert ns2b.n_decrease == 1
    assert ns2b.n_mechanisms == 1


def test_chain_contrast_chain_rho_aggregated_not_constant(
    stride_table: pd.DataFrame,
) -> None:
    # NS3 chain-scale ρ differs across loci (0.90, 0.75, 0.78); the builder must
    # aggregate (mean/median) rather than assume a single region-constant value.
    cc = build_chain_contrast(stride_table).set_index(["serotype", "chain"])
    ns3 = cc.loc[("DENVA", "NS3")]
    expected_mean = round((0.90 + 0.75 + 0.78) / 3, 6)
    assert round(ns3.rho_chain_mean, 6) == expected_mean


# ---------------------------------------------------------------------------
# empty input tolerance
# ---------------------------------------------------------------------------
def test_builders_tolerate_empty() -> None:
    empty = pd.DataFrame()
    assert build_scale_curve(empty).empty
    assert build_resolution_gap(empty).empty
    assert build_monotonicity_audit(empty).empty
    assert build_chain_contrast(empty).empty
