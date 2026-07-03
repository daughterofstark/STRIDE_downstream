"""Tests for the four S4 uncertainty builders against hand-computed values."""
from __future__ import annotations

import math

import pandas as pd

from stride_s4.build import (
    build_domain_effect_summary,
    build_residue_variance,
    build_significance_screen,
    build_variance_budget,
)
from stride_s4.models.schema import (
    DOMAIN_EFFECT_SUMMARY_COLUMNS,
    RESIDUE_VARIANCE_COLUMNS,
    SIGNIFICANCE_SCREEN_COLUMNS,
    TIER_EXPLORATORY,
    TIER_LICENSED,
    VARIANCE_BUDGET_COLUMNS,
    VARIANCE_REGIME_BALANCED,
    VARIANCE_REGIME_REPLICATE,
    VARIANCE_REGIME_SAMPLING,
)


# ---------------------------------------------------------------------------
# variance budget
# ---------------------------------------------------------------------------
def test_variance_budget_schema_and_count(stride_table: pd.DataFrame) -> None:
    vb = build_variance_budget(stride_table)
    assert list(vb.columns) == list(VARIANCE_BUDGET_COLUMNS)
    # 2 serotypes × 3 distinct domains = 6 rows
    assert len(vb) == 6
    assert (vb.tier == TIER_LICENSED).all()


def test_variance_budget_regimes(stride_table: pd.DataFrame) -> None:
    vb = build_variance_budget(stride_table).set_index(
        ["serotype", "chain", "domain"]
    )
    # Catalytic Triad: τ²=0.42, σ̄²=0.18 → frac_tau2 = 0.7 → replicate_dominated
    ct = vb.loc[("DENVA", "NS3", "Catalytic Triad")]
    assert abs(ct.frac_tau2 - 0.7) < 1e-9
    assert abs(ct.frac_sigma2 - 0.3) < 1e-9
    assert ct.variance_regime == VARIANCE_REGIME_REPLICATE
    # C-Terminal Tail: τ²=0.10, σ̄²=0.40 → frac_tau2 = 0.2 → sampling_dominated
    ctt = vb.loc[("DENVA", "NS3", "C-Terminal Tail")]
    assert abs(ctt.frac_tau2 - 0.2) < 1e-9
    assert ctt.variance_regime == VARIANCE_REGIME_SAMPLING
    # Cofactor Interface: τ²=0.25, σ̄²=0.25 → balanced
    cof = vb.loc[("DENVA", "NS2B", "Cofactor Interface")]
    assert abs(cof.frac_tau2 - 0.5) < 1e-9
    assert cof.variance_regime == VARIANCE_REGIME_BALANCED


def test_variance_budget_total_and_ratio(stride_table: pd.DataFrame) -> None:
    vb = build_variance_budget(stride_table).set_index(
        ["serotype", "chain", "domain"]
    )
    ct = vb.loc[("DENVA", "NS3", "Catalytic Triad")]
    assert abs(ct.total_unreproduced - 0.6) < 1e-9
    assert abs(ct.tau2_sigma2_ratio - (0.42 / 0.18)) < 1e-6


# ---------------------------------------------------------------------------
# residue variance
# ---------------------------------------------------------------------------
def test_residue_variance_schema_and_count(stride_table: pd.DataFrame) -> None:
    rv = build_residue_variance(stride_table)
    assert list(rv.columns) == list(RESIDUE_VARIANCE_COLUMNS)
    # 2 serotypes × 4 loci = 8 rows
    assert len(rv) == 8
    assert (rv.tier == TIER_EXPLORATORY).all()


def test_residue_variance_tau2_rank(stride_table: pd.DataFrame) -> None:
    rv = build_residue_variance(stride_table)
    denva = rv[rv.serotype == "DENVA"].set_index("canon_label")
    # residue τ²: NS3:51=0.30, NS3:99=0.20, NS3:200=0.05, NS2B:-1=0.25
    # ranks by τ² desc: NS3:51=1, NS2B:-1=2, NS3:99=3, NS3:200=4
    assert denva.loc["NS3:51", "tau2_rank"] == 1
    assert denva.loc["NS2B:-1", "tau2_rank"] == 2
    assert denva.loc["NS3:99", "tau2_rank"] == 3
    assert denva.loc["NS3:200", "tau2_rank"] == 4


def test_residue_variance_rank_is_within_serotype(
    stride_table: pd.DataFrame,
) -> None:
    rv = build_residue_variance(stride_table)
    # each serotype independently has ranks 1..4
    for serotype in ("DENVA", "DENVB"):
        ranks = sorted(rv[rv.serotype == serotype]["tau2_rank"].tolist())
        assert ranks == [1, 2, 3, 4]


# ---------------------------------------------------------------------------
# significance screen
# ---------------------------------------------------------------------------
def test_significance_screen_schema_and_count(
    stride_table: pd.DataFrame,
) -> None:
    ss = build_significance_screen(stride_table)
    assert list(ss.columns) == list(SIGNIFICANCE_SCREEN_COLUMNS)
    # one row per gated mechanism = 4 loci × 2 serotypes = 8
    assert len(ss) == 8


def test_significance_screen_signed_increase(stride_table: pd.DataFrame) -> None:
    ss = build_significance_screen(stride_table).set_index(
        ["serotype", "canon_label"]
    )
    row = ss.loc[("DENVA", "NS3:51")]
    assert row.is_signed
    assert row.direction == "increase"
    assert row.ci_excludes_zero  # [0.02, 0.22]
    assert abs(row.z_score - 2.4) < 1e-6  # 0.12 / 0.05
    assert abs(row.p_value - 0.016395) < 1e-5
    assert row.significant_raw
    # gated at domain → licensed tier
    assert row.tier == TIER_LICENSED


def test_significance_screen_ci_touches_zero(stride_table: pd.DataFrame) -> None:
    ss = build_significance_screen(stride_table).set_index(
        ["serotype", "canon_label"]
    )
    row = ss.loc[("DENVA", "NS3:99")]
    assert row.is_signed
    assert not row.ci_excludes_zero  # [-0.01, 0.21] touches 0
    assert not row.significant_raw
    # gated at residue → exploratory tier
    assert row.tier == TIER_EXPLORATORY


def test_significance_screen_mixed_mechanism(stride_table: pd.DataFrame) -> None:
    ss = build_significance_screen(stride_table).set_index(
        ["serotype", "canon_label"]
    )
    row = ss.loc[("DENVA", "NS2B:-1")]
    assert not row.is_signed
    assert row.direction == "mixed"
    assert not row.ci_excludes_zero
    assert math.isnan(row.p_value)
    assert math.isnan(row.p_value_bh)
    assert not row.significant_fdr


def test_significance_screen_bh_within_serotype(
    stride_table: pd.DataFrame,
) -> None:
    ss = build_significance_screen(stride_table)
    # DENVA signed mechanisms: NS3:51 (p≈0.0164), NS3:200, NS3:99
    denva_signed = ss[(ss.serotype == "DENVA") & (ss.is_signed)]
    # every signed row has a BH p-value; every mixed row does not
    assert denva_signed["p_value_bh"].notna().all()
    mixed = ss[~ss.is_signed]
    assert mixed["p_value_bh"].isna().all()
    # BH adjusted >= raw for each signed row
    for _, r in denva_signed.iterrows():
        assert r.p_value_bh + 1e-9 >= r.p_value


# ---------------------------------------------------------------------------
# domain effect summary
# ---------------------------------------------------------------------------
def test_domain_effect_summary_schema(stride_table: pd.DataFrame) -> None:
    ss = build_significance_screen(stride_table)
    des = build_domain_effect_summary(ss)
    assert list(des.columns) == list(DOMAIN_EFFECT_SUMMARY_COLUMNS)
    # domains with ≥1 mechanism: Catalytic Triad, C-Terminal Tail,
    # Cofactor Interface × 2 serotypes = 6
    assert len(des) == 6
    assert (des.tier == TIER_LICENSED).all()


def test_domain_effect_summary_catalytic_triad(
    stride_table: pd.DataFrame,
) -> None:
    ss = build_significance_screen(stride_table)
    des = build_domain_effect_summary(ss).set_index(
        ["serotype", "chain", "domain"]
    )
    # Catalytic Triad has 2 signed mechanisms (NS3:51 and NS3:99)
    ct = des.loc[("DENVA", "NS3", "Catalytic Triad")]
    assert ct.n_mechanisms == 2
    assert ct.n_signed == 2
    assert ct.n_mixed == 0
    assert ct.n_ci_excludes_zero == 1  # only NS3:51 excludes 0
    assert abs(ct.frac_ci_excludes_zero - 0.5) < 1e-9
    # weighted mean of 0.12 (se 0.05) and 0.10 (se 0.055)
    w1, w2 = 1 / 0.05**2, 1 / 0.055**2
    expected = (w1 * 0.12 + w2 * 0.10) / (w1 + w2)
    assert abs(ct.beta_weighted_mean - round(expected, 6)) < 1e-5
    assert abs(ct.beta_unweighted_mean - 0.11) < 1e-9  # (0.12+0.10)/2


def test_domain_effect_summary_mixed_only_domain(
    stride_table: pd.DataFrame,
) -> None:
    ss = build_significance_screen(stride_table)
    des = build_domain_effect_summary(ss).set_index(
        ["serotype", "chain", "domain"]
    )
    # Cofactor Interface has 1 mixed mechanism, 0 signed
    cof = des.loc[("DENVA", "NS2B", "Cofactor Interface")]
    assert cof.n_mechanisms == 1
    assert cof.n_signed == 0
    assert cof.n_mixed == 1
    assert math.isnan(cof.frac_ci_excludes_zero)
    assert math.isnan(cof.beta_weighted_mean)


# ---------------------------------------------------------------------------
# empty input tolerance
# ---------------------------------------------------------------------------
def test_builders_tolerate_empty() -> None:
    empty = pd.DataFrame()
    assert build_variance_budget(empty).empty
    assert build_residue_variance(empty).empty
    assert build_significance_screen(empty).empty
    assert build_domain_effect_summary(empty).empty
