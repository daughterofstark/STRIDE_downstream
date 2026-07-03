"""Tests for the five S2 reduction builders against hand-computed expectations."""
from __future__ import annotations

import pandas as pd

from stride_s2.build import (
    build_domain_reproducibility,
    build_residue_landscape,
    build_resolution_census,
    build_serotype_summary,
    build_signed_screen,
)
from stride_s2.models.schema import (
    DOMAIN_REPRODUCIBILITY_COLUMNS,
    RESIDUE_LANDSCAPE_COLUMNS,
    RESOLUTION_CENSUS_COLUMNS,
    SEROTYPE_SUMMARY_COLUMNS,
    SIGNED_SCREEN_COLUMNS,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)

_BAND = (0.5, 0.7, 0.9)
_DECIMALS = 4


# ---------------------------------------------------------------------------
# resolution census
# ---------------------------------------------------------------------------
def test_census_schema_and_totals(stride_table: pd.DataFrame) -> None:
    census = build_resolution_census(stride_table, _BAND, _DECIMALS)
    assert list(census.columns) == list(RESOLUTION_CENSUS_COLUMNS)
    # at every (serotype, ρ*) the counts sum to the 3 loci
    totals = census.groupby(["serotype", "rho_star"])["n_loci"].sum()
    assert (totals == 3).all()


def test_census_regate_at_provisional(stride_table: pd.DataFrame) -> None:
    census = build_resolution_census(stride_table, _BAND, _DECIMALS)
    at_prov = census[
        (census.serotype == "DENVA") & (census.rho_star == 0.5)
    ].set_index("gated_scale_level")["n_loci"].to_dict()
    # NS3:51 → residue; NS3:200 and NS3:99 → domain (finest ρ≥0.5 is domain)
    assert at_prov == {"residue": 1, "domain": 2}


def test_census_unresolved_at_high_rho_star(stride_table: pd.DataFrame) -> None:
    census = build_resolution_census(stride_table, _BAND, _DECIMALS)
    at_high = census[
        (census.serotype == "DENVA") & (census.rho_star == 0.9)
    ].set_index("gated_scale_level")["n_loci"].to_dict()
    # NS3:51 → chain (0.90); NS3:200 → unresolved (max 0.85); NS3:99 → unresolved
    assert at_high.get("unresolved") == 2
    assert at_high.get("chain") == 1


def test_census_tier_labels(stride_table: pd.DataFrame) -> None:
    census = build_resolution_census(stride_table, _BAND, _DECIMALS)
    residue_rows = census[census.gated_scale_level == "residue"]
    assert (residue_rows.tier == TIER_EXPLORATORY).all()
    domain_rows = census[census.gated_scale_level == "domain"]
    assert (domain_rows.tier == TIER_LICENSED).all()


# ---------------------------------------------------------------------------
# residue landscape
# ---------------------------------------------------------------------------
def test_residue_landscape_schema_and_key(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    rl = build_residue_landscape(stride_table, residue_annotation)
    assert list(rl.columns) == list(RESIDUE_LANDSCAPE_COLUMNS)
    # 2 serotypes × 3 loci = 6 rows, unique key
    assert len(rl) == 6
    assert not rl.duplicated(["serotype", "canon_label"]).any()


def test_residue_landscape_values(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    rl = build_residue_landscape(stride_table, residue_annotation)
    ns51 = rl[(rl.serotype == "DENVA") & (rl.canon_label == "NS3:51")].iloc[0]
    assert ns51.rho_residue == 0.80
    assert ns51.gates_at_residue_provisional
    assert ns51.gated_scale_level_provisional == "residue"
    assert ns51.tier == TIER_EXPLORATORY

    ns200 = rl[(rl.serotype == "DENVA") & (rl.canon_label == "NS3:200")].iloc[0]
    assert ns200.rho_residue == 0.30
    assert not ns200.gates_at_residue_provisional
    assert ns200.gated_scale_level_provisional == "domain"


def test_residue_landscape_is_all_exploratory(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    rl = build_residue_landscape(stride_table, residue_annotation)
    assert (rl.tier == TIER_EXPLORATORY).all()


def test_residue_landscape_falls_back_without_annotation(
    stride_table: pd.DataFrame,
) -> None:
    # empty annotation → chain/domain fall back to the STRIDE hierarchy columns
    rl = build_residue_landscape(stride_table, pd.DataFrame())
    assert len(rl) == 6
    ns51 = rl[(rl.serotype == "DENVA") & (rl.canon_label == "NS3:51")].iloc[0]
    assert ns51.chain == "NS3"
    assert ns51.domain == "Catalytic Triad"


# ---------------------------------------------------------------------------
# domain reproducibility
# ---------------------------------------------------------------------------
def test_domain_reproducibility_schema(
    stride_table: pd.DataFrame, domain_annotation: pd.DataFrame
) -> None:
    dr = build_domain_reproducibility(stride_table, domain_annotation)
    assert list(dr.columns) == list(DOMAIN_REPRODUCIBILITY_COLUMNS)
    # 2 serotypes × 3 domains (Catalytic Triad, C-Terminal Tail, Gly45 Turn) = 6
    assert len(dr) == 6
    assert (dr.tier == TIER_LICENSED).all()


def test_domain_reproducibility_values(
    stride_table: pd.DataFrame, domain_annotation: pd.DataFrame
) -> None:
    dr = build_domain_reproducibility(stride_table, domain_annotation)
    ct = dr[
        (dr.serotype == "DENVA") & (dr.domain == "Catalytic Triad")
    ].iloc[0]
    assert ct.rho_domain == 0.88  # NS3:51's domain-scale ρ
    assert ct.is_coherent  # coherence 0.9 ≥ 0.6

    ctail = dr[
        (dr.serotype == "DENVA") & (dr.domain == "C-Terminal Tail")
    ].iloc[0]
    assert ctail.rho_domain == 0.70
    assert not ctail.is_coherent  # coherence 0.5 < 0.6


def test_domain_reproducibility_detects_nonconstant_region(
    stride_table: pd.DataFrame, domain_annotation: pd.DataFrame
) -> None:
    # Force two loci in the SAME domain to carry different domain-scale ρ: this
    # is a genuine region-constant violation the builder must catch.
    import pytest

    from stride_s2.models.errors import ConsistencyError

    broken = stride_table.copy()
    # relabel NS3:99's domain to "Catalytic Triad" so it collides with NS3:51,
    # which has a different domain-scale ρ (0.88 vs 0.55)
    mask = broken["canon_label"] == "NS3:99"
    broken.loc[mask, "h_domain"] = "Catalytic Triad"

    with pytest.raises(ConsistencyError, match="non-constant"):
        build_domain_reproducibility(broken, domain_annotation)


# ---------------------------------------------------------------------------
# signed screen
# ---------------------------------------------------------------------------
def test_signed_screen_schema(stride_table: pd.DataFrame) -> None:
    ss = build_signed_screen(stride_table, _BAND, _DECIMALS)
    assert list(ss.columns) == list(SIGNED_SCREEN_COLUMNS)
    # 2 serotypes × 3 mechanisms × 3 ρ* = 18 rows
    assert len(ss) == 18


def test_signed_screen_increase_passes_below_rho(
    stride_table: pd.DataFrame,
) -> None:
    ss = build_signed_screen(stride_table, _BAND, _DECIMALS)
    ns51 = ss[(ss.serotype == "DENVA") & (ss.canon_label == "NS3:51")]
    # ρ=0.80: passes at ρ*=0.5 and 0.7, fails at 0.9
    passing = ns51.set_index("rho_star")["passes_screen"].to_dict()
    assert passing == {0.5: True, 0.7: True, 0.9: False}


def test_signed_screen_mixed_never_passes(stride_table: pd.DataFrame) -> None:
    ss = build_signed_screen(stride_table, _BAND, _DECIMALS)
    ns200 = ss[(ss.serotype == "DENVA") & (ss.canon_label == "NS3:200")]
    assert not ns200.passes_screen.any()
    assert not ns200.is_signed.any()


def test_signed_screen_ci_touching_zero_never_passes(
    stride_table: pd.DataFrame,
) -> None:
    ss = build_signed_screen(stride_table, _BAND, _DECIMALS)
    ns99 = ss[(ss.serotype == "DENVA") & (ss.canon_label == "NS3:99")]
    assert ns99.is_signed.all()  # direction is "decrease"
    assert not ns99.ci_excludes_zero.any()  # CI touches 0
    assert not ns99.passes_screen.any()


# ---------------------------------------------------------------------------
# serotype summary
# ---------------------------------------------------------------------------
def _summary(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> pd.DataFrame:
    census = build_resolution_census(stride_table, _BAND, _DECIMALS)
    rl = build_residue_landscape(stride_table, residue_annotation)
    ss = build_signed_screen(stride_table, _BAND, _DECIMALS)
    return build_serotype_summary(
        stride_table, census, rl, ss, _BAND, _DECIMALS
    )


def test_serotype_summary_schema(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    sm = _summary(stride_table, residue_annotation)
    assert list(sm.columns) == list(SEROTYPE_SUMMARY_COLUMNS)
    # 2 serotypes × 3 ρ* = 6 rows
    assert len(sm) == 6


def test_serotype_summary_direction_partition(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    sm = _summary(stride_table, residue_annotation)
    for row in sm.itertuples(index=False):
        assert row.n_signed + row.n_mixed == row.n_mechanisms
    # 3 mechanisms: NS3:51 (increase), NS3:99 (decrease), NS3:200 (mixed)
    a = sm[sm.serotype == "DENVA"].iloc[0]
    assert a.n_mechanisms == 3
    assert a.n_signed == 2
    assert a.n_mixed == 1


def test_serotype_summary_signed_significant_over_band(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    sm = _summary(stride_table, residue_annotation).set_index(
        ["serotype", "rho_star"]
    )
    # only NS3:51 passes (NS3:99 CI touches 0, NS3:200 mixed):
    # 1 at ρ*≤0.8, 0 at ρ*=0.9
    assert sm.loc[("DENVA", 0.5), "n_signed_significant"] == 1
    assert sm.loc[("DENVA", 0.7), "n_signed_significant"] == 1
    assert sm.loc[("DENVA", 0.9), "n_signed_significant"] == 0


def test_serotype_summary_rho_distribution(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    sm = _summary(stride_table, residue_annotation)
    a = sm[sm.serotype == "DENVA"].iloc[0]
    # residue-scale ρ are 0.80, 0.30, 0.20 → median 0.30
    assert a.rho_residue_median == 0.30
    assert a.rho_residue_min == 0.20
    assert a.rho_residue_max == 0.80
