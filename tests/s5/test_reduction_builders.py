"""Tests for the four S5 cross-serotype builders against hand-computed values."""
from __future__ import annotations

import pandas as pd

from stride_s5.build import (
    build_cross_serotype_scorecard,
    build_direction_concordance,
    build_domain_serotype_matrix,
    build_position_conservation,
)
from stride_s5.models.schema import (
    CONCORDANCE_AGREE,
    CONCORDANCE_CONFLICT,
    CONCORDANCE_MAJORITY,
    CONSERVATION_ALL,
    CONSERVATION_MAJORITY,
    CONSERVATION_NONE,
    CONSERVATION_SOME,
    CROSS_SEROTYPE_SCORECARD_COLUMNS,
    DIRECTION_CONCORDANCE_COLUMNS,
    DIRECTION_DECREASE,
    DIRECTION_INCREASE,
    DOMAIN_SEROTYPE_MATRIX_COLUMNS,
    POSITION_CONSERVATION_COLUMNS,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)


def _row(df: pd.DataFrame, canon: str) -> pd.Series:
    return df[df["canon_label"] == canon].iloc[0]


# ---------------------------------------------------------------------------
# position_conservation
# ---------------------------------------------------------------------------
def test_position_conservation_schema_and_key(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc = build_position_conservation(stride_table, conservation_table)
    assert list(pc.columns) == list(POSITION_CONSERVATION_COLUMNS)
    assert not pc.duplicated(["canon_label"]).any()
    assert len(pc) == 6  # one per shared position
    assert (pc["tier"] == TIER_EXPLORATORY).all()
    assert (pc["is_provisional_rho_star"]).all()


def test_position_conservation_classes(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc = build_position_conservation(stride_table, conservation_table)
    assert _row(pc, "NS3:51")["conservation_class"] == CONSERVATION_ALL
    assert _row(pc, "NS3:75")["conservation_class"] == CONSERVATION_ALL
    assert _row(pc, "NS3:135")["conservation_class"] == CONSERVATION_ALL
    assert _row(pc, "NS3:200")["conservation_class"] == CONSERVATION_NONE
    assert _row(pc, "NS2B:-1")["conservation_class"] == CONSERVATION_MAJORITY
    assert _row(pc, "NS3:250")["conservation_class"] == CONSERVATION_SOME


def test_position_conservation_counts_and_flags(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc = build_position_conservation(stride_table, conservation_table)
    ns2b = _row(pc, "NS2B:-1")
    assert int(ns2b["n_serotypes_present"]) == 3
    assert int(ns2b["n_serotypes_reproducible"]) == 2
    assert int(ns2b["n_serotypes_signed_reproducible"]) == 1
    assert bool(ns2b["is_serotype_divergent"]) is True
    assert bool(ns2b["in_all_serotypes"]) is False
    # catalytic triad flags
    assert bool(_row(pc, "NS3:51")["is_catalytic_triad"]) is True
    assert bool(_row(pc, "NS3:200")["is_catalytic_triad"]) is False
    # divergent only for the two positions signed in some-but-not-all serotypes
    divergent = set(pc[pc["is_serotype_divergent"]]["canon_label"])
    assert divergent == {"NS2B:-1", "NS3:250"}


def test_position_conservation_rho_stats(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc = build_position_conservation(stride_table, conservation_table)
    ns351 = _row(pc, "NS3:51")
    # ρ values 0.80, 0.82, 0.78, 0.85 → min 0.78, max 0.85, median 0.81
    assert abs(float(ns351["rho_residue_min"]) - 0.78) < 1e-9
    assert abs(float(ns351["rho_residue_max"]) - 0.85) < 1e-9
    assert abs(float(ns351["rho_residue_median"]) - 0.81) < 1e-9


# ---------------------------------------------------------------------------
# direction_concordance
# ---------------------------------------------------------------------------
def test_direction_concordance_schema_and_key(
    stride_table: pd.DataFrame,
) -> None:
    dc = build_direction_concordance(stride_table)
    assert list(dc.columns) == list(DIRECTION_CONCORDANCE_COLUMNS)
    assert not dc.duplicated(["canon_label"]).any()
    # only positions signed in ≥2 serotypes: NS3:51/75/135/250 (NS2B:-1 excluded)
    assert set(dc["canon_label"]) == {"NS3:51", "NS3:75", "NS3:135", "NS3:250"}


def test_direction_concordance_classes(stride_table: pd.DataFrame) -> None:
    dc = build_direction_concordance(stride_table)
    assert _row(dc, "NS3:51")["concordance_class"] == CONCORDANCE_AGREE
    assert _row(dc, "NS3:51")["majority_direction"] == DIRECTION_INCREASE
    assert _row(dc, "NS3:75")["concordance_class"] == CONCORDANCE_CONFLICT
    assert _row(dc, "NS3:135")["concordance_class"] == CONCORDANCE_MAJORITY
    assert _row(dc, "NS3:135")["majority_direction"] == DIRECTION_INCREASE
    dc250 = _row(dc, "NS3:250")
    assert dc250["concordance_class"] == CONCORDANCE_AGREE
    assert dc250["majority_direction"] == DIRECTION_DECREASE
    assert int(dc250["n_decrease"]) == 2 and int(dc250["n_increase"]) == 0


def test_direction_concordance_counts_partition(
    stride_table: pd.DataFrame,
) -> None:
    dc = build_direction_concordance(stride_table)
    for row in dc.itertuples(index=False):
        assert int(row.n_increase) + int(row.n_decrease) == int(
            row.n_serotypes_signed
        )


# ---------------------------------------------------------------------------
# domain_serotype_matrix
# ---------------------------------------------------------------------------
def test_domain_serotype_matrix_schema_and_key(
    stride_table: pd.DataFrame,
) -> None:
    m = build_domain_serotype_matrix(stride_table)
    assert list(m.columns) == list(DOMAIN_SEROTYPE_MATRIX_COLUMNS)
    assert not m.duplicated(["serotype", "chain", "domain"]).any()
    assert (m["tier"] == TIER_LICENSED).all()
    # 15 tidy-long cells: CatTriad×4 + C-Term×4 + Cofactor×3 + Oxyanion×4
    assert len(m) == 15


def test_domain_serotype_matrix_catalytic_flag(
    stride_table: pd.DataFrame,
) -> None:
    m = build_domain_serotype_matrix(stride_table)
    assert int(m["is_catalytic_domain"].sum()) == 8  # Triad×4 + Oxyanion×4
    cat = set(m[m["is_catalytic_domain"]]["domain"])
    assert cat == {"Catalytic Triad", "Oxyanion Loop"}
    # ρ is read from the domain-scale region-constant row
    triad = m[(m["serotype"] == "DENV1") & (m["domain"] == "Catalytic Triad")]
    assert abs(float(triad["rho_domain"].iloc[0]) - 0.88) < 1e-9


# ---------------------------------------------------------------------------
# cross_serotype_scorecard
# ---------------------------------------------------------------------------
def test_scorecard_schema_and_key(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    sc = build_cross_serotype_scorecard(stride_table, conservation_table)
    assert list(sc.columns) == list(CROSS_SEROTYPE_SCORECARD_COLUMNS)
    assert not sc.duplicated(["serotype"]).any()
    assert len(sc) == 4
    assert (sc["tier"] == TIER_EXPLORATORY).all()


def test_scorecard_values(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    sc = build_cross_serotype_scorecard(stride_table, conservation_table)
    d1 = sc[sc["serotype"] == "DENV1"].iloc[0]
    assert int(d1["n_loci"]) == 6
    assert int(d1["n_reproducible_residue"]) == 5
    assert int(d1["n_mechanisms"]) == 6
    assert int(d1["n_signed"]) == 5
    assert int(d1["n_mixed"]) == 1
    assert int(d1["n_shared_positions"]) == 5
    assert int(d1["n_shared_reproducible"]) == 4
    # DENV4 lacks NS2B:-1 → only 5 loci
    d4 = sc[sc["serotype"] == "DENV4"].iloc[0]
    assert int(d4["n_loci"]) == 5


def test_scorecard_partition(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    sc = build_cross_serotype_scorecard(stride_table, conservation_table)
    for row in sc.itertuples(index=False):
        assert int(row.n_signed) + int(row.n_mixed) == int(row.n_mechanisms)
        assert int(row.n_reproducible_residue) <= int(row.n_loci)
        assert int(row.n_shared_reproducible) <= int(row.n_shared_positions)


# ---------------------------------------------------------------------------
# purity / determinism at the builder level
# ---------------------------------------------------------------------------
def test_builders_do_not_mutate_inputs(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    st_before = stride_table.copy()
    ct_before = conservation_table.copy()
    build_position_conservation(stride_table, conservation_table)
    build_direction_concordance(stride_table)
    build_domain_serotype_matrix(stride_table)
    build_cross_serotype_scorecard(stride_table, conservation_table)
    pd.testing.assert_frame_equal(stride_table, st_before)
    pd.testing.assert_frame_equal(conservation_table, ct_before)


def test_builders_on_empty_inputs() -> None:
    from stride_s5.models.schema import (
        CONSERVATION_TABLE_REQUIRED,
        STRIDE_TABLE_REQUIRED,
    )

    empty_st = pd.DataFrame(columns=list(STRIDE_TABLE_REQUIRED))
    empty_ct = pd.DataFrame(columns=list(CONSERVATION_TABLE_REQUIRED))
    assert build_position_conservation(empty_st, empty_ct).empty
    assert build_direction_concordance(empty_st).empty
    assert build_domain_serotype_matrix(empty_st).empty
    assert build_cross_serotype_scorecard(empty_st, empty_ct).empty
