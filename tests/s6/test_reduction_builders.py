"""Tests for the four S6 replicate-layer builders."""
from __future__ import annotations

import math

import pandas as pd

from stride_s6.build import (
    build_replicate_blocked_analyses,
    build_replicate_concordance,
    build_replicate_effect_spread,
    build_replicate_regime,
    per_run_effects,
)
from stride_s6.models.schema import (
    ANALYSIS_LORO_STABILITY,
    CONCORDANCE_CLASSES,
    CONCORDANCE_STRONG,
    REPLICATE_BLOCKED_ANALYSES_COLUMNS,
    REPLICATE_CONCORDANCE_COLUMNS,
    REPLICATE_EFFECT_SPREAD_COLUMNS,
    REPLICATE_REGIME_COLUMNS,
    STATUS_BLOCKED,
    TIER_EXPLORATORY,
)
from tests.s6.fixtures import make_licensed_inventory


# --------------------------------------------------------------------------- #
# replicate_regime
# --------------------------------------------------------------------------- #
def test_regime_columns_and_one_row_per_serotype(
    replicate_inventory: pd.DataFrame, effects: pd.DataFrame
) -> None:
    out = build_replicate_regime(replicate_inventory, effects)
    assert list(out.columns) == list(REPLICATE_REGIME_COLUMNS)
    assert list(out["serotype"]) == ["DENV1", "DENV2", "DENV3", "DENV4"]


def test_regime_completeness_and_availability(
    replicate_inventory: pd.DataFrame, effects: pd.DataFrame
) -> None:
    out = build_replicate_regime(replicate_inventory, effects)
    denv1 = out[out["serotype"] == "DENV1"].iloc[0]
    assert denv1["n_positions"] == 6
    assert denv1["frac_complete"] == 1.0
    assert bool(denv1["per_replicate_effects_available"]) is True
    assert denv1["n_replicates_with_effects"] == 3
    # DENV4 loses one position from the complete set → completeness < 1
    denv4 = out[out["serotype"] == "DENV4"].iloc[0]
    assert denv4["n_positions_in_all_replicates"] == 5
    assert denv4["frac_complete"] < 1.0


def test_regime_residue_license_off_at_k3_on_at_k5(
    effects: pd.DataFrame,
) -> None:
    off = build_replicate_regime(make_licensed_inventory(3), effects)
    assert not off["residue_claims_licensed"].any()
    on = build_replicate_regime(make_licensed_inventory(5), effects)
    assert on["residue_claims_licensed"].all()


def test_regime_blocked_when_no_effects(empty_inventory: pd.DataFrame) -> None:
    out = build_replicate_regime(empty_inventory, per_run_effects(None))
    assert not out["per_replicate_effects_available"].any()
    assert (out["n_replicates_with_effects"] == 0).all()


def test_regime_does_not_mutate_inputs(
    replicate_inventory: pd.DataFrame, effects: pd.DataFrame
) -> None:
    before = replicate_inventory.copy()
    build_replicate_regime(replicate_inventory, effects)
    pd.testing.assert_frame_equal(replicate_inventory, before)


# --------------------------------------------------------------------------- #
# replicate_effect_spread
# --------------------------------------------------------------------------- #
def test_spread_columns_and_tier(
    effects: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
    replicate_table: pd.DataFrame,
) -> None:
    out = build_replicate_effect_spread(effects, replicate_inventory, replicate_table)
    assert list(out.columns) == list(REPLICATE_EFFECT_SPREAD_COLUMNS)
    assert (out["tier"] == TIER_EXPLORATORY).all()


def test_spread_statistics_are_ordered(
    effects: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
    replicate_table: pd.DataFrame,
) -> None:
    out = build_replicate_effect_spread(effects, replicate_inventory, replicate_table)
    assert (out["theta_max"] >= out["theta_min"]).all()
    assert (out["theta_range"] >= 0).all()
    assert (out["max_pairwise_abs_diff"] >= 0).all()
    assert (out["abs_theta_mean"] >= 0).all()
    # every default-scenario position is observed (DENV4/NS3:250 in 2 runs still
    # yields a row) → one row per (serotype, position)
    assert len(out) == 4 * 6


def test_spread_empty_when_no_effects(replicate_inventory: pd.DataFrame) -> None:
    out = build_replicate_effect_spread(
        per_run_effects(None), replicate_inventory, None
    )
    assert out.empty
    assert list(out.columns) == list(REPLICATE_EFFECT_SPREAD_COLUMNS)


# --------------------------------------------------------------------------- #
# replicate_concordance
# --------------------------------------------------------------------------- #
def test_concordance_columns_one_row_per_serotype(effects: pd.DataFrame) -> None:
    out = build_replicate_concordance(effects)
    assert list(out.columns) == list(REPLICATE_CONCORDANCE_COLUMNS)
    assert list(out["serotype"]) == ["DENV1", "DENV2", "DENV3", "DENV4"]


def test_concordance_strong_for_order_preserving_runs(
    effects: pd.DataFrame,
) -> None:
    out = build_replicate_concordance(effects)
    denv1 = out[out["serotype"] == "DENV1"].iloc[0]
    assert denv1["n_positions_complete"] == 6
    assert denv1["kendalls_w"] == 1.0
    assert denv1["concordance_class"] == CONCORDANCE_STRONG


def test_concordance_class_tracks_coefficient(effects: pd.DataFrame) -> None:
    out = build_replicate_concordance(effects)
    # DENV3 reverses one run → strictly less concordant than DENV1's perfect W
    denv1_w = out[out["serotype"] == "DENV1"].iloc[0]["kendalls_w"]
    denv3_w = out[out["serotype"] == "DENV3"].iloc[0]["kendalls_w"]
    assert denv3_w < denv1_w
    assert set(out["concordance_class"]).issubset(set(CONCORDANCE_CLASSES))


def test_concordance_empty_when_no_effects() -> None:
    out = build_replicate_concordance(per_run_effects(None))
    assert out.empty


def test_concordance_insufficient_for_two_positions() -> None:
    # only two positions per serotype → below the min for a defined W
    df = pd.DataFrame(
        {
            "serotype": ["DENV1"] * 6,
            "replicate_index": [1, 2, 3, 1, 2, 3],
            "canon_label": ["A", "A", "A", "B", "B", "B"],
            "r": [0.1, 0.2, 0.3, 0.9, 0.8, 0.7],
        }
    )
    out = build_replicate_concordance(per_run_effects(df))
    row = out.iloc[0]
    assert row["concordance_class"] == "insufficient"
    assert math.isnan(row["kendalls_w"])


# --------------------------------------------------------------------------- #
# replicate_blocked_analyses
# --------------------------------------------------------------------------- #
def test_blocked_ledger_columns_and_loro_always_blocked() -> None:
    for available in (True, False):
        out = build_replicate_blocked_analyses(available)
        assert list(out.columns) == list(REPLICATE_BLOCKED_ANALYSES_COLUMNS)
        loro = out[out["analysis_id"] == ANALYSIS_LORO_STABILITY].iloc[0]
        assert loro["status"] == STATUS_BLOCKED
        assert bool(loro["available"]) is False


def test_blocked_ledger_per_run_tracks_availability() -> None:
    avail = build_replicate_blocked_analyses(True)
    per_run = avail[avail["analysis_id"] != ANALYSIS_LORO_STABILITY]
    assert per_run["available"].all()

    blocked = build_replicate_blocked_analyses(False)
    per_run_b = blocked[blocked["analysis_id"] != ANALYSIS_LORO_STABILITY]
    assert not per_run_b["available"].any()


def test_blocked_ledger_all_fields_populated() -> None:
    out = build_replicate_blocked_analyses(False)
    for col in ("description", "reason", "required_input", "design_ref"):
        assert out[col].map(lambda s: bool(str(s).strip())).all()
