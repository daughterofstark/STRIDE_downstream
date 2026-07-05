"""Tests that the S6 structural validators accept clean tables and catch tampering."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s6.build import (
    build_replicate_blocked_analyses,
    build_replicate_concordance,
    build_replicate_effect_spread,
    build_replicate_regime,
)
from stride_s6.models import S6Report
from stride_s6.models.errors import ConsistencyError
from stride_s6.validation import (
    validate_all,
    validate_replicate_blocked_analyses,
    validate_replicate_concordance,
    validate_replicate_effect_spread,
    validate_replicate_regime,
    validate_unique_keys,
)

_Tables = tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]


@pytest.fixture
def tables(
    replicate_inventory: pd.DataFrame,
    replicate_table: pd.DataFrame,
    effects: pd.DataFrame,
) -> _Tables:
    regime = build_replicate_regime(replicate_inventory, effects)
    spread = build_replicate_effect_spread(
        effects, replicate_inventory, replicate_table
    )
    concordance = build_replicate_concordance(effects)
    blocked = build_replicate_blocked_analyses(True)
    return regime, spread, concordance, blocked


def test_validate_all_passes_on_clean_tables(tables: _Tables) -> None:
    regime, spread, concordance, blocked = tables
    report = S6Report()
    validate_all(regime, spread, concordance, blocked, report)
    assert report.all_passed
    assert len(report.checks) == 5


def test_duplicate_key_rejected(tables: _Tables) -> None:
    regime, spread, concordance, blocked = tables
    dup = pd.concat([regime, regime.head(1)], ignore_index=True)
    with pytest.raises(ConsistencyError):
        validate_unique_keys(dup, spread, concordance, blocked, S6Report())


def test_regime_bad_license_flag_rejected(tables: _Tables) -> None:
    regime, *_ = tables
    bad = regime.copy()
    bad.loc[0, "residue_claims_licensed"] = True  # K=3 → must be False
    with pytest.raises(ConsistencyError):
        validate_replicate_regime(bad, S6Report())


def test_regime_bad_availability_flag_rejected(tables: _Tables) -> None:
    regime, *_ = tables
    bad = regime.copy()
    bad.loc[0, "n_replicates_with_effects"] = 0  # contradicts available=True
    with pytest.raises(ConsistencyError):
        validate_replicate_regime(bad, S6Report())


def test_spread_bad_order_rejected(tables: _Tables) -> None:
    _, spread, *_ = tables
    bad = spread.copy()
    bad.loc[0, "theta_min"] = bad.loc[0, "theta_max"] + 1.0
    with pytest.raises(ConsistencyError):
        validate_replicate_effect_spread(bad, S6Report())


def test_concordance_out_of_range_w_rejected(tables: _Tables) -> None:
    _, _, concordance, _ = tables
    bad = concordance.copy()
    bad.loc[0, "kendalls_w"] = 1.5
    bad.loc[0, "concordance_class"] = "strong"
    with pytest.raises(ConsistencyError):
        validate_replicate_concordance(bad, S6Report())


def test_concordance_class_mismatch_rejected(tables: _Tables) -> None:
    _, _, concordance, _ = tables
    bad = concordance.copy()
    # DENV1 has W=1.0 (strong); relabel to weak → mismatch
    bad.loc[0, "concordance_class"] = "weak"
    with pytest.raises(ConsistencyError):
        validate_replicate_concordance(bad, S6Report())


def test_blocked_ledger_loro_unblocked_rejected(tables: _Tables) -> None:
    *_, blocked = tables
    bad = blocked.copy()
    row = bad["analysis_id"] == "leave_one_replicate_out_stability"
    bad.loc[row, "status"] = "available"
    bad.loc[row, "available"] = True
    with pytest.raises(ConsistencyError):
        validate_replicate_blocked_analyses(bad, S6Report())


def test_blocked_ledger_status_flag_mismatch_rejected(tables: _Tables) -> None:
    *_, blocked = tables
    bad = blocked.copy()
    bad.loc[0, "available"] = not bool(bad.loc[0, "available"])
    with pytest.raises(ConsistencyError):
        validate_replicate_blocked_analyses(bad, S6Report())
