"""Tests for the S5 structural validation checks."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s5.build import (
    build_cross_serotype_scorecard,
    build_direction_concordance,
    build_domain_serotype_matrix,
    build_position_conservation,
)
from stride_s5.models import S5Report
from stride_s5.models.errors import ConsistencyError
from stride_s5.validation import (
    validate_cross_serotype_scorecard,
    validate_direction_concordance,
    validate_domain_serotype_matrix,
    validate_position_conservation,
    validate_unique_keys,
)


def _tables(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pc = build_position_conservation(stride_table, conservation_table)
    dc = build_direction_concordance(stride_table)
    m = build_domain_serotype_matrix(stride_table)
    sc = build_cross_serotype_scorecard(stride_table, conservation_table)
    return pc, dc, m, sc


# ---------------------------------------------------------------------------
# happy path
# ---------------------------------------------------------------------------
def test_all_checks_pass_on_valid_tables(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc, dc, m, sc = _tables(stride_table, conservation_table)
    report = S5Report()
    validate_unique_keys(pc, dc, m, sc, report)
    validate_position_conservation(pc, report)
    validate_direction_concordance(dc, report)
    validate_domain_serotype_matrix(m, report)
    validate_cross_serotype_scorecard(sc, report)
    assert report.all_passed
    assert len(report.checks) == 5


# ---------------------------------------------------------------------------
# unique keys
# ---------------------------------------------------------------------------
def test_duplicate_key_detected(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc, dc, m, sc = _tables(stride_table, conservation_table)
    dup = pd.concat([pc, pc.head(1)], ignore_index=True)
    with pytest.raises(ConsistencyError, match="not unique"):
        validate_unique_keys(dup, dc, m, sc, S5Report())


# ---------------------------------------------------------------------------
# position_conservation
# ---------------------------------------------------------------------------
def test_position_conservation_counts_unordered(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc, *_ = _tables(stride_table, conservation_table)
    pc = pc.copy()
    # signed > reproducible violates 0 <= signed <= reproducible <= present
    pc.loc[pc["canon_label"] == "NS3:200", "n_serotypes_signed_reproducible"] = 2
    with pytest.raises(ConsistencyError, match="counts not ordered"):
        validate_position_conservation(pc, S5Report())


def test_position_conservation_bad_class(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc, *_ = _tables(stride_table, conservation_table)
    pc = pc.copy()
    pc.loc[pc["canon_label"] == "NS3:51", "conservation_class"] = "reproducible_none"
    with pytest.raises(ConsistencyError, match="conservation_class"):
        validate_position_conservation(pc, S5Report())


def test_position_conservation_bad_divergent(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc, *_ = _tables(stride_table, conservation_table)
    pc = pc.copy()
    pc.loc[pc["canon_label"] == "NS3:51", "is_serotype_divergent"] = True
    with pytest.raises(ConsistencyError, match="is_serotype_divergent"):
        validate_position_conservation(pc, S5Report())


def test_position_conservation_bad_catalytic_flag(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc, *_ = _tables(stride_table, conservation_table)
    pc = pc.copy()
    pc.loc[pc["canon_label"] == "NS3:200", "is_catalytic_triad"] = True
    with pytest.raises(ConsistencyError, match="is_catalytic_triad"):
        validate_position_conservation(pc, S5Report())


def test_position_conservation_bad_fraction(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    pc, *_ = _tables(stride_table, conservation_table)
    pc = pc.copy()
    pc.loc[pc["canon_label"] == "NS3:51", "frac_reproducible"] = 0.123
    with pytest.raises(ConsistencyError, match="frac_reproducible"):
        validate_position_conservation(pc, S5Report())


# ---------------------------------------------------------------------------
# direction_concordance
# ---------------------------------------------------------------------------
def test_direction_concordance_counts_not_partition(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    _, dc, *_ = _tables(stride_table, conservation_table)
    dc = dc.copy()
    dc.loc[dc["canon_label"] == "NS3:51", "n_increase"] = 1  # 1 + 0 != 4
    with pytest.raises(ConsistencyError, match="!= n_serotypes_signed"):
        validate_direction_concordance(dc, S5Report())


def test_direction_concordance_below_min(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    _, dc, *_ = _tables(stride_table, conservation_table)
    dc = dc.copy()
    dc.loc[dc["canon_label"] == "NS3:250", "n_serotypes_signed"] = 1
    dc.loc[dc["canon_label"] == "NS3:250", "n_decrease"] = 1
    with pytest.raises(ConsistencyError, match="below the minimum"):
        validate_direction_concordance(dc, S5Report())


def test_direction_concordance_bad_class(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    _, dc, *_ = _tables(stride_table, conservation_table)
    dc = dc.copy()
    dc.loc[dc["canon_label"] == "NS3:51", "concordance_class"] = "conflict"
    with pytest.raises(ConsistencyError, match="concordance_class"):
        validate_direction_concordance(dc, S5Report())


def test_direction_concordance_bad_majority(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    _, dc, *_ = _tables(stride_table, conservation_table)
    dc = dc.copy()
    dc.loc[dc["canon_label"] == "NS3:51", "majority_direction"] = "decrease"
    with pytest.raises(ConsistencyError, match="majority_direction"):
        validate_direction_concordance(dc, S5Report())


# ---------------------------------------------------------------------------
# domain_serotype_matrix
# ---------------------------------------------------------------------------
def test_matrix_rho_out_of_range(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    _, _, m, _ = _tables(stride_table, conservation_table)
    m = m.copy()
    m.loc[0, "rho_domain"] = 1.5
    with pytest.raises(ConsistencyError, match="rho_domain"):
        validate_domain_serotype_matrix(m, S5Report())


def test_matrix_bad_tier(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    _, _, m, _ = _tables(stride_table, conservation_table)
    m = m.copy()
    m.loc[0, "tier"] = "exploratory"
    with pytest.raises(ConsistencyError, match="is not 'licensed'"):
        validate_domain_serotype_matrix(m, S5Report())


def test_matrix_bad_catalytic_flag(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    _, _, m, _ = _tables(stride_table, conservation_table)
    m = m.copy()
    idx = m[m["domain"] == "C-Terminal Tail"].index[0]
    m.loc[idx, "is_catalytic_domain"] = True
    with pytest.raises(ConsistencyError, match="is_catalytic_domain"):
        validate_domain_serotype_matrix(m, S5Report())


# ---------------------------------------------------------------------------
# cross_serotype_scorecard
# ---------------------------------------------------------------------------
def test_scorecard_not_partition(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    *_, sc = _tables(stride_table, conservation_table)
    sc = sc.copy()
    sc.loc[sc["serotype"] == "DENV1", "n_signed"] = 99
    with pytest.raises(ConsistencyError, match="!= n_mechanisms"):
        validate_cross_serotype_scorecard(sc, S5Report())


def test_scorecard_reproducible_exceeds_loci(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    *_, sc = _tables(stride_table, conservation_table)
    sc = sc.copy()
    sc.loc[sc["serotype"] == "DENV1", "n_reproducible_residue"] = 99
    with pytest.raises(ConsistencyError, match="n_reproducible_residue"):
        validate_cross_serotype_scorecard(sc, S5Report())


def test_scorecard_shared_reproducible_exceeds_shared(
    stride_table: pd.DataFrame, conservation_table: pd.DataFrame
) -> None:
    *_, sc = _tables(stride_table, conservation_table)
    sc = sc.copy()
    sc.loc[sc["serotype"] == "DENV1", "n_shared_reproducible"] = 99
    with pytest.raises(ConsistencyError, match="n_shared_reproducible"):
        validate_cross_serotype_scorecard(sc, S5Report())


def test_validators_pass_on_empty() -> None:
    empty = pd.DataFrame()
    report = S5Report()
    validate_position_conservation(empty, report)
    validate_direction_concordance(empty, report)
    validate_domain_serotype_matrix(empty, report)
    validate_cross_serotype_scorecard(empty, report)
    assert report.all_passed
