"""Tests for the S4 structural validation checks."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s4.build import (
    build_domain_effect_summary,
    build_residue_variance,
    build_significance_screen,
    build_variance_budget,
)
from stride_s4.models import S4Report
from stride_s4.models.errors import ConsistencyError
from stride_s4.validation import (
    validate_domain_effect_totals,
    validate_significance_screen,
    validate_unique_keys,
    validate_variance_fractions,
)


def _all_tables(
    stride_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    vb = build_variance_budget(stride_table)
    rv = build_residue_variance(stride_table)
    ss = build_significance_screen(stride_table)
    des = build_domain_effect_summary(ss)
    return vb, rv, ss, des


def test_all_checks_pass_on_valid_tables(stride_table: pd.DataFrame) -> None:
    vb, rv, ss, des = _all_tables(stride_table)
    report = S4Report()
    validate_unique_keys(vb, rv, ss, des, report)
    validate_variance_fractions(vb, rv, report)
    validate_significance_screen(ss, report)
    validate_domain_effect_totals(des, report)
    assert report.all_passed
    assert len(report.checks) == 4


def test_duplicate_key_raises(stride_table: pd.DataFrame) -> None:
    vb, rv, ss, des = _all_tables(stride_table)
    dup = pd.concat([vb, vb.iloc[[0]]], ignore_index=True)
    with pytest.raises(ConsistencyError, match="variance_budget key"):
        validate_unique_keys(dup, rv, ss, des, S4Report())


def test_variance_fraction_out_of_range_raises(
    stride_table: pd.DataFrame,
) -> None:
    vb, rv, _ss, _des = _all_tables(stride_table)
    bad = vb.copy()
    bad.loc[bad.index[0], "frac_tau2"] = 1.5
    with pytest.raises(ConsistencyError, match="outside \\[0, 1\\]"):
        validate_variance_fractions(bad, rv, S4Report())


def test_variance_fraction_sum_violation_raises(
    stride_table: pd.DataFrame,
) -> None:
    vb, rv, _ss, _des = _all_tables(stride_table)
    bad = vb.copy()
    bad.loc[bad.index[0], "frac_sigma2"] = 0.9  # no longer sums to 1 with frac_tau2
    with pytest.raises(ConsistencyError, match="!= 1|outside"):
        validate_variance_fractions(bad, rv, S4Report())


def test_variance_regime_mismatch_raises(stride_table: pd.DataFrame) -> None:
    vb, rv, _ss, _des = _all_tables(stride_table)
    bad = vb.copy()
    # flip a replicate_dominated label to sampling_dominated
    idx = bad[bad.variance_regime == "replicate_dominated"].index[0]
    bad.loc[idx, "variance_regime"] = "sampling_dominated"
    with pytest.raises(ConsistencyError, match="regime"):
        validate_variance_fractions(bad, rv, S4Report())


def test_significance_signed_mismatch_raises(stride_table: pd.DataFrame) -> None:
    _vb, _rv, ss, _des = _all_tables(stride_table)
    bad = ss.copy()
    # mark a mixed row as signed without giving it a direction
    idx = bad[~bad.is_signed].index[0]
    bad.loc[idx, "is_signed"] = True
    with pytest.raises(ConsistencyError, match="is_signed"):
        validate_significance_screen(bad, S4Report())


def test_significance_unsigned_with_ci_raises(stride_table: pd.DataFrame) -> None:
    _vb, _rv, ss, _des = _all_tables(stride_table)
    bad = ss.copy()
    idx = bad[~bad.is_signed].index[0]
    bad.loc[idx, "ci_excludes_zero"] = True
    with pytest.raises(ConsistencyError, match="unsigned row has ci_excludes_zero"):
        validate_significance_screen(bad, S4Report())


def test_significance_raw_mismatch_raises(stride_table: pd.DataFrame) -> None:
    _vb, _rv, ss, _des = _all_tables(stride_table)
    bad = ss.copy()
    idx = bad[bad.ci_excludes_zero].index[0]
    bad.loc[idx, "significant_raw"] = False
    with pytest.raises(ConsistencyError, match="significant_raw"):
        validate_significance_screen(bad, S4Report())


def test_domain_effect_partition_violation_raises(
    stride_table: pd.DataFrame,
) -> None:
    _vb, _rv, _ss, des = _all_tables(stride_table)
    bad = des.copy()
    bad.loc[bad.index[0], "n_mixed"] = 99
    with pytest.raises(ConsistencyError, match="n_mechanisms"):
        validate_domain_effect_totals(bad, S4Report())


def test_domain_effect_ci_exceeds_signed_raises(
    stride_table: pd.DataFrame,
) -> None:
    _vb, _rv, _ss, des = _all_tables(stride_table)
    bad = des.copy()
    bad.loc[bad.index[0], "n_ci_excludes_zero"] = 999
    with pytest.raises(ConsistencyError, match="n_ci_excludes_zero"):
        validate_domain_effect_totals(bad, S4Report())


# ---------------------------------------------------------------------------
# empty-input tolerance
# ---------------------------------------------------------------------------
def test_checks_tolerate_empty_tables() -> None:
    empty = pd.DataFrame()
    report = S4Report()
    validate_variance_fractions(empty, empty, report)
    validate_significance_screen(empty, report)
    validate_domain_effect_totals(empty, report)
    assert report.all_passed
