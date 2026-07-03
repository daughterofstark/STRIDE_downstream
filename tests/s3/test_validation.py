"""Tests for the S3 structural validation checks."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s3.build import (
    build_chain_contrast,
    build_monotonicity_audit,
    build_resolution_gap,
    build_scale_curve,
)
from stride_s3.models import S3Report
from stride_s3.models.errors import ConsistencyError
from stride_s3.validation import (
    validate_chain_contrast_totals,
    validate_gap_consistency,
    validate_monotonicity_audit_consistency,
    validate_scale_curve_completeness,
    validate_unique_keys,
)


def _all_tables(
    stride_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sc = build_scale_curve(stride_table)
    rg = build_resolution_gap(stride_table)
    ma = build_monotonicity_audit(stride_table)
    cc = build_chain_contrast(stride_table)
    return sc, rg, ma, cc


def test_all_checks_pass_on_valid_tables(stride_table: pd.DataFrame) -> None:
    sc, rg, ma, cc = _all_tables(stride_table)
    report = S3Report()
    validate_unique_keys(sc, rg, ma, cc, report)
    validate_scale_curve_completeness(sc, rg, report)
    validate_gap_consistency(rg, report)
    validate_monotonicity_audit_consistency(ma, report)
    validate_chain_contrast_totals(cc, report)
    assert report.all_passed
    assert len(report.checks) == 5


def test_duplicate_key_raises(stride_table: pd.DataFrame) -> None:
    sc, rg, ma, cc = _all_tables(stride_table)
    dup = pd.concat([rg, rg.iloc[[0]]], ignore_index=True)
    with pytest.raises(ConsistencyError, match="resolution_gap key"):
        validate_unique_keys(sc, dup, ma, cc, S3Report())


def test_scale_curve_incomplete_raises(stride_table: pd.DataFrame) -> None:
    sc, rg, _ma, _cc = _all_tables(stride_table)
    # drop one scale row for a locus → incomplete curve
    bad = sc.drop(sc[(sc.canon_label == "NS3:51") & (sc.scale_index == 6)].index)
    with pytest.raises(ConsistencyError, match="scale rows|scale indices"):
        validate_scale_curve_completeness(bad, rg, S3Report())


def test_gap_arithmetic_violation_raises(stride_table: pd.DataFrame) -> None:
    _sc, rg, _ma, _cc = _all_tables(stride_table)
    bad = rg.copy()
    bad.loc[bad.index[0], "delta_rho_domain_residue"] = 999.0
    with pytest.raises(ConsistencyError, match="delta"):
        validate_gap_consistency(bad, S3Report())


def test_monotonicity_flag_mismatch_raises(stride_table: pd.DataFrame) -> None:
    _sc, _rg, ma, _cc = _all_tables(stride_table)
    bad = ma.copy()
    # claim monotone while leaving a positive violation count
    idx = bad[~bad.is_monotone].index[0]
    bad.loc[idx, "is_monotone"] = True
    with pytest.raises(ConsistencyError, match="disagrees with"):
        validate_monotonicity_audit_consistency(bad, S3Report())


def test_chain_contrast_partition_violation_raises(
    stride_table: pd.DataFrame,
) -> None:
    _sc, _rg, _ma, cc = _all_tables(stride_table)
    bad = cc.copy()
    bad.loc[bad.index[0], "n_increase"] = 999
    with pytest.raises(ConsistencyError, match="n_mechanisms"):
        validate_chain_contrast_totals(bad, S3Report())


def test_scale_curve_gap_mismatch_raises(stride_table: pd.DataFrame) -> None:
    sc, rg, _ma, _cc = _all_tables(stride_table)
    bad_gap = rg.copy()
    bad_gap.loc[bad_gap.index[0], "rho_residue"] = 0.123456
    with pytest.raises(ConsistencyError, match="rho_residue"):
        validate_scale_curve_completeness(sc, bad_gap, S3Report())


# ---------------------------------------------------------------------------
# empty-input tolerance
# ---------------------------------------------------------------------------
def test_checks_tolerate_empty_tables() -> None:
    empty = pd.DataFrame()
    report = S3Report()
    validate_scale_curve_completeness(empty, empty, report)
    validate_gap_consistency(empty, report)
    validate_monotonicity_audit_consistency(empty, report)
    validate_chain_contrast_totals(empty, report)
    assert report.all_passed
