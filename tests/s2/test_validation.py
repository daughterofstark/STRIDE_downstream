"""Tests for the S2 structural validation checks."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s2.build import (
    build_domain_reproducibility,
    build_residue_landscape,
    build_resolution_census,
    build_serotype_summary,
    build_signed_screen,
)
from stride_s2.models import S2Report
from stride_s2.models.errors import ConsistencyError
from stride_s2.validation import (
    validate_census_totals,
    validate_regating_monotonicity,
    validate_serotype_summary_consistency,
    validate_tiers,
    validate_unique_keys,
)

_BAND = (0.5, 0.7, 0.9)
_DECIMALS = 4


def _all_tables(
    stride_table: pd.DataFrame, ra: pd.DataFrame, da: pd.DataFrame
) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    census = build_resolution_census(stride_table, _BAND, _DECIMALS)
    rl = build_residue_landscape(stride_table, ra)
    dr = build_domain_reproducibility(stride_table, da)
    ss = build_signed_screen(stride_table, _BAND, _DECIMALS)
    sm = build_serotype_summary(stride_table, census, rl, ss, _BAND, _DECIMALS)
    return census, rl, dr, ss, sm


def test_all_checks_pass_on_valid_tables(
    stride_table: pd.DataFrame,
    residue_annotation: pd.DataFrame,
    domain_annotation: pd.DataFrame,
) -> None:
    census, rl, dr, ss, sm = _all_tables(
        stride_table, residue_annotation, domain_annotation
    )
    report = S2Report()
    validate_unique_keys(census, rl, dr, ss, sm, report)
    validate_census_totals(census, rl, report)
    validate_regating_monotonicity(census, report)
    validate_serotype_summary_consistency(sm, census, report)
    validate_tiers(rl, dr, report)
    assert report.all_passed
    assert len(report.checks) == 5


def test_duplicate_key_raises(
    stride_table: pd.DataFrame,
    residue_annotation: pd.DataFrame,
    domain_annotation: pd.DataFrame,
) -> None:
    census, rl, dr, ss, sm = _all_tables(
        stride_table, residue_annotation, domain_annotation
    )
    dup_rl = pd.concat([rl, rl.iloc[[0]]], ignore_index=True)
    with pytest.raises(ConsistencyError, match="residue_landscape key"):
        validate_unique_keys(census, dup_rl, dr, ss, sm, S2Report())


def test_census_total_mismatch_raises(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    census = build_resolution_census(stride_table, _BAND, _DECIMALS)
    rl = build_residue_landscape(stride_table, residue_annotation)
    # corrupt a census count so the total no longer matches the locus count
    bad = census.copy()
    bad.loc[bad.index[0], "n_loci"] = 999
    with pytest.raises(ConsistencyError, match="sums to"):
        validate_census_totals(bad, rl, S2Report())


def test_regating_monotonicity_violation_raises(
    stride_table: pd.DataFrame,
) -> None:
    census = build_resolution_census(stride_table, _BAND, _DECIMALS)
    # inject a residue-gated count that increases with ρ*: set ρ*=0.5 residue
    # count to 0 and ρ*=0.9 residue count high
    bad = census.copy()
    residue = bad["gated_scale_level"] == "residue"
    bad.loc[residue & (bad.rho_star == 0.5), "n_loci"] = 0
    bad.loc[residue & (bad.rho_star == 0.9), "n_loci"] = 5
    with pytest.raises(ConsistencyError, match="monotone"):
        validate_regating_monotonicity(bad, S2Report())


def test_serotype_summary_partition_violation_raises(
    stride_table: pd.DataFrame, residue_annotation: pd.DataFrame
) -> None:
    census, rl, _dr, ss, sm = _all_tables(
        stride_table, residue_annotation, make_domain_annotation_ok()
    )
    bad = sm.copy()
    bad.loc[bad.index[0], "n_signed"] = 999
    with pytest.raises(ConsistencyError, match="n_signed"):
        validate_serotype_summary_consistency(bad, census, S2Report())


def test_tier_violation_raises(
    stride_table: pd.DataFrame,
    residue_annotation: pd.DataFrame,
    domain_annotation: pd.DataFrame,
) -> None:
    _c, rl, dr, _s, _m = _all_tables(
        stride_table, residue_annotation, domain_annotation
    )
    bad_rl = rl.copy()
    bad_rl.loc[bad_rl.index[0], "tier"] = "licensed"
    with pytest.raises(ConsistencyError, match="non-exploratory"):
        validate_tiers(bad_rl, dr, S2Report())


def make_domain_annotation_ok() -> pd.DataFrame:
    from tests.s2.fixtures import make_domain_annotation

    return make_domain_annotation()


# ---------------------------------------------------------------------------
# empty-input tolerance
# ---------------------------------------------------------------------------
def test_checks_tolerate_empty_tables() -> None:
    empty = pd.DataFrame()
    report = S2Report()
    validate_census_totals(empty, empty, report)
    validate_regating_monotonicity(empty, report)
    validate_serotype_summary_consistency(empty, empty, report)
    validate_tiers(empty, empty, report)
    assert report.all_passed
