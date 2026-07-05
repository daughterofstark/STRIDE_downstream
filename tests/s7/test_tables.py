"""Tests for the manuscript-table builders (build_t1..build_t5)."""
from __future__ import annotations

import pandas as pd

from stride_s7.build import build_all_tables
from stride_s7.build.tables import (
    build_t1,
    build_t2,
    build_t3,
    build_t4,
    build_t5,
)
from stride_s7.models.schema import (
    CATALYTIC_DOMAINS,
    IN_DOMAIN_EFFECT_SUMMARY,
    IN_DOMAIN_REPRODUCIBILITY,
    IN_DOMAIN_SEROTYPE_MATRIX,
    IN_POSITION_CONSERVATION,
    IN_SEROTYPE_SUMMARY,
    IN_VARIANCE_BUDGET,
    T4_TOP_N,
    TABLE_COLUMNS,
    TABLE_IDS,
)


def test_build_all_tables_has_every_id(inputs: dict[str, pd.DataFrame]) -> None:
    tables = build_all_tables(inputs)
    assert set(tables) == set(TABLE_IDS)


def test_every_table_has_declared_columns(inputs: dict[str, pd.DataFrame]) -> None:
    tables = build_all_tables(inputs)
    for tid in TABLE_IDS:
        assert list(tables[tid].columns) == list(TABLE_COLUMNS[tid])


def test_tables_are_deterministic(inputs: dict[str, pd.DataFrame]) -> None:
    a = build_all_tables(inputs)
    b = build_all_tables(inputs)
    for tid in TABLE_IDS:
        pd.testing.assert_frame_equal(a[tid], b[tid])


def test_t1_one_row_per_serotype(inputs: dict[str, pd.DataFrame]) -> None:
    df = build_t1(inputs[IN_SEROTYPE_SUMMARY])
    assert len(df) == inputs[IN_SEROTYPE_SUMMARY]["serotype"].nunique()
    # only the provisional rho* rows survive
    assert (df["rho_star"] == 0.5).all()


def test_t2_join_brings_in_signed_effect(inputs: dict[str, pd.DataFrame]) -> None:
    df = build_t2(
        inputs[IN_DOMAIN_REPRODUCIBILITY], inputs[IN_DOMAIN_EFFECT_SUMMARY]
    )
    assert "beta_weighted_mean" in df.columns
    assert "rho_domain" in df.columns
    assert len(df) == len(inputs[IN_DOMAIN_REPRODUCIBILITY])
    assert df["beta_weighted_mean"].notna().all()


def test_t3_only_catalytic(inputs: dict[str, pd.DataFrame]) -> None:
    df = build_t3(inputs[IN_DOMAIN_SEROTYPE_MATRIX])
    assert not df.empty
    assert set(df["domain"].unique()) <= set(CATALYTIC_DOMAINS)


def test_t4_ranked_and_capped(inputs: dict[str, pd.DataFrame]) -> None:
    df = build_t4(inputs[IN_POSITION_CONSERVATION])
    assert len(df) <= T4_TOP_N
    # ranked by n_serotypes_signed_reproducible descending
    vals = df["n_serotypes_signed_reproducible"].tolist()
    assert vals == sorted(vals, reverse=True)


def test_t5_variance_budget_passthrough(inputs: dict[str, pd.DataFrame]) -> None:
    df = build_t5(inputs[IN_VARIANCE_BUDGET])
    assert len(df) == len(inputs[IN_VARIANCE_BUDGET])
    assert "variance_regime" in df.columns


def test_tables_do_not_mutate_inputs(inputs: dict[str, pd.DataFrame]) -> None:
    before = {k: v.copy() for k, v in inputs.items()}
    build_all_tables(inputs)
    for k in before:
        pd.testing.assert_frame_equal(inputs[k], before[k])


def test_tables_on_empty_inputs_have_columns() -> None:
    from tests.s7.fixtures import make_empty_inputs

    tables = build_all_tables(make_empty_inputs())
    for tid in TABLE_IDS:
        assert list(tables[tid].columns) == list(TABLE_COLUMNS[tid])
        assert tables[tid].empty
