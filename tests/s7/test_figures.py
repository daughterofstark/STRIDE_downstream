"""Tests for the figure-data builders (prepare_f1..prepare_f8)."""
from __future__ import annotations

import pandas as pd

from stride_s7.build import build_all_figures
from stride_s7.build.figures import (
    prepare_f1,
    prepare_f2,
    prepare_f3,
    prepare_f4,
    prepare_f5,
    prepare_f6,
    prepare_f7,
    prepare_f8,
)
from stride_s7.models.schema import (
    CATALYTIC_DOMAINS,
    FIGURE_DATA_COLUMNS,
    FIGURE_IDS,
    IN_DOMAIN_REPRODUCIBILITY,
    IN_DOMAIN_SEROTYPE_MATRIX,
    IN_POSITION_CONSERVATION,
    IN_RESIDUE_LANDSCAPE,
    IN_RESOLUTION_CENSUS,
    IN_SCALE_CURVE,
    IN_SIGNIFICANCE_SCREEN,
    IN_VARIANCE_BUDGET,
)


def test_build_all_figures_has_every_id(inputs: dict[str, pd.DataFrame]) -> None:
    figs = build_all_figures(inputs)
    assert set(figs) == set(FIGURE_IDS)


def test_every_figure_has_declared_columns(inputs: dict[str, pd.DataFrame]) -> None:
    figs = build_all_figures(inputs)
    for fid in FIGURE_IDS:
        assert list(figs[fid].columns) == list(FIGURE_DATA_COLUMNS[fid])


def test_figures_are_deterministic(inputs: dict[str, pd.DataFrame]) -> None:
    a = build_all_figures(inputs)
    b = build_all_figures(inputs)
    for fid in FIGURE_IDS:
        pd.testing.assert_frame_equal(a[fid], b[fid])


def test_f1_sorted_and_nonempty(inputs: dict[str, pd.DataFrame]) -> None:
    df = prepare_f1(inputs[IN_RESIDUE_LANDSCAPE])
    assert not df.empty
    assert df["serotype"].is_monotonic_increasing


def test_f2_uses_provisional_rows_only(inputs: dict[str, pd.DataFrame]) -> None:
    df = prepare_f2(inputs[IN_RESOLUTION_CENSUS])
    # census has provisional and non-provisional rows; only provisional survive
    assert not df.empty
    assert (df["n_loci"] > 0).all()


def test_f3_is_domain_serotype_long(inputs: dict[str, pd.DataFrame]) -> None:
    df = prepare_f3(inputs[IN_DOMAIN_SEROTYPE_MATRIX])
    assert {"domain", "serotype", "rho_domain"} <= set(df.columns)
    assert len(df) == len(inputs[IN_DOMAIN_SEROTYPE_MATRIX])


def test_f4_keeps_only_signed(inputs: dict[str, pd.DataFrame]) -> None:
    df = prepare_f4(inputs[IN_SIGNIFICANCE_SCREEN])
    src = inputs[IN_SIGNIFICANCE_SCREEN]
    assert len(df) == int(src["is_signed"].astype(bool).sum())
    assert df["beta_signed"].notna().all()


def test_f5_conservation_columns(inputs: dict[str, pd.DataFrame]) -> None:
    df = prepare_f5(inputs[IN_POSITION_CONSERVATION])
    assert "conservation_class" in df.columns
    assert (df["frac_reproducible"] <= 1.0).all()


def test_f6_fractions_present(inputs: dict[str, pd.DataFrame]) -> None:
    df = prepare_f6(inputs[IN_VARIANCE_BUDGET])
    assert {"frac_tau2", "frac_sigma2"} <= set(df.columns)


def test_f7_only_catalytic_regions(inputs: dict[str, pd.DataFrame]) -> None:
    df = prepare_f7(inputs[IN_SCALE_CURVE])
    assert not df.empty
    assert set(df["domain"].unique()) <= set(CATALYTIC_DOMAINS)


def test_f8_domain_scale_scatter(inputs: dict[str, pd.DataFrame]) -> None:
    df = prepare_f8(inputs[IN_DOMAIN_REPRODUCIBILITY])
    assert {"rho_domain", "coherence_domain", "is_coherent"} <= set(df.columns)


def test_figures_do_not_mutate_inputs(inputs: dict[str, pd.DataFrame]) -> None:
    before = {k: v.copy() for k, v in inputs.items()}
    build_all_figures(inputs)
    for k in before:
        pd.testing.assert_frame_equal(inputs[k], before[k])


def test_figures_on_empty_inputs_have_columns() -> None:
    from tests.s7.fixtures import make_empty_inputs

    figs = build_all_figures(make_empty_inputs())
    for fid in FIGURE_IDS:
        assert list(figs[fid].columns) == list(FIGURE_DATA_COLUMNS[fid])
        assert figs[fid].empty
