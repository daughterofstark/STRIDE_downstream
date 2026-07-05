"""Prepare plotting-ready data for the publication figures F1–F8.

Each ``prepare_fN`` reads the relevant S2–S6 table(s), selects and orders the
columns the figure needs, rounds floats deterministically, and returns a
plotting-ready frame with exactly the declared schema. **No statistics are
computed** — every value is carried through from a prior stage. The plotting layer
(:mod:`stride_s7.plotting`) renders these frames; it never touches the source
tables.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    CATALYTIC_DOMAINS,
    F1_COLUMNS,
    F2_COLUMNS,
    F3_COLUMNS,
    F4_COLUMNS,
    F5_COLUMNS,
    F6_COLUMNS,
    F7_COLUMNS,
    F8_COLUMNS,
)
from ._common import round_floats, select_columns, sort_by


def prepare_f1(residue_landscape: pd.DataFrame) -> pd.DataFrame:
    """F1 — reproducibility landscape: ρ vs canonical residue, per serotype."""
    df = round_floats(residue_landscape, ["rho_residue"])
    df = sort_by(df, ["serotype", "chain", "canon_label"])
    return select_columns(df, F1_COLUMNS)


def prepare_f2(resolution_census: pd.DataFrame) -> pd.DataFrame:
    """F2 — achieved-resolution census: gated-scale loci counts per serotype.

    Uses the provisional-ρ\\* rows of the census (the single-threshold view), one
    row per ``(serotype, gated_scale_level)``.
    """
    from ._common import provisional_rows

    df = provisional_rows(resolution_census)
    df = df[df["n_loci"].astype(float) > 0] if not df.empty else df
    df = sort_by(df, ["serotype", "gated_scale_index"])
    return select_columns(df, F2_COLUMNS)


def prepare_f3(domain_serotype_matrix: pd.DataFrame) -> pd.DataFrame:
    """F3 — domain × serotype ρ heatmap (long form, ready to pivot)."""
    df = round_floats(domain_serotype_matrix, ["rho_domain"])
    df = sort_by(df, ["domain", "serotype"])
    return select_columns(df, F3_COLUMNS)


def prepare_f4(significance_screen: pd.DataFrame) -> pd.DataFrame:
    """F4 — signed-effect forest: β_signed ± CI for coherent (signed) mechanisms."""
    df = significance_screen
    df = df[df["is_signed"].astype(bool)] if not df.empty else df
    df = round_floats(df, ["beta_signed", "beta_ci_lower", "beta_ci_upper"])
    df = sort_by(df, ["serotype", "chain", "canon_label"])
    return select_columns(df, F4_COLUMNS)


def prepare_f5(position_conservation: pd.DataFrame) -> pd.DataFrame:
    """F5 — cross-serotype conservation map over shared positions."""
    df = round_floats(position_conservation, ["frac_reproducible"])
    df = sort_by(df, ["chain", "canon_label"])
    return select_columns(df, F5_COLUMNS)


def prepare_f6(variance_budget: pd.DataFrame) -> pd.DataFrame:
    """F6 — variance composition: τ²/σ̄² fractions stacked by domain."""
    df = round_floats(variance_budget, ["frac_tau2", "frac_sigma2"])
    df = sort_by(df, ["serotype", "chain", "domain"])
    return select_columns(df, F6_COLUMNS)


def prepare_f7(scale_curve: pd.DataFrame) -> pd.DataFrame:
    """F7 — ρ-vs-scale trajectories, restricted to the catalytic regions."""
    df = scale_curve
    if not df.empty:
        df = df[df["domain"].isin(CATALYTIC_DOMAINS)]
    df = round_floats(df, ["rho"])
    df = sort_by(df, ["serotype", "canon_label", "scale_index"])
    return select_columns(df, F7_COLUMNS)


def prepare_f8(domain_reproducibility: pd.DataFrame) -> pd.DataFrame:
    """F8 — coherence vs ρ scatter at the domain scale."""
    df = round_floats(domain_reproducibility, ["rho_domain", "coherence_domain"])
    df = sort_by(df, ["serotype", "chain", "domain"])
    return select_columns(df, F8_COLUMNS)


#: id → prepare function, for the orchestration layer.
FIGURE_BUILDERS = {
    "F1": prepare_f1,
    "F2": prepare_f2,
    "F3": prepare_f3,
    "F4": prepare_f4,
    "F5": prepare_f5,
    "F6": prepare_f6,
    "F7": prepare_f7,
    "F8": prepare_f8,
}


def build_all_figures(frames: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Prepare every figure's plotting-ready data from the loaded input frames."""
    from ..models.schema import (
        IN_DOMAIN_REPRODUCIBILITY,
        IN_DOMAIN_SEROTYPE_MATRIX,
        IN_POSITION_CONSERVATION,
        IN_RESIDUE_LANDSCAPE,
        IN_RESOLUTION_CENSUS,
        IN_SCALE_CURVE,
        IN_SIGNIFICANCE_SCREEN,
        IN_VARIANCE_BUDGET,
    )

    return {
        "F1": prepare_f1(frames[IN_RESIDUE_LANDSCAPE]),
        "F2": prepare_f2(frames[IN_RESOLUTION_CENSUS]),
        "F3": prepare_f3(frames[IN_DOMAIN_SEROTYPE_MATRIX]),
        "F4": prepare_f4(frames[IN_SIGNIFICANCE_SCREEN]),
        "F5": prepare_f5(frames[IN_POSITION_CONSERVATION]),
        "F6": prepare_f6(frames[IN_VARIANCE_BUDGET]),
        "F7": prepare_f7(frames[IN_SCALE_CURVE]),
        "F8": prepare_f8(frames[IN_DOMAIN_REPRODUCIBILITY]),
    }
