"""Assemble the manuscript tables T1–T5 from the S2–S6 outputs.

Each ``build_tN`` reads the relevant prior-stage table(s), selects / joins / orders
the columns the manuscript table needs, and returns a frame with exactly the
declared schema. **No statistics are recomputed** — the builders only assemble and
format values that S2–S6 already produced.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    CATALYTIC_DOMAINS,
    T1_COLUMNS,
    T2_COLUMNS,
    T3_COLUMNS,
    T4_COLUMNS,
    T4_TOP_N,
    T5_COLUMNS,
)
from ._common import provisional_rows, round_floats, select_columns, sort_by

_T1_FLOATS = [
    "frac_mixed",
    "rho_residue_median",
    "rho_residue_q1",
    "rho_residue_q3",
    "rho_residue_min",
    "rho_residue_max",
    "rho_star",
]
_T2_FLOATS = [
    "rho_domain",
    "beta_domain",
    "coherence_domain",
    "beta_weighted_mean",
    "beta_weighted_se",
]
_T3_FLOATS = ["rho_domain", "beta_domain"]
_T4_FLOATS = ["frac_reproducible", "rho_residue_median"]
_T5_FLOATS = ["tau2", "sigma2_bar", "frac_tau2", "frac_sigma2", "tau2_sigma2_ratio"]


def build_t1(serotype_summary: pd.DataFrame) -> pd.DataFrame:
    """T1 — per-serotype summary (loci, census, ρ median/IQR, %mixed, signed-sig)."""
    df = provisional_rows(serotype_summary)
    df = round_floats(df, _T1_FLOATS)
    df = sort_by(df, ["serotype"])
    return select_columns(df, T1_COLUMNS)


def build_t2(
    domain_reproducibility: pd.DataFrame,
    domain_effect_summary: pd.DataFrame,
) -> pd.DataFrame:
    """T2 — domain-level ρ and signed effect per serotype.

    Joins the S2 domain reproducibility (ρ, unsigned β, coherence) with the S4
    β_se-weighted signed-effect summary on ``(serotype, chain, domain)``.
    """
    left = domain_reproducibility.loc[
        :,
        [
            "serotype",
            "chain",
            "domain",
            "rho_domain",
            "beta_domain",
            "coherence_domain",
            "is_coherent",
            "tier",
        ],
    ]
    right = domain_effect_summary.loc[
        :,
        [
            "serotype",
            "chain",
            "domain",
            "beta_weighted_mean",
            "beta_weighted_se",
            "n_signed",
            "n_significant_fdr",
        ],
    ]
    merged = left.merge(right, on=["serotype", "chain", "domain"], how="left")
    merged = round_floats(merged, _T2_FLOATS)
    merged = sort_by(merged, ["serotype", "chain", "domain"])
    return select_columns(merged, T2_COLUMNS)


def build_t3(domain_serotype_matrix: pd.DataFrame) -> pd.DataFrame:
    """T3 — Catalytic Triad / Oxyanion Loop cross-serotype behaviour.

    The domain × serotype matrix restricted to the catalytic machinery.
    """
    df = domain_serotype_matrix
    if not df.empty:
        mask = df["is_catalytic_domain"].astype(bool) | df["domain"].isin(
            CATALYTIC_DOMAINS
        )
        df = df[mask]
    df = round_floats(df, _T3_FLOATS)
    df = sort_by(df, ["domain", "serotype"])
    return select_columns(df, T3_COLUMNS)


def build_t4(position_conservation: pd.DataFrame) -> pd.DataFrame:
    """T4 — top reproducible signed positions shared across serotypes.

    Ranks shared positions by the number of serotypes in which they are
    signed-and-reproducible, then by reproducible fraction and median ρ, and keeps
    the top ``T4_TOP_N``. Ranking is a sort over existing columns, not a new metric.
    """
    df = round_floats(position_conservation, _T4_FLOATS)
    if not df.empty:
        df = df.sort_values(
            [
                "n_serotypes_signed_reproducible",
                "n_serotypes_reproducible",
                "frac_reproducible",
                "rho_residue_median",
                "canon_label",
            ],
            ascending=[False, False, False, False, True],
            kind="stable",
        ).head(T4_TOP_N)
    return select_columns(df, T4_COLUMNS)


def build_t5(variance_budget: pd.DataFrame) -> pd.DataFrame:
    """T5 — variance-component budget per domain (τ²/σ̄² split + regime)."""
    df = round_floats(variance_budget, _T5_FLOATS)
    df = sort_by(df, ["serotype", "chain", "domain"])
    return select_columns(df, T5_COLUMNS)


#: id → (builder, tuple of source filenames), for the orchestration layer.
def build_all_tables(frames: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Build every manuscript table from the loaded input frames."""
    from ..models.schema import (
        IN_DOMAIN_EFFECT_SUMMARY,
        IN_DOMAIN_REPRODUCIBILITY,
        IN_DOMAIN_SEROTYPE_MATRIX,
        IN_POSITION_CONSERVATION,
        IN_SEROTYPE_SUMMARY,
        IN_VARIANCE_BUDGET,
    )

    return {
        "T1": build_t1(frames[IN_SEROTYPE_SUMMARY]),
        "T2": build_t2(
            frames[IN_DOMAIN_REPRODUCIBILITY],
            frames[IN_DOMAIN_EFFECT_SUMMARY],
        ),
        "T3": build_t3(frames[IN_DOMAIN_SEROTYPE_MATRIX]),
        "T4": build_t4(frames[IN_POSITION_CONSERVATION]),
        "T5": build_t5(frames[IN_VARIANCE_BUDGET]),
    }
