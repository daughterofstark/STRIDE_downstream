"""Build the per-serotype rank-concordance table (Tier B, exploratory).

Implements design §3.1's per-run rank concordance: for each serotype that carries
usable per-run effects, the Kendall's *W* and the mean pairwise Spearman
correlation of the per-run effect rankings across its runs, over the positions
observed in *every* run (complete cases). One row per serotype with effects; the
table is empty when no per-run effects are available (the analysis is recorded as
blocked in the ledger). A serotype with too few complete positions is reported
with an ``insufficient`` class and ``nan`` coefficients rather than a fabricated
number.

Reads inputs, never mutates them; deterministic (serotype-sorted output).
"""
from __future__ import annotations

import math

import pandas as pd

from ..models.schema import (
    CONCORDANCE_INSUFFICIENT,
    CONCORDANCE_MODERATE,
    CONCORDANCE_MODERATE_W,
    CONCORDANCE_STRONG,
    CONCORDANCE_STRONG_W,
    CONCORDANCE_WEAK,
    MIN_POSITIONS_FOR_CONCORDANCE,
    MIN_REPLICATES_FOR_CONCORDANCE,
    REPLICATE_CONCORDANCE_COLUMNS,
    ROUND_DECIMALS,
    TIER_EXPLORATORY,
)
from ._concordance import kendalls_w, mean_pairwise_spearman
from ._frames import serotype_effect_matrix


def _concordance_class(w: float, n_positions: int, n_runs: int) -> str:
    """Label the concordance from *W* and the sample support."""
    if (
        n_runs < MIN_REPLICATES_FOR_CONCORDANCE
        or n_positions < MIN_POSITIONS_FOR_CONCORDANCE
        or math.isnan(w)
    ):
        return CONCORDANCE_INSUFFICIENT
    if w >= CONCORDANCE_STRONG_W:
        return CONCORDANCE_STRONG
    if w >= CONCORDANCE_MODERATE_W:
        return CONCORDANCE_MODERATE
    return CONCORDANCE_WEAK


def _round_or_nan(value: float) -> float:
    return value if math.isnan(value) else round(value, ROUND_DECIMALS)


def build_replicate_concordance(effects: pd.DataFrame) -> pd.DataFrame:
    """Return the per-serotype rank-concordance table.

    Parameters
    ----------
    effects
        Tidied per-run effects (empty when the replicate table is absent).

    Returns
    -------
    DataFrame
        Columns :data:`REPLICATE_CONCORDANCE_COLUMNS`, one row per serotype that
        carries per-run effects (empty when none do).
    """
    if effects.empty:
        return pd.DataFrame(columns=list(REPLICATE_CONCORDANCE_COLUMNS))

    rows = []
    for serotype in sorted(effects["serotype"].astype(str).unique()):
        positions, runs, matrix = serotype_effect_matrix(effects, serotype)
        n_runs = len(runs)
        n_positions = len(positions)
        if (
            n_runs >= MIN_REPLICATES_FOR_CONCORDANCE
            and n_positions >= MIN_POSITIONS_FOR_CONCORDANCE
        ):
            w = kendalls_w(matrix)
            rho = mean_pairwise_spearman(matrix)
        else:
            w = math.nan
            rho = math.nan
        rows.append(
            {
                "serotype": serotype,
                "n_replicates_with_effects": int(n_runs),
                "n_positions_complete": int(n_positions),
                "kendalls_w": _round_or_nan(w),
                "mean_pairwise_spearman": _round_or_nan(rho),
                "concordance_class": _concordance_class(w, n_positions, n_runs),
                "tier": TIER_EXPLORATORY,
            }
        )

    out = pd.DataFrame(rows, columns=list(REPLICATE_CONCORDANCE_COLUMNS))
    return out.reset_index(drop=True)
