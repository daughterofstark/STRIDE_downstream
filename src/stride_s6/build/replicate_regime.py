"""Build the per-serotype replicate-regime ledger.

One row per serotype summarising the replicate structure from the S1A
``replicate_inventory``: the replicate count *K*, how many shared positions are
observed and how many are present in every run (completeness), whether
residue-scale claims are licensed (K ≥ 5, design §0.1), and — from the optional
S0 ``replicate_table`` — whether genuine per-run effects are available and in how
many runs.

Reads inputs, never mutates them; deterministic (serotype-sorted output).
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    MIN_K_FOR_RESIDUE_LICENSE,
    MIN_REPLICATES_FOR_CONCORDANCE,
    REPLICATE_REGIME_COLUMNS,
    ROUND_DECIMALS,
)
from ._frames import runs_with_effects


def build_replicate_regime(
    replicate_inventory: pd.DataFrame,
    effects: pd.DataFrame,
) -> pd.DataFrame:
    """Return the per-serotype replicate-regime table.

    Parameters
    ----------
    replicate_inventory
        The S1A replicate inventory (one row per ``(serotype, canon_label)``).
    effects
        The tidied per-run effects from :func:`stride_s6.build._frames.per_run_effects`
        (empty when the replicate table is absent).

    Returns
    -------
    DataFrame
        Columns :data:`REPLICATE_REGIME_COLUMNS`, one row per serotype.
    """
    if replicate_inventory.empty:
        return pd.DataFrame(columns=list(REPLICATE_REGIME_COLUMNS))

    inv = replicate_inventory.copy()
    inv["serotype"] = inv["serotype"].astype(str)
    inv["n_replicates"] = pd.to_numeric(inv["n_replicates"], errors="coerce")
    inv["in_all_replicates"] = inv["in_all_replicates"].astype(bool)

    rows = []
    for serotype in sorted(inv["serotype"].unique()):
        sub = inv[inv["serotype"] == serotype]
        n_positions = int(len(sub))
        # K: the replicate count for the serotype (max over positions; positions
        # normally share it — S0 enforces equal replicate counts by default).
        k_series = sub["n_replicates"].dropna()
        n_replicates = int(k_series.max()) if not k_series.empty else 0
        n_in_all = int(sub["in_all_replicates"].sum())
        frac_complete = (
            round(n_in_all / n_positions, ROUND_DECIMALS)
            if n_positions
            else 0.0
        )
        runs = runs_with_effects(effects, serotype)
        n_runs_effects = len(runs)
        available = n_runs_effects >= MIN_REPLICATES_FOR_CONCORDANCE
        rows.append(
            {
                "serotype": serotype,
                "n_replicates": n_replicates,
                "n_positions": n_positions,
                "n_positions_in_all_replicates": n_in_all,
                "frac_complete": frac_complete,
                "residue_claims_licensed": bool(
                    n_replicates >= MIN_K_FOR_RESIDUE_LICENSE
                ),
                "per_replicate_effects_available": bool(available),
                "n_replicates_with_effects": int(n_runs_effects),
            }
        )

    out = pd.DataFrame(rows, columns=list(REPLICATE_REGIME_COLUMNS))
    return out.reset_index(drop=True)
