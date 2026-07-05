"""Build the descriptive per-run effect-spread table (Tier B, exploratory).

One row per observed ``(serotype, canon_label)`` summarising the spread of that
position's per-run θ across the runs it appears in: the count, mean, sample sd,
min / max / range, mean absolute effect, and the largest pairwise absolute
difference between any two runs. This is a *descriptive* summary of the Level-1
observations that are actually present — it makes no inferential claim and, being
residue-scale, is labelled ``exploratory``. When no per-run effects are available
the table is empty (the analysis is recorded as blocked in the ledger).

Reads inputs, never mutates them; deterministic (``(serotype, canon_label)``-sorted).
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    REPLICATE_EFFECT_SPREAD_COLUMNS,
    ROUND_DECIMALS,
    TIER_EXPLORATORY,
)
from ._frames import chain_from_canon


def build_replicate_effect_spread(
    effects: pd.DataFrame,
    replicate_inventory: pd.DataFrame,
    replicate_table: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Return the per-position per-run effect-spread table.

    Parameters
    ----------
    effects
        Tidied per-run effects (empty when the replicate table is absent).
    replicate_inventory
        The S1A inventory, used for the ``in_all_replicates`` flag per position.
    replicate_table
        The raw S0 replicate table, used only to recover the optional
        ``domain_label`` context column when present.
    """
    if effects.empty:
        return pd.DataFrame(columns=list(REPLICATE_EFFECT_SPREAD_COLUMNS))

    in_all = _in_all_lookup(replicate_inventory)
    domain = _domain_lookup(replicate_table)

    rows = []
    grouped = effects.groupby(["serotype", "canon_label"], sort=True)
    for (serotype, canon), sub in grouped:
        thetas = [float(v) for v in sub["r"].tolist()]
        n_obs = len(thetas)
        theta_mean = sum(thetas) / n_obs
        theta_min = min(thetas)
        theta_max = max(thetas)
        theta_range = theta_max - theta_min
        abs_mean = sum(abs(t) for t in thetas) / n_obs
        if n_obs >= 2:
            var = sum((t - theta_mean) ** 2 for t in thetas) / (n_obs - 1)
            theta_sd = var**0.5
            max_pair = max(
                abs(thetas[i] - thetas[j])
                for i in range(n_obs)
                for j in range(i + 1, n_obs)
            )
        else:
            theta_sd = 0.0
            max_pair = 0.0
        rows.append(
            {
                "serotype": str(serotype),
                "canon_label": str(canon),
                "chain": chain_from_canon(canon),
                "domain_label": domain.get((str(serotype), str(canon)), ""),
                "n_obs": int(n_obs),
                "theta_mean": round(theta_mean, ROUND_DECIMALS),
                "theta_sd": round(theta_sd, ROUND_DECIMALS),
                "theta_min": round(theta_min, ROUND_DECIMALS),
                "theta_max": round(theta_max, ROUND_DECIMALS),
                "theta_range": round(theta_range, ROUND_DECIMALS),
                "abs_theta_mean": round(abs_mean, ROUND_DECIMALS),
                "max_pairwise_abs_diff": round(max_pair, ROUND_DECIMALS),
                "in_all_replicates": bool(
                    in_all.get((str(serotype), str(canon)), False)
                ),
                "tier": TIER_EXPLORATORY,
            }
        )

    out = pd.DataFrame(rows, columns=list(REPLICATE_EFFECT_SPREAD_COLUMNS))
    return out.reset_index(drop=True)


def _in_all_lookup(replicate_inventory: pd.DataFrame) -> dict[tuple[str, str], bool]:
    if replicate_inventory.empty:
        return {}
    inv = replicate_inventory
    return {
        (str(s), str(c)): bool(a)
        for s, c, a in zip(
            inv["serotype"],
            inv["canon_label"],
            inv["in_all_replicates"],
            strict=True,
        )
    }


def _domain_lookup(
    replicate_table: pd.DataFrame | None,
) -> dict[tuple[str, str], str]:
    if replicate_table is None or replicate_table.empty:
        return {}
    if "domain_label" not in replicate_table.columns:
        return {}
    out: dict[tuple[str, str], str] = {}
    for s, c, d in zip(
        replicate_table["serotype"],
        replicate_table["canon_label"],
        replicate_table["domain_label"],
        strict=True,
    ):
        key = (str(s), str(c))
        if key not in out and pd.notna(d):
            out[key] = str(d)
    return out
