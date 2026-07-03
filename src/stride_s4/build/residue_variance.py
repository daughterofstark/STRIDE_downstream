"""Per-residue variance decomposition and replicate-disagreement ranking (Tier B).

One row per (serotype, canon_label): the residue-scale variance components τ² and
σ̄² split into fractions with a regime label, plus ``tau2_rank`` — the position's
rank within its serotype by τ² descending. This is the design's
replicate-disagreement mapping ("rank positions by tau2 to find where the three
runs most disagree", §3.1) at the exploratory residue scale.

Pure: no IO, no mutation of inputs. The variance components are read from the
residue-scale profile row, never recomputed.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    RESIDUE_VARIANCE_COLUMNS,
    RHO_DECIMALS,
    TIER_EXPLORATORY,
)
from ._frames import residue_slice
from ._stats import variance_fractions


def build_residue_variance(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the per-residue variance-decomposition table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).

    Returns
    -------
    DataFrame
        One row per (serotype, canon_label), column-ordered per
        :data:`~stride_s4.models.schema.RESIDUE_VARIANCE_COLUMNS`.
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(RESIDUE_VARIANCE_COLUMNS))

    res = residue_slice(stride_table)
    if res.empty:
        return pd.DataFrame(columns=list(RESIDUE_VARIANCE_COLUMNS))

    records = []
    for row in res.itertuples(index=False):
        total, frac_tau2, frac_sigma2, ratio, regime = variance_fractions(
            float(row.tau2), float(row.sigma2_bar)
        )
        records.append(
            {
                "serotype": str(row.serotype),
                "canon_label": str(row.canon_label),
                "chain": str(row.h_chain),
                "domain": str(row.h_domain),
                "rho_residue": _round(float(row.rho)),
                "beta_residue": _round(float(row.beta)),
                "beta_se_residue": _round(float(row.beta_se)),
                "tau2": _round(float(row.tau2)),
                "sigma2_bar": _round(float(row.sigma2_bar)),
                "total_unreproduced": _round(total),
                "frac_tau2": _round(frac_tau2),
                "frac_sigma2": _round(frac_sigma2),
                "tau2_sigma2_ratio": _round(ratio),
                "variance_regime": regime,
                "tier": TIER_EXPLORATORY,
            }
        )

    out = pd.DataFrame.from_records(records)
    # rank positions within a serotype by τ² descending (rank 1 = largest τ²).
    # 'min' method keeps ties sharing a rank; deterministic given the stable sort.
    out["tau2_rank"] = (
        out.groupby("serotype")["tau2"]
        .rank(ascending=False, method="min")
        .astype("int64")
    )
    out = out[list(RESIDUE_VARIANCE_COLUMNS)]
    out = out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
    if out.duplicated(["serotype", "canon_label"]).any():
        raise ConsistencyError(
            "residue variance has duplicate (serotype, canon_label) rows"
        )
    return out


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
