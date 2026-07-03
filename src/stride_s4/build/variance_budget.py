"""Variance-component budget per domain (Tier A — licensed).

One row per (serotype, chain, domain): the domain-scale variance components τ²
(replicate disagreement) and σ̄² (sampling noise), split into fractions of the
unreproduced variance with a regime label. This is the design's
variance-component budgeting and τ²/σ̄² ratio diagnostic at the licensed domain
scale (§3.1, §3.5, table T5).

Pure: no IO, no mutation of inputs. The variance components are read from the
domain region's domain-scale profile row, never recomputed.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    RHO_DECIMALS,
    TIER_LICENSED,
    VARIANCE_BUDGET_COLUMNS,
)
from ._frames import domain_regions
from ._stats import variance_fractions


def build_variance_budget(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the per-domain variance-component budget table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).

    Returns
    -------
    DataFrame
        One row per (serotype, chain, domain), column-ordered per
        :data:`~stride_s4.models.schema.VARIANCE_BUDGET_COLUMNS`.
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(VARIANCE_BUDGET_COLUMNS))

    regions = domain_regions(stride_table)
    if regions.empty:
        return pd.DataFrame(columns=list(VARIANCE_BUDGET_COLUMNS))

    records = []
    for row in regions.itertuples(index=False):
        total, frac_tau2, frac_sigma2, ratio, regime = variance_fractions(
            float(row.tau2), float(row.sigma2_bar)
        )
        records.append(
            {
                "serotype": row.serotype,
                "chain": row.chain,
                "domain": row.domain,
                "region_id": row.region_id,
                "rho_domain": _round(float(row.rho)),
                "beta_domain": _round(float(row.beta)),
                "beta_se_domain": _round(float(row.beta_se)),
                "tau2": _round(float(row.tau2)),
                "sigma2_bar": _round(float(row.sigma2_bar)),
                "total_unreproduced": _round(total),
                "frac_tau2": _round(frac_tau2),
                "frac_sigma2": _round(frac_sigma2),
                "tau2_sigma2_ratio": _round(ratio),
                "variance_regime": regime,
                "tier": TIER_LICENSED,
            }
        )

    out = pd.DataFrame.from_records(records)[list(VARIANCE_BUDGET_COLUMNS)]
    out = out.sort_values(["serotype", "chain", "domain"]).reset_index(drop=True)
    if out.duplicated(["serotype", "chain", "domain"]).any():
        raise ConsistencyError(
            "variance budget has duplicate (serotype, chain, domain) rows"
        )
    return out


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
