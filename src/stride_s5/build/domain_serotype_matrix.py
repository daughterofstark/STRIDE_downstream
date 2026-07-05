"""Domain × serotype ρ matrix (Tier A — licensed).

One row per (serotype, chain, domain): the domain-scale ρ (and β, β_se, and the
variance components), read from the domain region's region-constant domain-scale
profile row. This is the tidy-long form of the design's ρ(domain × serotype)
heatmap over the NS3 domains + NS2B (§3.3, figure F3); a later stage pivots it to
the matrix/heatmap. The catalytic-machinery domains (Catalytic Triad, Oxyanion
Loop) are flagged so the conserved-catalytic-machinery check (§3.3) can read them
off directly.

Domain scale is the licensed claim level at K = 3, so every row is Tier A. ρ and
the variance components are read, never recomputed.

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DOMAIN_SEROTYPE_MATRIX_COLUMNS,
    RHO_DECIMALS,
    TIER_LICENSED,
)
from ._classify import is_catalytic_domain
from ._frames import domain_regions


def build_domain_serotype_matrix(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the tidy-long domain × serotype ρ matrix.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).

    Returns
    -------
    DataFrame
        One row per (serotype, chain, domain), column-ordered per
        :data:`~stride_s5.models.schema.DOMAIN_SEROTYPE_MATRIX_COLUMNS`.
    """
    regions = domain_regions(stride_table)
    if regions.empty:
        return pd.DataFrame(columns=list(DOMAIN_SEROTYPE_MATRIX_COLUMNS))

    records = []
    for row in regions.itertuples(index=False):
        records.append(
            {
                "serotype": str(row.serotype),
                "chain": str(row.chain),
                "domain": str(row.domain),
                "region_id": str(row.region_id),
                "rho_domain": _round(float(row.rho)),
                "beta_domain": _round(float(row.beta)),
                "beta_se_domain": _round(float(row.beta_se)),
                "tau2_domain": _round(float(row.tau2)),
                "sigma2_bar_domain": _round(float(row.sigma2_bar)),
                "is_catalytic_domain": is_catalytic_domain(str(row.domain)),
                "tier": TIER_LICENSED,
            }
        )

    out = pd.DataFrame.from_records(records)[
        list(DOMAIN_SEROTYPE_MATRIX_COLUMNS)
    ]
    out = out.sort_values(["serotype", "chain", "domain"]).reset_index(drop=True)
    if out.duplicated(["serotype", "chain", "domain"]).any():
        raise ConsistencyError(
            "domain serotype matrix has duplicate (serotype, chain, domain) rows"
        )
    return out


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
