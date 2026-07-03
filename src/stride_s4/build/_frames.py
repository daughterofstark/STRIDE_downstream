"""Internal frame-extraction helpers shared by the S4 builders.

These pull the reusable views S4 needs out of the S0 STRIDE table — the
residue-scale slice (one row per locus), the domain-scale slice (region-level
constants shared by member loci), and the gated-mechanism slice (one row per
emitted mechanism) — with the input never mutated. Pure functions; no IO.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DOMAIN_SCALE_LEVEL,
    RESIDUE_SCALE_LEVEL,
    SCALE_LEVEL_TO_INDEX,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)


def residue_slice(stride_table: pd.DataFrame) -> pd.DataFrame:
    """The residue-scale rows (one per locus), sorted deterministically.

    Raises
    ------
    ConsistencyError
        If a (serotype, canon_label) carries more than one residue-scale row.
    """
    res = stride_table[
        stride_table["scale_level"] == RESIDUE_SCALE_LEVEL
    ].copy()
    dup = res.duplicated(["serotype", "canon_label"], keep=False)
    if dup.any():
        examples = (
            res.loc[dup, ["serotype", "canon_label"]]
            .drop_duplicates()
            .head(3)
            .to_dict("records")
        )
        raise ConsistencyError(
            "STRIDE table has duplicate residue-scale rows for "
            f"(serotype, canon_label); examples: {examples}"
        )
    return res.sort_values(["serotype", "canon_label"]).reset_index(drop=True)


def domain_slice(stride_table: pd.DataFrame) -> pd.DataFrame:
    """The domain-scale rows, sorted deterministically.

    At domain scale a region's ρ/β/variance components are region-level constants
    shared by every member locus, so multiple loci may carry identical domain
    rows; :func:`domain_regions` collapses them.
    """
    dom = stride_table[
        stride_table["scale_level"] == DOMAIN_SCALE_LEVEL
    ].copy()
    return dom.sort_values(
        ["serotype", "h_chain", "h_domain", "canon_label"]
    ).reset_index(drop=True)


def domain_regions(stride_table: pd.DataFrame) -> pd.DataFrame:
    """One row per (serotype, chain, domain) with region-constant domain values.

    Collapses the domain-scale slice to a single row per domain region, verifying
    that ρ, β, β_se, τ² and σ̄² are constant across the region's member loci (the
    region-level invariant the design relies on, §1.4).

    Raises
    ------
    ConsistencyError
        If any of those quantities is not constant within a domain region.
    """
    dom = domain_slice(stride_table)
    if dom.empty:
        return dom
    records = []
    for keys, grp in dom.groupby(["serotype", "h_chain", "h_domain"], sort=True):
        serotype, chain, domain = (str(k) for k in keys)
        for col in ("rho", "beta", "beta_se", "tau2", "sigma2_bar"):
            _assert_region_constant(grp, serotype, chain, domain, col)
        first = grp.iloc[0]
        records.append(
            {
                "serotype": serotype,
                "chain": chain,
                "domain": domain,
                "region_id": str(first["region_id"]),
                "rho": float(first["rho"]),
                "beta": float(first["beta"]),
                "beta_se": float(first["beta_se"]),
                "tau2": float(first["tau2"]),
                "sigma2_bar": float(first["sigma2_bar"]),
            }
        )
    return pd.DataFrame.from_records(records)


def gated_slice(stride_table: pd.DataFrame) -> pd.DataFrame:
    """The gated rows (one per emitted mechanism), sorted deterministically."""
    gated = stride_table[stride_table["is_gated_scale"].astype(bool)].copy()
    return gated.sort_values(["serotype", "canon_label"]).reset_index(drop=True)


def scale_tier(scale_level: str) -> str:
    """Classify a scale level into the licensed / exploratory tier.

    Domain scale and coarser is *licensed*; anything finer than domain is
    *exploratory* (design §5.3–5.4).
    """
    domain_index = SCALE_LEVEL_TO_INDEX[DOMAIN_SCALE_LEVEL]
    idx = SCALE_LEVEL_TO_INDEX.get(scale_level)
    if idx is None:
        return TIER_EXPLORATORY
    return TIER_LICENSED if idx >= domain_index else TIER_EXPLORATORY


def _assert_region_constant(
    grp: pd.DataFrame, serotype: str, chain: str, domain: str, col: str
) -> None:
    vals = grp[col].dropna().unique()
    if len(vals) > 1:
        raise ConsistencyError(
            f"domain region ({serotype}, {chain}, {domain}) has non-constant "
            f"{col} at domain scale: {sorted(vals.tolist())}"
        )
