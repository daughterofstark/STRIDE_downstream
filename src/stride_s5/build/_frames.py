r"""Internal frame-extraction helpers shared by the S5 builders.

These pull the reusable views S5 needs out of the S0 STRIDE table and reduce
them, per the design's anti-pseudoreplication rule, to **one value per serotype
per position/region first** (§5.2) — before any cross-serotype comparison:

- :func:`residue_slice` — the residue-scale rows (one per locus), carrying ρ and
  the structural labels;
- :func:`gated_directions` — the gated row of each locus (one per emitted
  mechanism), carrying the mechanism direction;
- :func:`position_frame` — the per-(serotype, canon_label) cross-serotype-ready
  frame: residue-scale ρ, the descriptive "reproducible" flag (ρ ≥ ρ\*), and the
  signed direction of the position's gated mechanism;
- :func:`domain_regions` — one row per (serotype, chain, domain) with the
  region-constant domain-scale ρ/β/variance components.

Pure functions; the input is never mutated and no IO happens here. ρ and the
variance components are read, never recomputed (§1.2).
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DIRECTION_DECREASE,
    DIRECTION_INCREASE,
    DIRECTION_MIXED,
    DIRECTION_NONE,
    DOMAIN_SCALE_INDEX,
    DOMAIN_SCALE_LEVEL,
    RESIDUE_SCALE_LEVEL,
    SCALE_LEVEL_TO_INDEX,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)

_SIGNED_DIRECTIONS = frozenset({DIRECTION_INCREASE, DIRECTION_DECREASE})


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


def gated_directions(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Per (serotype, canon_label), the gated row's scale and mechanism direction.

    Exactly one row per locus carries ``is_gated_scale == True`` (the locus's
    gated scale ℓ̂*); this returns its ``scale_level`` and ``mech_direction``
    (``mixed`` when the mechanism carries no coherent sign).

    Raises
    ------
    ConsistencyError
        If a (serotype, canon_label) carries more than one gated row.
    """
    gated = stride_table[stride_table["is_gated_scale"].astype(bool)].copy()
    if gated.empty:
        return pd.DataFrame(
            columns=["serotype", "canon_label", "gated_scale_level", "direction"]
        )
    dup = gated.duplicated(["serotype", "canon_label"], keep=False)
    if dup.any():
        examples = (
            gated.loc[dup, ["serotype", "canon_label"]]
            .drop_duplicates()
            .head(3)
            .to_dict("records")
        )
        raise ConsistencyError(
            "STRIDE table has multiple gated rows for "
            f"(serotype, canon_label); examples: {examples}"
        )
    records = []
    for row in gated.itertuples(index=False):
        direction = (
            str(row.mech_direction)
            if pd.notna(row.mech_direction)
            else DIRECTION_MIXED
        )
        records.append(
            {
                "serotype": str(row.serotype),
                "canon_label": str(row.canon_label),
                "gated_scale_level": str(row.scale_level),
                "direction": direction,
            }
        )
    return pd.DataFrame.from_records(records)


def position_frame(
    stride_table: pd.DataFrame, rho_star: float
) -> pd.DataFrame:
    """One row per (serotype, canon_label): the cross-serotype-ready position view.

    Combines the residue-scale ρ (read, never recomputed) with the position's
    gated mechanism direction, and derives the descriptive reproducibility and
    signed flags at the provisional gate. This is the per-serotype aggregation the
    cross-serotype builders compare across serotypes (§5.2 — aggregate first).

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).
    rho_star
        The provisional gate threshold; a position is "reproducible" at the
        residue scale when its residue-scale ρ ≥ ``rho_star`` (equivalently, it is
        gated at the residue scale).

    Returns
    -------
    DataFrame
        Columns ``serotype, canon_label, chain, domain, rho_residue,
        reproducible, gated_scale_level, direction, is_signed``.
    """
    res = residue_slice(stride_table)
    if res.empty:
        return pd.DataFrame(
            columns=[
                "serotype",
                "canon_label",
                "chain",
                "domain",
                "rho_residue",
                "reproducible",
                "gated_scale_level",
                "direction",
                "is_signed",
            ]
        )
    directions = gated_directions(stride_table)
    dir_by_key: dict[tuple[str, str], tuple[str, str]] = {}
    for row in directions.itertuples(index=False):
        dir_by_key[(str(row.serotype), str(row.canon_label))] = (
            str(row.gated_scale_level),
            str(row.direction),
        )

    records = []
    for row in res.itertuples(index=False):
        serotype = str(row.serotype)
        canon = str(row.canon_label)
        rho_residue = float(row.rho)
        reproducible = rho_residue >= rho_star
        gated_scale_level, direction = dir_by_key.get(
            (serotype, canon), ("", DIRECTION_MIXED)
        )
        # A signed residue-scale claim requires reproducibility at the residue
        # scale (so the gated mechanism is the residue-scale mechanism) and a
        # coherent direction. Positions that gate coarser than residue are not
        # reproducible here, so they carry no residue-scale signed claim.
        if reproducible and direction in _SIGNED_DIRECTIONS:
            is_signed = True
            signed_direction = direction
        else:
            is_signed = False
            signed_direction = (
                direction if reproducible else DIRECTION_NONE
            )
        records.append(
            {
                "serotype": serotype,
                "canon_label": canon,
                "chain": str(row.h_chain),
                "domain": str(row.h_domain),
                "rho_residue": rho_residue,
                "reproducible": bool(reproducible),
                "gated_scale_level": gated_scale_level,
                "direction": signed_direction,
                "is_signed": bool(is_signed),
            }
        )
    out = pd.DataFrame.from_records(records)
    return out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)


def domain_slice(stride_table: pd.DataFrame) -> pd.DataFrame:
    """The domain-scale rows, sorted deterministically."""
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
        return pd.DataFrame(
            columns=[
                "serotype",
                "chain",
                "domain",
                "region_id",
                "rho",
                "beta",
                "beta_se",
                "tau2",
                "sigma2_bar",
            ]
        )
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


def scale_tier(scale_level: str) -> str:
    """Classify a scale level into the licensed / exploratory tier.

    Domain scale and coarser is *licensed*; anything finer than domain is
    *exploratory* (design §5.3–5.4).
    """
    idx = SCALE_LEVEL_TO_INDEX.get(scale_level)
    if idx is None:
        return TIER_EXPLORATORY
    return TIER_LICENSED if idx >= DOMAIN_SCALE_INDEX else TIER_EXPLORATORY


def _assert_region_constant(
    grp: pd.DataFrame, serotype: str, chain: str, domain: str, col: str
) -> None:
    vals = grp[col].dropna().unique()
    if len(vals) > 1:
        raise ConsistencyError(
            f"domain region ({serotype}, {chain}, {domain}) has non-constant "
            f"{col} at domain scale: {sorted(vals.tolist())}"
        )
