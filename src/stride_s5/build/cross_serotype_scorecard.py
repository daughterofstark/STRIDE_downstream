"""Per-serotype cross-serotype scorecard (per serotype).

One row per serotype: the cross-serotype reproducibility scorecard of the design
(§3.6) — the locus count, the residue-gated fraction, the signed / mixed
mechanism composition, the residue-scale ρ median/min/max, and how many of the
pan-serotype shared positions are reproducible in this serotype. This is the
single per-serotype table that summarises each system's cross-serotype-relevant
profile.

Its core quantities are residue-scale (or per-mechanism) descriptions, so every
row is Tier B (exploratory) and stamped with the provisional ρ\\*. Mechanisms are
counted **one per gated region** (loci sharing a coarse gated region collapse to
one mechanism), mirroring the mechanism payload; ρ is read, never recomputed.

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    CROSS_SEROTYPE_SCORECARD_COLUMNS,
    DIRECTION_DECREASE,
    DIRECTION_INCREASE,
    PROVISIONAL_RHO_STAR,
    RHO_DECIMALS,
    TIER_EXPLORATORY,
)
from ._frames import position_frame, residue_slice

_SIGNED = frozenset({DIRECTION_INCREASE, DIRECTION_DECREASE})


def build_cross_serotype_scorecard(
    stride_table: pd.DataFrame,
    conservation_table: pd.DataFrame,
    rho_star: float = PROVISIONAL_RHO_STAR,
) -> pd.DataFrame:
    """Build the per-serotype cross-serotype scorecard.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).
    conservation_table
        The S1A ``conservation_table`` — used to identify the pan-serotype shared
        positions (``in_all_serotypes``).
    rho_star
        The provisional gate threshold (default :data:`PROVISIONAL_RHO_STAR`).

    Returns
    -------
    DataFrame
        One row per serotype, column-ordered per
        :data:`~stride_s5.models.schema.CROSS_SEROTYPE_SCORECARD_COLUMNS`.
    """
    pf = position_frame(stride_table, rho_star)
    if pf.empty:
        return pd.DataFrame(columns=list(CROSS_SEROTYPE_SCORECARD_COLUMNS))

    res = residue_slice(stride_table)
    mechanisms = _mechanism_frame(stride_table)
    shared_labels = _shared_labels(conservation_table)

    records = []
    for serotype, grp in pf.groupby("serotype", sort=True):
        serotype = str(serotype)
        n_loci = int(len(grp))
        n_repro = int(grp["reproducible"].astype(bool).sum())
        frac_repro = (n_repro / n_loci) if n_loci > 0 else float("nan")

        mech = mechanisms[mechanisms["serotype"] == serotype]
        n_mech = int(len(mech))
        n_signed = int(mech["direction"].isin(_SIGNED).sum())
        n_mixed = n_mech - n_signed
        frac_signed = (n_signed / n_mech) if n_mech > 0 else float("nan")
        frac_mixed = (n_mixed / n_mech) if n_mech > 0 else float("nan")

        rho_vals = res[res["serotype"] == serotype]["rho"].astype(float)

        shared_here = grp[grp["canon_label"].isin(shared_labels)]
        n_shared = int(len(shared_here))
        n_shared_repro = int(shared_here["reproducible"].astype(bool).sum())

        records.append(
            {
                "serotype": serotype,
                "n_loci": n_loci,
                "n_reproducible_residue": n_repro,
                "frac_reproducible_residue": _round(frac_repro),
                "n_mechanisms": n_mech,
                "n_signed": n_signed,
                "n_mixed": n_mixed,
                "frac_signed": _round(frac_signed),
                "frac_mixed": _round(frac_mixed),
                "rho_residue_median": _round(_series_stat(rho_vals, "median")),
                "rho_residue_min": _round(_series_stat(rho_vals, "min")),
                "rho_residue_max": _round(_series_stat(rho_vals, "max")),
                "n_shared_positions": n_shared,
                "n_shared_reproducible": n_shared_repro,
                "rho_star": rho_star,
                "is_provisional_rho_star": True,
                "tier": TIER_EXPLORATORY,
            }
        )

    out = pd.DataFrame.from_records(records)[
        list(CROSS_SEROTYPE_SCORECARD_COLUMNS)
    ]
    out = out.sort_values(["serotype"]).reset_index(drop=True)
    if out.duplicated(["serotype"]).any():
        raise ConsistencyError(
            "cross serotype scorecard has duplicate serotype rows"
        )
    return out


def _mechanism_frame(stride_table: pd.DataFrame) -> pd.DataFrame:
    """One row per gated region (mechanism): (serotype, region_id, direction).

    Loci that gate into the same coarse region share one mechanism, so the gated
    rows are de-duplicated by (serotype, region_id).
    """
    gated = stride_table[stride_table["is_gated_scale"].astype(bool)]
    if gated.empty:
        return pd.DataFrame(columns=["serotype", "region_id", "direction"])
    records = []
    seen: set[tuple[str, str]] = set()
    for row in gated.sort_values(["serotype", "region_id"]).itertuples(
        index=False
    ):
        key = (str(row.serotype), str(row.region_id))
        if key in seen:
            continue
        seen.add(key)
        direction = (
            str(row.mech_direction)
            if pd.notna(row.mech_direction)
            else "mixed"
        )
        records.append(
            {
                "serotype": key[0],
                "region_id": key[1],
                "direction": direction,
            }
        )
    return pd.DataFrame.from_records(records)


def _shared_labels(conservation_table: pd.DataFrame) -> set[str]:
    if conservation_table.empty:
        return set()
    shared = conservation_table[
        conservation_table["in_all_serotypes"].astype(bool)
    ]
    return {str(x) for x in shared["canon_label"]}


def _series_stat(series: pd.Series, stat: str) -> float:
    if series.empty:
        return float("nan")
    if stat == "median":
        return float(series.median())
    if stat == "min":
        return float(series.min())
    return float(series.max())


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
