"""Per-serotype scorecard, swept over the ρ* band.

One row per (serotype, rho_star): the reproducibility scorecard the design's
per-serotype summary table (T1) and cross-serotype scorecard (§3.6) roll up —
locus count, resolution-census rollups (how many loci gate at residue vs
domain-or-coarser vs nowhere), mechanism direction counts (signed vs mixed),
the signed-significant count at ρ*, and the ρ distribution (median/IQR over
residue-scale ρ).

Everything is reported as a function of ρ* because the gate is uncalibrated
(§0.1, §5.3). Mechanism direction counts (signed/mixed) are ρ*-independent — a
mechanism's own direction does not change with ρ* — but the signed-significant
count does, because the screen includes ρ ≥ ρ*.

Pure: no IO, no mutation of inputs. Built from the already-constructed census,
residue-landscape, and signed-screen tables so the scorecard cannot disagree
with them.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DIRECTION_MIXED,
    DOMAIN_SCALE_LEVEL,
    PROVISIONAL_RHO_STAR,
    RESIDUE_SCALE_LEVEL,
    RHO_STAR_DECIMALS,
    SCALE_LEVEL_TO_INDEX,
    SEROTYPE_SUMMARY_COLUMNS,
)
from ._frames import gated_slice


def build_serotype_summary(
    stride_table: pd.DataFrame,
    resolution_census: pd.DataFrame,
    residue_landscape: pd.DataFrame,
    signed_screen: pd.DataFrame,
    rho_star_band: tuple[float, ...],
    rho_star_decimals: int = RHO_STAR_DECIMALS,
) -> pd.DataFrame:
    """Build the per-serotype scorecard over the ρ* band.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (used only for the gated-mechanism direction counts).
    resolution_census
        The S2 resolution census (per serotype, ρ*, gated scale).
    residue_landscape
        The S2 residue landscape (per serotype residue, for the ρ distribution).
    signed_screen
        The S2 signed screen (per serotype, residue, ρ*), for the
        signed-significant count.
    rho_star_band
        The already-normalised, ascending ρ* band.
    rho_star_decimals
        Decimals ρ* is rounded to when emitted.

    Returns
    -------
    DataFrame
        One row per (serotype, rho_star), column-ordered per
        :data:`~stride_s2.models.schema.SEROTYPE_SUMMARY_COLUMNS`.
    """
    if stride_table.empty or not rho_star_band:
        return pd.DataFrame(columns=list(SEROTYPE_SUMMARY_COLUMNS))

    provisional = round(PROVISIONAL_RHO_STAR, rho_star_decimals)
    domain_index = SCALE_LEVEL_TO_INDEX[DOMAIN_SCALE_LEVEL]

    serotypes = sorted(stride_table["serotype"].astype(str).unique().tolist())

    # ρ-distribution facts per serotype (ρ*-independent): residue-scale ρ.
    rho_stats = _residue_rho_stats(residue_landscape)

    # direction counts per serotype (ρ*-independent): from the gated mechanisms.
    gated = gated_slice(stride_table)
    dir_counts = _direction_counts(gated)

    # census rollups keyed by (serotype, rho_star)
    census_roll = _census_rollups(resolution_census, domain_index)

    # signed-significant count keyed by (serotype, rho_star)
    signif = _signed_significant_counts(signed_screen)

    records = []
    for serotype in serotypes:
        n_signed, n_mixed = dir_counts.get(serotype, (0, 0))
        n_mech = n_signed + n_mixed
        median, q1, q3, rmin, rmax = rho_stats.get(
            serotype, (float("nan"),) * 5
        )
        for rho_star in rho_star_band:
            roll = census_roll.get((serotype, rho_star), (0, 0, 0, 0))
            n_loci, n_res_gated, n_dom_gated, n_unresolved = roll
            records.append(
                {
                    "serotype": serotype,
                    "rho_star": rho_star,
                    "is_provisional_rho_star": rho_star == provisional,
                    "n_loci": n_loci,
                    "n_gated_residue": n_res_gated,
                    "n_gated_domain_or_coarser": n_dom_gated,
                    "n_unresolved": n_unresolved,
                    "frac_gated_residue": (
                        n_res_gated / n_loci if n_loci else 0.0
                    ),
                    "n_mechanisms": n_mech,
                    "n_signed": n_signed,
                    "n_mixed": n_mixed,
                    "n_signed_significant": signif.get(
                        (serotype, rho_star), 0
                    ),
                    "frac_mixed": n_mixed / n_mech if n_mech else 0.0,
                    "rho_residue_median": median,
                    "rho_residue_q1": q1,
                    "rho_residue_q3": q3,
                    "rho_residue_min": rmin,
                    "rho_residue_max": rmax,
                }
            )

    out = pd.DataFrame.from_records(records)[list(SEROTYPE_SUMMARY_COLUMNS)]
    out = out.sort_values(["serotype", "rho_star"]).reset_index(drop=True)
    if out.duplicated(["serotype", "rho_star"]).any():
        raise ConsistencyError(
            "serotype summary has duplicate (serotype, rho_star) rows"
        )
    return out


def _residue_rho_stats(
    residue_landscape: pd.DataFrame,
) -> dict[str, tuple[float, float, float, float, float]]:
    """Median / Q1 / Q3 / min / max of residue-scale ρ per serotype."""
    stats: dict[str, tuple[float, float, float, float, float]] = {}
    if residue_landscape.empty:
        return stats
    for serotype, grp in residue_landscape.groupby("serotype", sort=True):
        rho = grp["rho_residue"].dropna()
        if rho.empty:
            stats[str(serotype)] = (float("nan"),) * 5
            continue
        stats[str(serotype)] = (
            float(rho.median()),
            float(rho.quantile(0.25)),
            float(rho.quantile(0.75)),
            float(rho.min()),
            float(rho.max()),
        )
    return stats


def _direction_counts(gated: pd.DataFrame) -> dict[str, tuple[int, int]]:
    """(n_signed, n_mixed) mechanism counts per serotype."""
    counts: dict[str, tuple[int, int]] = {}
    if gated.empty:
        return counts
    for serotype, grp in gated.groupby("serotype", sort=True):
        directions = grp["mech_direction"].astype("string")
        n_mixed = int((directions == DIRECTION_MIXED).sum())
        n_signed = int(len(grp) - n_mixed)
        counts[str(serotype)] = (n_signed, n_mixed)
    return counts


def _census_rollups(
    resolution_census: pd.DataFrame, domain_index: int
) -> dict[tuple[str, float], tuple[int, int, int, int]]:
    """(n_loci, n_residue_gated, n_domain_or_coarser, n_unresolved) per key."""
    roll: dict[tuple[str, float], tuple[int, int, int, int]] = {}
    if resolution_census.empty:
        return roll
    residue_index = SCALE_LEVEL_TO_INDEX[RESIDUE_SCALE_LEVEL]
    grouped = resolution_census.groupby(["serotype", "rho_star"], sort=True)
    for (serotype, rho_star), grp in grouped:
        n_loci = int(grp["n_loci"].sum())
        n_res = int(
            grp.loc[grp["gated_scale_index"] == residue_index, "n_loci"].sum()
        )
        n_dom = int(
            grp.loc[grp["gated_scale_index"] >= domain_index, "n_loci"].sum()
        )
        n_unres = int(grp.loc[grp["gated_scale_index"] < 0, "n_loci"].sum())
        roll[(str(serotype), float(rho_star))] = (n_loci, n_res, n_dom, n_unres)
    return roll


def _signed_significant_counts(
    signed_screen: pd.DataFrame,
) -> dict[tuple[str, float], int]:
    """Count of mechanisms passing the screen per (serotype, rho_star)."""
    signif: dict[tuple[str, float], int] = {}
    if signed_screen.empty:
        return signif
    passed = signed_screen[signed_screen["passes_screen"].astype(bool)]
    grouped = passed.groupby(["serotype", "rho_star"], sort=True).size()
    for (serotype, rho_star), n in grouped.items():
        signif[(str(serotype), float(rho_star))] = int(n)
    return signif
