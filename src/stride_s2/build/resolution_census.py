"""Achieved-resolution census, swept over the ρ* band.

For each ρ* in the band, every locus is re-gated to its finest scale at which
ρ ≥ ρ* (design §1.5, §5.3), and the loci are counted per gated scale level. The
result is the achieved-resolution census as a function of ρ* — never at a single
threshold, because the gate is uncalibrated (§0.1).

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import (
    PROVISIONAL_RHO_STAR,
    RESOLUTION_CENSUS_COLUMNS,
    SCALE_LEVEL_TO_INDEX,
    SCALE_UNRESOLVED,
)
from ._frames import rho_by_scale_map
from ._screens import regate_scale, scale_tier


def build_resolution_census(
    stride_table: pd.DataFrame,
    rho_star_band: tuple[float, ...],
    rho_star_decimals: int,
) -> pd.DataFrame:
    """Build the resolution census over the ρ* band.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).
    rho_star_band
        The already-normalised, ascending ρ* band to sweep.
    rho_star_decimals
        Decimals ρ* values are rounded to when emitted (kept consistent with the
        band so the ``is_provisional_rho_star`` comparison is exact).

    Returns
    -------
    DataFrame
        One row per (serotype, rho_star, gated_scale_level) with a non-zero
        locus count, column-ordered per
        :data:`~stride_s2.models.schema.RESOLUTION_CENSUS_COLUMNS`.
    """
    if stride_table.empty or not rho_star_band:
        return pd.DataFrame(columns=list(RESOLUTION_CENSUS_COLUMNS))

    rho_maps = rho_by_scale_map(stride_table)
    provisional = round(PROVISIONAL_RHO_STAR, rho_star_decimals)

    # (serotype, rho_star, scale_level) -> count
    counts: dict[tuple[str, float, str], int] = {}
    for (serotype, _canon_label), rho_by_scale in rho_maps.items():
        for rho_star in rho_star_band:
            level, _idx = regate_scale(rho_by_scale, rho_star)
            key = (serotype, rho_star, level)
            counts[key] = counts.get(key, 0) + 1

    records = []
    for (serotype, rho_star, level), n in counts.items():
        idx = -1 if level == SCALE_UNRESOLVED else SCALE_LEVEL_TO_INDEX[level]
        records.append(
            {
                "serotype": serotype,
                "rho_star": rho_star,
                "is_provisional_rho_star": rho_star == provisional,
                "gated_scale_level": level,
                "gated_scale_index": idx,
                "tier": scale_tier(level),
                "n_loci": int(n),
            }
        )

    out = pd.DataFrame.from_records(records)[list(RESOLUTION_CENSUS_COLUMNS)]
    return out.sort_values(
        ["serotype", "rho_star", "gated_scale_index"]
    ).reset_index(drop=True)
