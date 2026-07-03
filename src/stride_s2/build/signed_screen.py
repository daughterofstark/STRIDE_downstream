"""Signed / significant screen over the gated mechanisms, swept over the band.

One row per (serotype, canon_label, rho_star) over the emitted mechanisms (the
gated rows). A mechanism yields a signed, significant, reproducible claim when
its direction is not ``mixed`` **and** its β confidence interval excludes 0
**and** its ρ meets ρ* (design §2.2, §3.5). The screen is evaluated across the
whole ρ* band so significance is reported as a function of ρ*, never at a single
threshold (§5.3).

Mixed mechanisms still carry a real unsigned magnitude but no signed direction,
so they never pass the screen — this is the common case in the real data and is
handled explicitly (§1.3).

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    PROVISIONAL_RHO_STAR,
    RHO_STAR_DECIMALS,
    SIGNED_SCREEN_COLUMNS,
)
from ._frames import gated_slice
from ._screens import ci_excludes_zero, is_signed, passes_signed_screen, scale_tier


def build_signed_screen(
    stride_table: pd.DataFrame,
    rho_star_band: tuple[float, ...],
    rho_star_decimals: int = RHO_STAR_DECIMALS,
) -> pd.DataFrame:
    """Build the signed/significant screen over the ρ* band.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales); the gated rows carry the mechanisms.
    rho_star_band
        The already-normalised, ascending ρ* band to sweep.
    rho_star_decimals
        Decimals ρ* is rounded to when emitted.

    Returns
    -------
    DataFrame
        One row per (serotype, canon_label, rho_star), column-ordered per
        :data:`~stride_s2.models.schema.SIGNED_SCREEN_COLUMNS`.
    """
    if stride_table.empty or not rho_star_band:
        return pd.DataFrame(columns=list(SIGNED_SCREEN_COLUMNS))

    gated = gated_slice(stride_table)
    if gated.empty:
        return pd.DataFrame(columns=list(SIGNED_SCREEN_COLUMNS))

    provisional = round(PROVISIONAL_RHO_STAR, rho_star_decimals)

    records = []
    for row in gated.itertuples(index=False):
        direction = _s(row.mech_direction)
        beta_signed = _opt_f(row.mech_beta_signed)
        ci_lo = _opt_f(row.mech_beta_ci_lower)
        ci_hi = _opt_f(row.mech_beta_ci_upper)
        rho = float(row.rho)
        signed = is_signed(direction)
        excludes_zero = ci_excludes_zero(ci_lo, ci_hi)
        gated_level = str(row.scale_level)
        for rho_star in rho_star_band:
            records.append(
                {
                    "serotype": row.serotype,
                    "canon_label": row.canon_label,
                    "rho_star": rho_star,
                    "is_provisional_rho_star": rho_star == provisional,
                    "gated_scale_level": gated_level,
                    "tier": scale_tier(gated_level),
                    "rho": rho,
                    "direction": direction if direction is not None else "",
                    "is_signed": signed,
                    "beta_signed": (
                        beta_signed
                        if beta_signed is not None
                        else float("nan")
                    ),
                    "beta_ci_lower": (
                        ci_lo if ci_lo is not None else float("nan")
                    ),
                    "beta_ci_upper": (
                        ci_hi if ci_hi is not None else float("nan")
                    ),
                    "ci_excludes_zero": excludes_zero,
                    "meets_rho_star": bool(rho >= rho_star),
                    "passes_screen": passes_signed_screen(
                        direction, ci_lo, ci_hi, rho, rho_star
                    ),
                }
            )

    out = pd.DataFrame.from_records(records)[list(SIGNED_SCREEN_COLUMNS)]
    out = out.sort_values(
        ["serotype", "canon_label", "rho_star"]
    ).reset_index(drop=True)
    if out.duplicated(["serotype", "canon_label", "rho_star"]).any():
        raise ConsistencyError(
            "signed screen has duplicate (serotype, canon_label, rho_star) rows"
        )
    return out


def _s(value: object) -> str | None:
    """Coerce a string cell to str, mapping pandas NA/NaN to None."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if pd.isna(value):
        return None
    return str(value)


def _opt_f(value: object) -> float | None:
    """Coerce a numeric cell to float, mapping NaN/None to None."""
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    if pd.isna(value):
        return None
    return float(value)  # type: ignore[arg-type]
