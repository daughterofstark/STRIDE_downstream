"""Shared, pure helpers for the S7 builders.

All helpers read their inputs and return new frames without mutating them. They
select / order / round only — they never compute a new statistic.
"""
from __future__ import annotations

import pandas as pd

from ..models.schema import PROVISIONAL_RHO_STAR, ROUND_DECIMALS


def provisional_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Rows evaluated at the provisional ρ\\* (the single-threshold view).

    Prefers the explicit ``is_provisional_rho_star`` flag; falls back to matching
    ``rho_star`` against :data:`PROVISIONAL_RHO_STAR` when only that column exists.
    """
    if "is_provisional_rho_star" in df.columns:
        out = df[df["is_provisional_rho_star"].astype(bool)]
        if not out.empty:
            return out.copy()
    if "rho_star" in df.columns:
        return df[df["rho_star"] == PROVISIONAL_RHO_STAR].copy()
    return df.copy()


def round_floats(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Round the named float columns to :data:`ROUND_DECIMALS` (out-of-place)."""
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").round(
                ROUND_DECIMALS
            )
    return out


def select_columns(df: pd.DataFrame, columns: tuple[str, ...]) -> pd.DataFrame:
    """Return exactly ``columns`` in order, preserving row order."""
    return df.loc[:, list(columns)].reset_index(drop=True)


def sort_by(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Stable-sort by ``keys`` (only those present), for deterministic output."""
    present = [k for k in keys if k in df.columns]
    if not present:
        return df.reset_index(drop=True)
    return df.sort_values(present, kind="stable").reset_index(drop=True)
