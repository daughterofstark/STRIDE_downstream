"""Loaders for the S1A tables.

S1B consumes **only** the four S1A parquet tables. It never re-reads the S0
canonical tables, the raw STRIDE CSVs/JSONs, or MD trajectories. These loaders
read the parquet, assert the columns S1B depends on are present, and return the
DataFrames unchanged otherwise.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..models.errors import InputError
from ..models.schema import (
    CANONICAL_RESIDUES_REQUIRED,
    CONSERVATION_TABLE_REQUIRED,
    DOMAIN_TABLE_REQUIRED,
    REPLICATE_INVENTORY_REQUIRED,
)


def _read_parquet(path: Path, what: str) -> pd.DataFrame:
    if not path.is_file():
        raise InputError(f"{what} not found: {path}")
    try:
        return pd.read_parquet(path)
    except Exception as exc:  # unreadable / not parquet
        raise InputError(f"could not read {what} at {path}: {exc}") from exc


def _require_columns(
    df: pd.DataFrame, required: tuple[str, ...], what: str
) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise InputError(f"{what} is missing required column(s): {missing}")


def load_canonical_residues(path: str | Path) -> pd.DataFrame:
    """Load the S1A canonical residues table, asserting required columns."""
    df = _read_parquet(Path(path), "canonical_residues")
    _require_columns(df, CANONICAL_RESIDUES_REQUIRED, "canonical_residues")
    return df


def load_domain_table(path: str | Path) -> pd.DataFrame:
    """Load the S1A domain table, asserting required columns."""
    df = _read_parquet(Path(path), "domain_table")
    _require_columns(df, DOMAIN_TABLE_REQUIRED, "domain_table")
    return df


def load_replicate_inventory(path: str | Path) -> pd.DataFrame:
    """Load the S1A replicate inventory, asserting required columns."""
    df = _read_parquet(Path(path), "replicate_inventory")
    _require_columns(df, REPLICATE_INVENTORY_REQUIRED, "replicate_inventory")
    return df


def load_conservation_table(path: str | Path) -> pd.DataFrame:
    """Load the S1A conservation table, asserting required columns."""
    df = _read_parquet(Path(path), "conservation_table")
    _require_columns(df, CONSERVATION_TABLE_REQUIRED, "conservation_table")
    return df
