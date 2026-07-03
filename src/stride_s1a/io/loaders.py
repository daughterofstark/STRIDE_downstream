"""Loaders for the S0 canonical tables.

S1A consumes **only** ``stride_table.parquet`` and ``replicate_table.parquet``
(or paths pointing at equivalent canonical outputs). It never re-reads the raw
STRIDE CSVs/JSONs. These loaders read the parquet, assert the columns S1A
depends on are present, and return DataFrames unchanged otherwise.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..models.errors import InputError
from ..models.schema import (
    REPLICATE_REQUIRED_COLUMNS,
    STRIDE_REQUIRED_COLUMNS,
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


def load_stride_table(path: str | Path) -> pd.DataFrame:
    """Load the S0 STRIDE table, asserting the columns S1A needs are present."""
    df = _read_parquet(Path(path), "stride_table")
    _require_columns(df, STRIDE_REQUIRED_COLUMNS, "stride_table")
    return df


def load_replicate_table(path: str | Path) -> pd.DataFrame:
    """Load the S0 replicate table, asserting the columns S1A needs are present.

    S0 only writes ``replicate_table.parquet`` when it is non-empty (a
    summaries-only S0 run omits it). A **missing** file is therefore treated as
    a valid "no replicates" state and yields an empty DataFrame; an existing but
    unreadable file is still an error. An empty-but-well-formed table is returned
    as-is for the caller to handle.
    """
    p = Path(path)
    if not p.is_file():
        return pd.DataFrame(columns=list(REPLICATE_REQUIRED_COLUMNS))
    try:
        df = pd.read_parquet(p)
    except Exception as exc:
        raise InputError(f"could not read replicate_table at {p}: {exc}") from exc
    if not df.empty:
        _require_columns(df, REPLICATE_REQUIRED_COLUMNS, "replicate_table")
    return df
