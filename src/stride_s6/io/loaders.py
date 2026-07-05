"""Loaders for the S6 inputs.

S6 consumes the S1A ``replicate_inventory`` (the replicate-structure index —
**required**; S1A always emits it) and, when it exists, the S0 ``replicate_table``
(the Level-1 per-run observations — **optional**). The replicate table is only
written by S0 when the per-run correlation CSVs were supplied; its absence is the
design's anticipated "per-run θ unavailable" state (§4.1), not an error, so
:func:`load_replicate_table` returns ``None`` when the file is missing rather than
raising.

These loaders read the parquet, assert the columns S6 depends on are present, and
return the DataFrames unchanged otherwise. :func:`file_digest` computes the
SHA-256 of an input file for the provenance header the design requires (§5.4).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from ..models.errors import InputError
from ..models.schema import (
    REPLICATE_INVENTORY_REQUIRED,
    REPLICATE_TABLE_REQUIRED,
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


def load_replicate_inventory(path: str | Path) -> pd.DataFrame:
    """Load the S1A replicate inventory, asserting the columns S6 depends on.

    This input is required; a missing or malformed inventory is an error.
    """
    df = _read_parquet(Path(path), "replicate_inventory")
    _require_columns(df, REPLICATE_INVENTORY_REQUIRED, "replicate_inventory")
    return df


def load_replicate_table(path: str | Path) -> pd.DataFrame | None:
    """Load the S0 replicate table if it exists, else return ``None``.

    The replicate table carries the Level-1 per-run observations (per-run θ). It is
    optional: when the per-run correlation CSVs were not supplied to S0 the file is
    absent, which is the design's anticipated blocked state (§4.1) rather than an
    error. When the file *is* present it must carry the columns S6 depends on.
    """
    p = Path(path)
    if not p.is_file():
        return None
    df = _read_parquet(p, "replicate_table")
    _require_columns(df, REPLICATE_TABLE_REQUIRED, "replicate_table")
    return df


def file_digest(path: str | Path) -> str:
    """Return the SHA-256 hex digest of a file, for the provenance header.

    Returns the empty string if the path does not exist, so provenance stamping
    never fails a run (e.g. when the optional replicate table is absent).
    """
    p = Path(path)
    if not p.is_file():
        return ""
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
