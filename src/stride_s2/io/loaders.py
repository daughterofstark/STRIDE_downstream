"""Loaders for the S2 inputs.

S2 consumes the S0 STRIDE table (``stride_table.parquet``) and the two S1B
annotation tables (``residue_annotation.parquet``, ``domain_annotation.parquet``).
It never re-reads the raw STRIDE CSV/JSON files or MD trajectories. These loaders
read the parquet, assert the columns S2 depends on are present, and return the
DataFrames unchanged otherwise.

:func:`file_digest` computes the SHA-256 of an input file for the provenance
header the design requires on every output (§5.4).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from ..models.errors import InputError
from ..models.schema import (
    DOMAIN_ANNOTATION_REQUIRED,
    RESIDUE_ANNOTATION_REQUIRED,
    STRIDE_TABLE_REQUIRED,
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
    """Load the S0 STRIDE table, asserting the columns S2 depends on."""
    df = _read_parquet(Path(path), "stride_table")
    _require_columns(df, STRIDE_TABLE_REQUIRED, "stride_table")
    return df


def load_residue_annotation(path: str | Path) -> pd.DataFrame:
    """Load the S1B residue annotation, asserting required columns."""
    df = _read_parquet(Path(path), "residue_annotation")
    _require_columns(df, RESIDUE_ANNOTATION_REQUIRED, "residue_annotation")
    return df


def load_domain_annotation(path: str | Path) -> pd.DataFrame:
    """Load the S1B domain annotation, asserting required columns."""
    df = _read_parquet(Path(path), "domain_annotation")
    _require_columns(df, DOMAIN_ANNOTATION_REQUIRED, "domain_annotation")
    return df


def file_digest(path: str | Path) -> str:
    """Return the SHA-256 hex digest of a file, for the provenance header.

    Returns the empty string if the path does not exist (a partial/optional
    input), so provenance stamping never fails a run.
    """
    p = Path(path)
    if not p.is_file():
        return ""
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
