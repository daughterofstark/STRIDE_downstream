"""Level-1 (replicate) schema validation.

Validates a ``*_correlations_v5.csv`` replicate table. Because STRIDE's engine
appends columns across milestones, the schema is validated as a *required core*
(strict) plus *known-optional* columns (type-checked when present, tolerated
when absent). Unknown extra columns are allowed but reported — a replicate table
is an observation record we ingest faithfully, not a frozen contract like the
Level-2 summaries.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..models.errors import SchemaError
from ..models.schema import (
    CORRELATIONS_FLOAT_COLUMNS,
    CORRELATIONS_INT_COLUMNS,
    CORRELATIONS_KNOWN_OPTIONAL_COLUMNS,
    CORRELATIONS_REQUIRED_COLUMNS,
    CORRELATIONS_STR_COLUMNS,
)


def validate_correlations_schema(
    df: pd.DataFrame, serotype: str, replicate: str
) -> list[str]:
    """Validate a replicate table; return any unknown extra column names.

    Raises :class:`SchemaError` on: missing required columns, wrong dtypes on
    required/known columns, blank required identifiers, non-finite effect
    values, or duplicate ``canon_label`` rows (a replicate must have one row per
    residue).
    """
    where = f"[{serotype}/{replicate} correlations]"

    # -- required columns ----------------------------------------------------
    missing = [c for c in CORRELATIONS_REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise SchemaError(f"{where} missing required columns: {missing}")

    # -- required identifiers non-null/blank ---------------------------------
    for col in ("name", "label"):
        if df[col].isna().any():
            raise SchemaError(f"{where} column {col!r} has missing value(s)")
        if (df[col].astype(str).str.strip() == "").any():
            raise SchemaError(f"{where} column {col!r} has blank value(s)")

    # -- dtypes on required + known-optional columns -------------------------
    known = set(CORRELATIONS_REQUIRED_COLUMNS) | set(
        CORRELATIONS_KNOWN_OPTIONAL_COLUMNS
    )
    for col in df.columns:
        if col not in known:
            continue
        if col in CORRELATIONS_FLOAT_COLUMNS:
            coerced = pd.to_numeric(df[col], errors="coerce")
            bad = coerced.isna() & df[col].notna()
            if bad.any():
                idx = list(df.index[bad][:5])
                raise SchemaError(
                    f"{where} column {col!r} has non-numeric value(s) at rows "
                    f"{idx}"
                )
        elif col in CORRELATIONS_INT_COLUMNS:
            coerced = pd.to_numeric(df[col], errors="coerce")
            if coerced.isna().any() or (coerced % 1 != 0).any():
                raise SchemaError(
                    f"{where} column {col!r} must be integer-valued"
                )
        elif col in CORRELATIONS_STR_COLUMNS:
            # allow NaN in optional string columns, but not blanks in present ones
            non_null = df[col].dropna()
            if (non_null.astype(str).str.strip() == "").any():
                raise SchemaError(
                    f"{where} column {col!r} has blank value(s)"
                )

    # -- effect field must be finite -----------------------------------------
    for col in ("r", "abs_r"):
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(vals).all():
            idx = list(df.index[~np.isfinite(vals)][:5])
            raise SchemaError(
                f"{where} effect column {col!r} has non-finite value(s) at "
                f"rows {idx}"
            )

    # -- one row per residue (canon_label unique within a replicate) ---------
    dup = df.duplicated(subset=["label"], keep=False)
    if dup.any():
        examples = df.loc[dup, "label"].drop_duplicates().head(5).tolist()
        raise SchemaError(
            f"{where} duplicate residue label(s) within one replicate: "
            f"{examples}"
        )

    return [c for c in df.columns if c not in known]
