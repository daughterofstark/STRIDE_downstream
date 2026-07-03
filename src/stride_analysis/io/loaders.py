"""Raw loaders.

Thin readers that turn bytes into a DataFrame / validated model and do *no*
consistency checking (that lives in :mod:`stride_analysis.validation`).
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from ..models import MechanismFile
from ..models.errors import SchemaError


def load_correlations(path: str | Path) -> pd.DataFrame:
    """Load a Level-1 ``*_correlations_v5.csv`` replicate table."""
    p = Path(path)
    try:
        df = pd.read_csv(p)
    except Exception as exc:
        raise SchemaError(f"could not parse correlations CSV {p}: {exc}") from exc
    if df.empty:
        raise SchemaError(f"correlations CSV {p} is empty (no data rows)")
    return df


def load_profile(path: str | Path) -> pd.DataFrame:
    """Load a Level-2 ``*_profile.csv``."""
    p = Path(path)
    try:
        df = pd.read_csv(p)
    except Exception as exc:
        raise SchemaError(f"could not parse profile CSV {p}: {exc}") from exc
    if df.empty:
        raise SchemaError(f"profile CSV {p} is empty (no data rows)")
    return df


def load_mechanism(path: str | Path) -> MechanismFile:
    """Load and structurally validate a Level-2 ``*_mechanism.json``."""
    p = Path(path)
    try:
        raw = json.loads(p.read_text())
    except Exception as exc:
        raise SchemaError(f"could not parse mechanism JSON {p}: {exc}") from exc
    try:
        return MechanismFile.model_validate(raw)
    except Exception as exc:
        raise SchemaError(
            f"mechanism JSON {p} failed schema validation: {exc}"
        ) from exc
