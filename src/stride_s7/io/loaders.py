"""Loaders for the S7 inputs.

S7 reads the S2–S6 reduction tables — never the raw STRIDE files and never the
S0/S1 canonical layer. Each input is read from its producing stage's output
directory, and the columns S7 depends on are asserted present. :func:`file_digest`
computes the SHA-256 of an input for the provenance header the design requires
(§5.4).

The :data:`STAGE_INPUTS` registry is the authoritative map of *which stage
directory* each input table is read from, and :data:`INPUT_REQUIRED_COLUMNS` maps
each table to the columns S7 needs; both keep the loader declarative.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from ..models.errors import InputError
from ..models.schema import (
    DOMAIN_EFFECT_SUMMARY_REQUIRED,
    DOMAIN_REPRODUCIBILITY_REQUIRED,
    DOMAIN_SEROTYPE_MATRIX_REQUIRED,
    IN_DOMAIN_EFFECT_SUMMARY,
    IN_DOMAIN_REPRODUCIBILITY,
    IN_DOMAIN_SEROTYPE_MATRIX,
    IN_POSITION_CONSERVATION,
    IN_REPLICATE_BLOCKED_ANALYSES,
    IN_REPLICATE_REGIME,
    IN_RESIDUE_LANDSCAPE,
    IN_RESOLUTION_CENSUS,
    IN_SCALE_CURVE,
    IN_SEROTYPE_SUMMARY,
    IN_SIGNIFICANCE_SCREEN,
    IN_VARIANCE_BUDGET,
    POSITION_CONSERVATION_REQUIRED,
    REPLICATE_BLOCKED_ANALYSES_REQUIRED,
    REPLICATE_REGIME_REQUIRED,
    RESIDUE_LANDSCAPE_REQUIRED,
    RESOLUTION_CENSUS_REQUIRED,
    SCALE_CURVE_REQUIRED,
    SEROTYPE_SUMMARY_REQUIRED,
    SIGNIFICANCE_SCREEN_REQUIRED,
    VARIANCE_BUDGET_REQUIRED,
)

#: input filename → the producing stage key ("s2".."s6")
STAGE_INPUTS: dict[str, str] = {
    IN_RESIDUE_LANDSCAPE: "s2",
    IN_RESOLUTION_CENSUS: "s2",
    IN_SEROTYPE_SUMMARY: "s2",
    IN_DOMAIN_REPRODUCIBILITY: "s2",
    IN_SCALE_CURVE: "s3",
    IN_SIGNIFICANCE_SCREEN: "s4",
    IN_VARIANCE_BUDGET: "s4",
    IN_DOMAIN_EFFECT_SUMMARY: "s4",
    IN_DOMAIN_SEROTYPE_MATRIX: "s5",
    IN_POSITION_CONSERVATION: "s5",
    IN_REPLICATE_REGIME: "s6",
    IN_REPLICATE_BLOCKED_ANALYSES: "s6",
}

#: input filename → the columns S7 depends on
INPUT_REQUIRED_COLUMNS: dict[str, tuple[str, ...]] = {
    IN_RESIDUE_LANDSCAPE: RESIDUE_LANDSCAPE_REQUIRED,
    IN_RESOLUTION_CENSUS: RESOLUTION_CENSUS_REQUIRED,
    IN_SEROTYPE_SUMMARY: SEROTYPE_SUMMARY_REQUIRED,
    IN_DOMAIN_REPRODUCIBILITY: DOMAIN_REPRODUCIBILITY_REQUIRED,
    IN_SCALE_CURVE: SCALE_CURVE_REQUIRED,
    IN_SIGNIFICANCE_SCREEN: SIGNIFICANCE_SCREEN_REQUIRED,
    IN_VARIANCE_BUDGET: VARIANCE_BUDGET_REQUIRED,
    IN_DOMAIN_EFFECT_SUMMARY: DOMAIN_EFFECT_SUMMARY_REQUIRED,
    IN_DOMAIN_SEROTYPE_MATRIX: DOMAIN_SEROTYPE_MATRIX_REQUIRED,
    IN_POSITION_CONSERVATION: POSITION_CONSERVATION_REQUIRED,
    IN_REPLICATE_REGIME: REPLICATE_REGIME_REQUIRED,
    IN_REPLICATE_BLOCKED_ANALYSES: REPLICATE_BLOCKED_ANALYSES_REQUIRED,
}


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


def resolve_input_paths(stage_dirs: dict[str, Path]) -> dict[str, Path]:
    """Map each required input filename to its path under the producing stage dir.

    Parameters
    ----------
    stage_dirs
        Mapping of stage key (``"s2"``..``"s6"``) → that stage's output directory.
    """
    paths: dict[str, Path] = {}
    for filename, stage in STAGE_INPUTS.items():
        base = stage_dirs[stage]
        paths[filename] = Path(base) / filename
    return paths


def load_inputs(paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    """Load every required input table, asserting the columns S7 depends on.

    A missing or malformed input is a hard failure (:class:`InputError`).
    """
    frames: dict[str, pd.DataFrame] = {}
    for filename, required in INPUT_REQUIRED_COLUMNS.items():
        df = _read_parquet(paths[filename], filename)
        _require_columns(df, required, filename)
        frames[filename] = df
    return frames


def file_digest(path: str | Path) -> str:
    """Return the SHA-256 hex digest of a file, for the provenance header.

    Returns the empty string if the path does not exist.
    """
    p = Path(path)
    if not p.is_file():
        return ""
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
