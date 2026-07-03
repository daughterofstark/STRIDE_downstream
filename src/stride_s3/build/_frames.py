"""Internal frame-extraction helpers shared by the S3 builders.

These pull the reusable per-locus views S3 needs out of the S0 STRIDE table — the
dense ρ-by-scale array used for the curve and the monotonicity audit, and the
per-scale slices used for the gap and chain contrast — with the input never
mutated. Pure functions; no IO.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DOMAIN_SCALE_LEVEL,
    N_SCALES,
    RESIDUE_SCALE_LEVEL,
    SCALE_LEVEL_TO_INDEX,
)


def locus_scale_frame(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Return the STRIDE rows S3 uses, sorted by locus then scale index.

    Only rows whose ``scale_level`` is a recognised scale are kept. The result is
    deterministically ordered so every downstream reduction is reproducible.
    """
    keep = stride_table["scale_level"].isin(SCALE_LEVEL_TO_INDEX)
    frame = stride_table[keep].copy()
    return frame.sort_values(
        ["serotype", "canon_label", "scale_index"]
    ).reset_index(drop=True)


def rho_by_index_per_locus(
    locus_scale: pd.DataFrame,
) -> dict[tuple[str, str], list[float]]:
    """Map ``(serotype, canon_label) -> [ρ at scale index 0 .. 6]``.

    Each locus must carry exactly one row per scale index (a dense curve). A
    missing or duplicated scale raises :class:`ConsistencyError`.
    """
    out: dict[tuple[str, str], list[float]] = {}
    for (serotype, canon_label), grp in locus_scale.groupby(
        ["serotype", "canon_label"], sort=True
    ):
        by_index: dict[int, float] = {}
        for scale_index, rho in zip(
            grp["scale_index"], grp["rho"], strict=True
        ):
            idx = int(scale_index)
            if idx in by_index:
                raise ConsistencyError(
                    f"locus ({serotype}, {canon_label}) has duplicate rows for "
                    f"scale_index {idx}"
                )
            by_index[idx] = float(rho) if pd.notna(rho) else float("nan")
        missing = [i for i in range(N_SCALES) if i not in by_index]
        if missing:
            raise ConsistencyError(
                f"locus ({serotype}, {canon_label}) is missing scale "
                f"index/indices {missing}; expected one row per scale 0..{N_SCALES - 1}"
            )
        out[(str(serotype), str(canon_label))] = [
            by_index[i] for i in range(N_SCALES)
        ]
    return out


def scale_slice(locus_scale: pd.DataFrame, scale_level: str) -> pd.DataFrame:
    """The rows at one scale level, indexed by (serotype, canon_label)."""
    return locus_scale[locus_scale["scale_level"] == scale_level].copy()


def residue_slice(locus_scale: pd.DataFrame) -> pd.DataFrame:
    """The residue-scale rows (one per locus)."""
    return scale_slice(locus_scale, RESIDUE_SCALE_LEVEL)


def domain_slice(locus_scale: pd.DataFrame) -> pd.DataFrame:
    """The domain-scale rows (one per locus)."""
    return scale_slice(locus_scale, DOMAIN_SCALE_LEVEL)


def gated_slice(locus_scale: pd.DataFrame) -> pd.DataFrame:
    """The gated rows (one per emitted mechanism), sorted deterministically."""
    gated = locus_scale[locus_scale["is_gated_scale"].astype(bool)].copy()
    return gated.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
