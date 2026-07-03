"""Internal frame-extraction helpers shared by the S2 builders.

These pull the reusable views S2 needs out of the S0 STRIDE table — the
per-(serotype, locus) ρ-by-scale map used for re-gating, the residue-scale
slice, the domain-scale slice, and the gated-mechanism slice — with the input
never mutated. Pure functions; no IO.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DOMAIN_SCALE_LEVEL,
    RESIDUE_SCALE_LEVEL,
    SCALE_LEVEL_TO_INDEX,
)


def rho_by_scale_map(
    stride_table: pd.DataFrame,
) -> dict[tuple[str, str], dict[str, float]]:
    """Map ``(serotype, canon_label) -> {scale_level: rho}`` for every locus.

    This is the input to re-gating: each locus's ρ at all seven scales. Only
    recognised scale levels are kept; ρ is read straight from the profile and
    never recomputed.
    """
    out: dict[tuple[str, str], dict[str, float]] = {}
    for serotype, canon_label, scale_level, rho in zip(
        stride_table["serotype"],
        stride_table["canon_label"],
        stride_table["scale_level"],
        stride_table["rho"],
        strict=True,
    ):
        if scale_level not in SCALE_LEVEL_TO_INDEX:
            continue
        key = (str(serotype), str(canon_label))
        out.setdefault(key, {})[str(scale_level)] = (
            float(rho) if pd.notna(rho) else float("nan")
        )
    return out


def residue_slice(stride_table: pd.DataFrame) -> pd.DataFrame:
    """The residue-scale rows (one per locus), sorted deterministically."""
    res = stride_table[
        stride_table["scale_level"] == RESIDUE_SCALE_LEVEL
    ].copy()
    dup = res.duplicated(["serotype", "canon_label"], keep=False)
    if dup.any():
        examples = (
            res.loc[dup, ["serotype", "canon_label"]]
            .drop_duplicates()
            .head(3)
            .to_dict("records")
        )
        raise ConsistencyError(
            "STRIDE table has duplicate residue-scale rows for "
            f"(serotype, canon_label); examples: {examples}"
        )
    return res.sort_values(["serotype", "canon_label"]).reset_index(drop=True)


def domain_slice(stride_table: pd.DataFrame) -> pd.DataFrame:
    """The domain-scale rows, sorted deterministically.

    At domain scale a region's ρ/β/coherence are region-level constants shared
    by every member locus, so multiple loci may carry identical domain rows.
    """
    dom = stride_table[
        stride_table["scale_level"] == DOMAIN_SCALE_LEVEL
    ].copy()
    return dom.sort_values(
        ["serotype", "h_chain", "h_domain", "canon_label"]
    ).reset_index(drop=True)


def gated_slice(stride_table: pd.DataFrame) -> pd.DataFrame:
    """The gated rows (one per emitted mechanism), sorted deterministically."""
    gated = stride_table[stride_table["is_gated_scale"].astype(bool)].copy()
    return gated.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
