"""Domain-vs-residue reproducibility gap: Δρ = ρ(domain) − ρ(residue).

One row per (serotype, canon_label): the locus's ρ at residue and domain scale,
the gap between them, and the ``is_distributed`` flag — a locus that fails the
provisional gate at residue scale but clears it at domain scale, i.e. an effect
that is reproducible only once aggregated (design §3.4). A large positive gap
marks a distributed effect; a near-zero gap marks a locally-resolved one.

Pure: no IO, no mutation of inputs. ρ is read from the profile, never recomputed.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    PROVISIONAL_RHO_STAR,
    RESOLUTION_GAP_COLUMNS,
    RHO_DECIMALS,
    SCALE_INDEX_TO_LEVEL,
    TIER_LICENSED,
)
from ._curves import is_distributed_effect
from ._frames import locus_scale_frame, rho_by_index_per_locus

_RESIDUE_INDEX = 0
_DOMAIN_INDEX = 3


def build_resolution_gap(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the per-locus domain-vs-residue reproducibility gap table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).

    Returns
    -------
    DataFrame
        One row per (serotype, canon_label), column-ordered per
        :data:`~stride_s3.models.schema.RESOLUTION_GAP_COLUMNS`.
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(RESOLUTION_GAP_COLUMNS))

    locus_scale = locus_scale_frame(stride_table)
    rho_curves = rho_by_index_per_locus(locus_scale)

    labels: dict[tuple[str, str], tuple[str, str]] = {}
    gated_index: dict[tuple[str, str], int] = {}
    for row in locus_scale.itertuples(index=False):
        key = (str(row.serotype), str(row.canon_label))
        labels.setdefault(key, (str(row.h_chain), str(row.h_domain)))
        if bool(row.is_gated_scale):
            gated_index[key] = int(row.scale_index)

    records = []
    for key, rho_by_index in rho_curves.items():
        serotype, canon_label = key
        chain, domain = labels[key]
        rho_residue = rho_by_index[_RESIDUE_INDEX]
        rho_domain = rho_by_index[_DOMAIN_INDEX]
        finite = [r for r in rho_by_index if r == r]  # drop NaN
        rho_min = min(finite) if finite else float("nan")
        rho_max = max(finite) if finite else float("nan")
        gated_idx = gated_index.get(key, -1)
        rho_at_gated = (
            rho_by_index[gated_idx]
            if 0 <= gated_idx < len(rho_by_index)
            else float("nan")
        )
        records.append(
            {
                "serotype": serotype,
                "canon_label": canon_label,
                "chain": chain,
                "domain": domain,
                "rho_residue": _round(rho_residue),
                "rho_domain": _round(rho_domain),
                "delta_rho_domain_residue": _round(rho_domain - rho_residue),
                "rho_min": _round(rho_min),
                "rho_max": _round(rho_max),
                "gated_scale_level": (
                    SCALE_INDEX_TO_LEVEL.get(gated_idx, "")
                    if gated_idx >= 0
                    else ""
                ),
                "gated_scale_index": gated_idx,
                "rho_at_gated": _round(rho_at_gated),
                "is_distributed": is_distributed_effect(
                    rho_residue, rho_domain, PROVISIONAL_RHO_STAR
                ),
                "domain_tier": TIER_LICENSED,
            }
        )

    out = pd.DataFrame.from_records(records)[list(RESOLUTION_GAP_COLUMNS)]
    out = out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
    if out.duplicated(["serotype", "canon_label"]).any():
        raise ConsistencyError(
            "resolution gap has duplicate (serotype, canon_label) rows"
        )
    return out


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
