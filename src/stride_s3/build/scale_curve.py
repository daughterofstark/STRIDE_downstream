"""Scale-of-resolution profiling: the ρ-vs-scale curve for every locus.

One row per (serotype, canon_label, scale_index): the locus's ρ at that scale,
the gain over the next-finer scale, the cumulative gain from residue scale, the
profile's own gated flag, and the scale's tier. The shape of this curve — how
fast ρ climbs as the region coarsens — is the mechanical basis for "how spread
out is the reproducible signal" (design §3.4, figure F7).

Pure: no IO, no mutation of inputs. ρ is read from the profile, never recomputed.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    RHO_DECIMALS,
    SCALE_CURVE_COLUMNS,
    SCALE_INDEX_TO_LEVEL,
)
from ._curves import scale_tier, step_gains
from ._frames import locus_scale_frame, rho_by_index_per_locus


def build_scale_curve(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the per-locus ρ-vs-scale curve table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).

    Returns
    -------
    DataFrame
        One row per (serotype, canon_label, scale_index), column-ordered per
        :data:`~stride_s3.models.schema.SCALE_CURVE_COLUMNS`.
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(SCALE_CURVE_COLUMNS))

    locus_scale = locus_scale_frame(stride_table)
    rho_curves = rho_by_index_per_locus(locus_scale)

    # per-locus static labels (chain/domain, gated index) from the profile
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
        gains = step_gains(rho_by_index)
        rho_residue = rho_by_index[0]
        gated_idx = gated_index.get(key, -1)
        for idx, rho in enumerate(rho_by_index):
            level = SCALE_INDEX_TO_LEVEL[idx]
            records.append(
                {
                    "serotype": serotype,
                    "canon_label": canon_label,
                    "chain": chain,
                    "domain": domain,
                    "scale_index": idx,
                    "scale_level": level,
                    "rho": _round(rho),
                    "rho_prev": _round(rho_by_index[idx - 1]) if idx > 0 else float("nan"),
                    "rho_step_gain": _round(gains[idx]),
                    "rho_cumulative_gain": _round(rho - rho_residue),
                    "is_gated_scale": bool(idx == gated_idx),
                    "tier": scale_tier(level),
                }
            )

    out = pd.DataFrame.from_records(records)[list(SCALE_CURVE_COLUMNS)]
    out = out.sort_values(
        ["serotype", "canon_label", "scale_index"]
    ).reset_index(drop=True)
    if out.duplicated(["serotype", "canon_label", "scale_index"]).any():
        raise ConsistencyError(
            "scale curve has duplicate (serotype, canon_label, scale_index) rows"
        )
    return out


def _round(value: float) -> float:
    """Round a ρ-derived value, preserving NaN."""
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
