"""Monotonicity / upward-closure (I2) audit.

One row per (serotype, canon_label): whether the locus's ρ is non-decreasing as
the scale coarsens (residue → complex), and — where it is not — how many steps
decrease, the largest drop, and where the first drop occurs (design §3.4). The
design measured 230/236 loci monotone in DENV1, so the audit is a diagnostic
that surfaces the handful of exceptions rather than a pass/fail gate.

Pure: no IO, no mutation of inputs. ρ is read from the profile, never recomputed.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    MONOTONICITY_AUDIT_COLUMNS,
    N_SCALES,
    RHO_DECIMALS,
)
from ._curves import audit_monotonicity
from ._frames import locus_scale_frame, rho_by_index_per_locus

_RESIDUE_INDEX = 0
_COMPLEX_INDEX = 6


def build_monotonicity_audit(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the per-locus upward-closure (monotonicity) audit table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).

    Returns
    -------
    DataFrame
        One row per (serotype, canon_label), column-ordered per
        :data:`~stride_s3.models.schema.MONOTONICITY_AUDIT_COLUMNS`.
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(MONOTONICITY_AUDIT_COLUMNS))

    locus_scale = locus_scale_frame(stride_table)
    rho_curves = rho_by_index_per_locus(locus_scale)

    labels: dict[tuple[str, str], tuple[str, str]] = {}
    for row in locus_scale.itertuples(index=False):
        key = (str(row.serotype), str(row.canon_label))
        labels.setdefault(key, (str(row.h_chain), str(row.h_domain)))

    records = []
    for key, rho_by_index in rho_curves.items():
        serotype, canon_label = key
        chain, domain = labels[key]
        is_monotone, n_violations, max_decrease, first_violation = (
            audit_monotonicity(rho_by_index)
        )
        records.append(
            {
                "serotype": serotype,
                "canon_label": canon_label,
                "chain": chain,
                "domain": domain,
                "n_scales": N_SCALES,
                "is_monotone": bool(is_monotone),
                "n_violations": int(n_violations),
                "max_decrease": _round(max_decrease),
                "first_violation_scale_index": int(first_violation),
                "rho_residue": _round(rho_by_index[_RESIDUE_INDEX]),
                "rho_complex": _round(rho_by_index[_COMPLEX_INDEX]),
            }
        )

    out = pd.DataFrame.from_records(records)[
        list(MONOTONICITY_AUDIT_COLUMNS)
    ]
    out = out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
    if out.duplicated(["serotype", "canon_label"]).any():
        raise ConsistencyError(
            "monotonicity audit has duplicate (serotype, canon_label) rows"
        )
    return out


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
