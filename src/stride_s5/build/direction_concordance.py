"""Direction concordance for shared signed positions (Tier B — exploratory).

One row per shared ``canon_label`` that is signed *and* reproducible in at least
``MIN_SEROTYPES_FOR_CONCORDANCE`` serotypes: do those serotypes agree on the sign
of the effect (increase vs decrease)? This is the design's direction-concordance
product for shared, signed positions (§3.3).

Serotype is the unit of replication: the per-serotype signed direction is
aggregated first (:func:`~stride_s5.build._frames.position_frame`), then the
directions are tallied across serotypes — never across residues (§5.2). The
concordance is a descriptive tally (agree / majority / conflict), not a p-value,
in keeping with the n = 4 guardrail.

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DIRECTION_CONCORDANCE_COLUMNS,
    DIRECTION_DECREASE,
    DIRECTION_INCREASE,
    MIN_SEROTYPES_FOR_CONCORDANCE,
    PROVISIONAL_RHO_STAR,
    RHO_DECIMALS,
    TIER_EXPLORATORY,
)
from ._classify import (
    concordance_class,
    is_catalytic_triad,
    majority_direction,
)
from ._frames import position_frame


def build_direction_concordance(
    stride_table: pd.DataFrame,
    rho_star: float = PROVISIONAL_RHO_STAR,
) -> pd.DataFrame:
    """Build the per-position direction-concordance table.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).
    rho_star
        The provisional gate threshold (default :data:`PROVISIONAL_RHO_STAR`).

    Returns
    -------
    DataFrame
        One row per shared signed ``canon_label`` (signed-and-reproducible in
        ≥ ``MIN_SEROTYPES_FOR_CONCORDANCE`` serotypes), column-ordered per
        :data:`~stride_s5.models.schema.DIRECTION_CONCORDANCE_COLUMNS`.
    """
    pf = position_frame(stride_table, rho_star)
    if pf.empty:
        return pd.DataFrame(columns=list(DIRECTION_CONCORDANCE_COLUMNS))

    signed = pf[pf["is_signed"].astype(bool)]
    if signed.empty:
        return pd.DataFrame(columns=list(DIRECTION_CONCORDANCE_COLUMNS))

    records = []
    for canon, grp in signed.groupby("canon_label", sort=True):
        n_signed = int(len(grp))
        if n_signed < MIN_SEROTYPES_FOR_CONCORDANCE:
            continue
        n_inc = int((grp["direction"] == DIRECTION_INCREASE).sum())
        n_dec = int((grp["direction"] == DIRECTION_DECREASE).sum())
        frac_majority = max(n_inc, n_dec) / n_signed if n_signed > 0 else float("nan")
        first = grp.iloc[0]
        records.append(
            {
                "canon_label": str(canon),
                "chain": str(first["chain"]),
                "domain": str(first["domain"]),
                "n_serotypes_signed": n_signed,
                "n_increase": n_inc,
                "n_decrease": n_dec,
                "majority_direction": majority_direction(n_inc, n_dec),
                "frac_majority": _round(frac_majority),
                "concordance_class": concordance_class(n_inc, n_dec),
                "is_catalytic_triad": is_catalytic_triad(str(canon)),
                "rho_star": rho_star,
                "is_provisional_rho_star": True,
                "tier": TIER_EXPLORATORY,
            }
        )

    if not records:
        return pd.DataFrame(columns=list(DIRECTION_CONCORDANCE_COLUMNS))

    out = pd.DataFrame.from_records(records)[
        list(DIRECTION_CONCORDANCE_COLUMNS)
    ]
    out = out.sort_values(["canon_label"]).reset_index(drop=True)
    if out.duplicated(["canon_label"]).any():
        raise ConsistencyError(
            "direction concordance has duplicate canon_label rows"
        )
    return out


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
