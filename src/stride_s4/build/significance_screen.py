"""CI-based significance screen over gated mechanisms, with FDR control.

One row per gated mechanism (keyed by the residue locus carrying the gated row):
its signed effect, the CI-exclusion test, a two-sided Wald p-value, and the
Benjamini–Hochberg adjusted p-value computed **within each serotype** over its
signed mechanisms. This realises the design's CI-based significance screen
(§3.5) with FDR control across the positions of a serotype (§5.2); the serotype
is the unit of replication, so it is the multiple-testing family.

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DIRECTION_MIXED,
    GATE_ALPHA,
    P_DECIMALS,
    RHO_DECIMALS,
    SIGNIFICANCE_SCREEN_COLUMNS,
)
from ._frames import gated_slice, scale_tier
from ._stats import (
    benjamini_hochberg,
    ci_excludes_zero,
    two_sided_wald_p,
    wald_z,
)


def build_significance_screen(stride_table: pd.DataFrame) -> pd.DataFrame:
    """Build the per-mechanism CI significance screen with BH-FDR.

    Parameters
    ----------
    stride_table
        The S0 STRIDE table (all scales).

    Returns
    -------
    DataFrame
        One row per gated mechanism, column-ordered per
        :data:`~stride_s4.models.schema.SIGNIFICANCE_SCREEN_COLUMNS`.
    """
    if stride_table.empty:
        return pd.DataFrame(columns=list(SIGNIFICANCE_SCREEN_COLUMNS))

    gated = gated_slice(stride_table)
    if gated.empty:
        return pd.DataFrame(columns=list(SIGNIFICANCE_SCREEN_COLUMNS))

    records = []
    for row in gated.itertuples(index=False):
        direction = (
            str(row.mech_direction)
            if pd.notna(row.mech_direction)
            else DIRECTION_MIXED
        )
        is_signed = direction != DIRECTION_MIXED
        beta_signed = _f(row.mech_beta_signed)
        beta_se = _f(row.mech_beta_se)
        ci_lower = _f(row.mech_beta_ci_lower)
        ci_upper = _f(row.mech_beta_ci_upper)

        excludes = is_signed and ci_excludes_zero(ci_lower, ci_upper)
        z = wald_z(beta_signed, beta_se) if is_signed else float("nan")
        p = two_sided_wald_p(beta_signed, beta_se) if is_signed else float("nan")

        records.append(
            {
                "serotype": str(row.serotype),
                "canon_label": str(row.canon_label),
                "chain": str(row.h_chain),
                "domain": str(row.h_domain),
                "gated_scale_level": str(row.scale_level),
                "tier": scale_tier(str(row.scale_level)),
                "direction": direction,
                "is_signed": bool(is_signed),
                "beta_signed": _round(beta_signed),
                "beta_se": _round(beta_se),
                "beta_ci_lower": _round(ci_lower),
                "beta_ci_upper": _round(ci_upper),
                "ci_excludes_zero": bool(excludes),
                "z_score": _round(z),
                "p_value": _round_p(p),
                # p_value_bh + significant_fdr filled in per serotype below
            }
        )

    out = pd.DataFrame.from_records(records)

    # Benjamini–Hochberg within each serotype over its signed mechanisms.
    out["p_value_bh"] = float("nan")
    for _serotype, grp in out.groupby("serotype", sort=False):
        p_bh = benjamini_hochberg(grp["p_value"].tolist())
        out.loc[grp.index, "p_value_bh"] = [_round_p(v) for v in p_bh]

    out["significant_raw"] = out["ci_excludes_zero"].astype(bool)
    out["significant_fdr"] = out["is_signed"] & out["p_value_bh"].apply(
        lambda v: bool(v == v and v <= GATE_ALPHA)  # v == v filters NaN
    )

    out = out[list(SIGNIFICANCE_SCREEN_COLUMNS)]
    out = out.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
    if out.duplicated(["serotype", "canon_label"]).any():
        raise ConsistencyError(
            "significance screen has duplicate (serotype, canon_label) rows"
        )
    return out


def _f(value: object) -> float:
    """Coerce to float, mapping missing/None to NaN."""
    if value is None:
        return float("nan")
    try:
        f = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")
    return f


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)


def _round_p(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), P_DECIMALS)
