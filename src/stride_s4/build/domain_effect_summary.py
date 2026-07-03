"""β_se-weighted effect summary per domain (Tier A — licensed).

One row per (serotype, chain, domain) carrying at least one gated mechanism: the
count composition of its mechanisms, the fraction of signed mechanisms whose CI
excludes 0, the number significant after FDR, and the 1/β_se²-weighted mean
signed effect. This realises the design's β_se-weighted summaries and the
per-domain CI-exclusion fraction (§3.5).

Gated mechanisms are grouped into their domain by the ``chain``/``domain`` labels
of the gated row. The builder consumes the significance screen so its CI and FDR
decisions stay consistent with the per-mechanism table.

Pure: no IO, no mutation of inputs.
"""
from __future__ import annotations

import pandas as pd

from ..models.errors import ConsistencyError
from ..models.schema import (
    DOMAIN_EFFECT_SUMMARY_COLUMNS,
    RHO_DECIMALS,
    TIER_LICENSED,
)
from ._stats import weighted_mean_se


def build_domain_effect_summary(
    significance_screen: pd.DataFrame,
) -> pd.DataFrame:
    """Build the per-domain β_se-weighted effect summary.

    Parameters
    ----------
    significance_screen
        The S4 significance-screen table (one row per gated mechanism), as built
        by :func:`~stride_s4.build.significance_screen.build_significance_screen`.

    Returns
    -------
    DataFrame
        One row per (serotype, chain, domain) with ≥1 mechanism, column-ordered
        per :data:`~stride_s4.models.schema.DOMAIN_EFFECT_SUMMARY_COLUMNS`.
    """
    if significance_screen.empty:
        return pd.DataFrame(columns=list(DOMAIN_EFFECT_SUMMARY_COLUMNS))

    records = []
    group_cols = ["serotype", "chain", "domain"]
    for keys, grp in significance_screen.groupby(group_cols, sort=True):
        serotype, chain, domain = (str(k) for k in keys)
        signed = grp[grp["is_signed"].astype(bool)]
        n_mechanisms = int(len(grp))
        n_signed = int(len(signed))
        n_mixed = n_mechanisms - n_signed
        n_ci = int(signed["ci_excludes_zero"].astype(bool).sum())
        n_fdr = int(signed["significant_fdr"].astype(bool).sum())
        frac_ci = (n_ci / n_signed) if n_signed > 0 else float("nan")

        w_mean, w_se, u_mean = weighted_mean_se(
            signed["beta_signed"].tolist(),
            signed["beta_se"].tolist(),
        )

        records.append(
            {
                "serotype": serotype,
                "chain": chain,
                "domain": domain,
                "n_mechanisms": n_mechanisms,
                "n_signed": n_signed,
                "n_mixed": n_mixed,
                "n_ci_excludes_zero": n_ci,
                "frac_ci_excludes_zero": _round(frac_ci),
                "n_significant_fdr": n_fdr,
                "beta_weighted_mean": _round(w_mean),
                "beta_weighted_se": _round(w_se),
                "beta_unweighted_mean": _round(u_mean),
                "tier": TIER_LICENSED,
            }
        )

    out = pd.DataFrame.from_records(records)[
        list(DOMAIN_EFFECT_SUMMARY_COLUMNS)
    ]
    out = out.sort_values(["serotype", "chain", "domain"]).reset_index(drop=True)
    if out.duplicated(["serotype", "chain", "domain"]).any():
        raise ConsistencyError(
            "domain effect summary has duplicate (serotype, chain, domain) rows"
        )
    return out


def _round(value: float) -> float:
    if value != value:  # NaN
        return float("nan")
    return round(float(value), RHO_DECIMALS)
