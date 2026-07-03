"""Structural validation of the S4 uncertainty layer.

Structural / arithmetic checks only — no statistical assertions and no biological
claims (the gate is uncalibrated, §0.1). Each check appends a
:class:`~stride_s4.models.ValidationCheck` to the report and raises
:class:`~stride_s4.models.errors.ConsistencyError` on failure, so later stages
can trust the uncertainty layer.

Checks:

- every output table's key is unique;
- the variance tables are *arithmetically consistent* — where the unreproduced
  variance is positive, ``frac_tau2 + frac_sigma2 == 1``, both fractions lie in
  ``[0, 1]``, and the regime label matches the recomputed τ² fraction;
- the significance screen is *self-consistent* — ``is_signed`` iff the direction
  is not mixed, ``ci_excludes_zero``/``significant_fdr`` only for signed rows, a
  BH-adjusted p-value is present exactly for signed rows, and the BH values are
  monotone in the raw-p order within each serotype;
- the domain effect summary *partitions* each domain's mechanisms
  (``n_signed + n_mixed == n_mechanisms``; ``n_ci_excludes_zero`` and
  ``n_significant_fdr`` never exceed ``n_signed``; the CI fraction lies in
  ``[0, 1]``).
"""
from __future__ import annotations

import math

import pandas as pd

from ..build._stats import classify_variance_regime
from ..models import S4Report
from ..models.errors import ConsistencyError
from ..models.schema import (
    DIRECTION_MIXED,
    RHO_DECIMALS,
    VARIANCE_REGIME_UNDEFINED,
)

_TOL = 10.0 ** (-(RHO_DECIMALS - 1))


def validate_unique_keys(
    variance_budget: pd.DataFrame,
    residue_variance: pd.DataFrame,
    significance_screen: pd.DataFrame,
    domain_effect_summary: pd.DataFrame,
    report: S4Report,
) -> None:
    """Every output table's declared key is unique."""
    _assert_unique(
        variance_budget, ["serotype", "chain", "domain"], "variance_budget"
    )
    _assert_unique(
        residue_variance, ["serotype", "canon_label"], "residue_variance"
    )
    _assert_unique(
        significance_screen,
        ["serotype", "canon_label"],
        "significance_screen",
    )
    _assert_unique(
        domain_effect_summary,
        ["serotype", "chain", "domain"],
        "domain_effect_summary",
    )
    report.add(
        "every output table key is unique", "global", True, "4 tables checked"
    )


def validate_variance_fractions(
    variance_budget: pd.DataFrame,
    residue_variance: pd.DataFrame,
    report: S4Report,
) -> None:
    """τ²/σ̄² fractions sum to 1, lie in [0, 1], and match the regime label."""
    for name, table in (
        ("variance_budget", variance_budget),
        ("residue_variance", residue_variance),
    ):
        if table.empty:
            continue
        for row in table.itertuples(index=False):
            frac_tau2 = float(row.frac_tau2)
            frac_sigma2 = float(row.frac_sigma2)
            key = _row_key(name, row)
            if math.isnan(frac_tau2) or math.isnan(frac_sigma2):
                # only allowed when the regime is undefined (zero total variance)
                if row.variance_regime != VARIANCE_REGIME_UNDEFINED:
                    raise ConsistencyError(
                        f"{name} {key}: NaN variance fraction with regime "
                        f"{row.variance_regime!r} (expected "
                        f"{VARIANCE_REGIME_UNDEFINED!r})"
                    )
                continue
            for label, val in (("frac_tau2", frac_tau2), ("frac_sigma2", frac_sigma2)):
                if not (-_TOL <= val <= 1.0 + _TOL):
                    raise ConsistencyError(
                        f"{name} {key}: {label}={val} outside [0, 1]"
                    )
            if abs((frac_tau2 + frac_sigma2) - 1.0) > _TOL:
                raise ConsistencyError(
                    f"{name} {key}: frac_tau2 + frac_sigma2 = "
                    f"{frac_tau2 + frac_sigma2} != 1"
                )
            expected = classify_variance_regime(frac_tau2)
            if row.variance_regime != expected:
                raise ConsistencyError(
                    f"{name} {key}: regime {row.variance_regime!r} disagrees "
                    f"with frac_tau2={frac_tau2} (expected {expected!r})"
                )
    report.add(
        "variance fractions sum to 1, lie in [0,1], and match the regime label",
        "global",
        True,
        f"variance_budget={len(variance_budget)}, "
        f"residue_variance={len(residue_variance)}",
    )


def validate_significance_screen(
    significance_screen: pd.DataFrame, report: S4Report
) -> None:
    """Signed/CI/p-value fields are mutually consistent; BH is monotone."""
    if significance_screen.empty:
        report.add("significance screen self-consistent", "global", True, "empty")
        return
    for row in significance_screen.itertuples(index=False):
        key = f"({row.serotype}, {row.canon_label})"
        signed = bool(row.is_signed)
        if signed != (row.direction != DIRECTION_MIXED):
            raise ConsistencyError(
                f"significance_screen {key}: is_signed={signed} disagrees with "
                f"direction={row.direction!r}"
            )
        if not signed:
            # mixed mechanisms carry no signed CI / p-value / FDR decision
            if bool(row.ci_excludes_zero):
                raise ConsistencyError(
                    f"significance_screen {key}: unsigned row has "
                    f"ci_excludes_zero=True"
                )
            if bool(row.significant_fdr):
                raise ConsistencyError(
                    f"significance_screen {key}: unsigned row has "
                    f"significant_fdr=True"
                )
            if not _isnan(row.p_value_bh):
                raise ConsistencyError(
                    f"significance_screen {key}: unsigned row has a BH p-value"
                )
        else:
            if _isnan(row.p_value_bh):
                raise ConsistencyError(
                    f"significance_screen {key}: signed row is missing its BH "
                    f"p-value"
                )
        if bool(row.significant_raw) != bool(row.ci_excludes_zero):
            raise ConsistencyError(
                f"significance_screen {key}: significant_raw != ci_excludes_zero"
            )

    # BH monotonicity within each serotype: sorted by raw p, adjusted p is
    # non-decreasing.
    signed_rows = significance_screen[
        significance_screen["is_signed"].astype(bool)
    ]
    for serotype, grp in signed_rows.groupby("serotype", sort=False):
        ordered = grp.sort_values("p_value")
        prev = -math.inf
        for adj in ordered["p_value_bh"]:
            a = float(adj)
            if a + _TOL < prev:
                raise ConsistencyError(
                    f"significance_screen (serotype {serotype}): BH-adjusted "
                    f"p-values are not monotone in raw-p order"
                )
            prev = max(prev, a)
    report.add(
        "significance screen self-consistent (signed/CI/p-value/BH)",
        "global",
        True,
        f"{len(significance_screen)} mechanisms, {int(signed_rows.shape[0])} signed",
    )


def validate_domain_effect_totals(
    domain_effect_summary: pd.DataFrame, report: S4Report
) -> None:
    """Domain counts partition mechanisms; CI/FDR counts bounded by n_signed."""
    if domain_effect_summary.empty:
        report.add("domain effect totals partition", "global", True, "empty")
        return
    for row in domain_effect_summary.itertuples(index=False):
        key = f"({row.serotype}, {row.chain}, {row.domain})"
        if int(row.n_signed) + int(row.n_mixed) != int(row.n_mechanisms):
            raise ConsistencyError(
                f"domain_effect_summary {key}: n_signed + n_mixed != "
                f"n_mechanisms ({row.n_signed} + {row.n_mixed} != "
                f"{row.n_mechanisms})"
            )
        if int(row.n_ci_excludes_zero) > int(row.n_signed):
            raise ConsistencyError(
                f"domain_effect_summary {key}: n_ci_excludes_zero "
                f"{row.n_ci_excludes_zero} > n_signed {row.n_signed}"
            )
        if int(row.n_significant_fdr) > int(row.n_signed):
            raise ConsistencyError(
                f"domain_effect_summary {key}: n_significant_fdr "
                f"{row.n_significant_fdr} > n_signed {row.n_signed}"
            )
        frac = float(row.frac_ci_excludes_zero)
        if not math.isnan(frac) and not (-_TOL <= frac <= 1.0 + _TOL):
            raise ConsistencyError(
                f"domain_effect_summary {key}: frac_ci_excludes_zero {frac} "
                f"outside [0, 1]"
            )
    report.add(
        "domain effect totals partition mechanisms (counts bounded by n_signed)",
        "global",
        True,
        f"{len(domain_effect_summary)} (serotype, chain, domain) rows checked",
    )


def _row_key(name: str, row: object) -> str:
    if name == "variance_budget":
        return f"({row.serotype}, {row.chain}, {row.domain})"  # type: ignore[attr-defined]
    return f"({row.serotype}, {row.canon_label})"  # type: ignore[attr-defined]


def _isnan(value: object) -> bool:
    try:
        return math.isnan(float(value))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return False


def _assert_unique(df: pd.DataFrame, key: list[str], what: str) -> None:
    if df.empty:
        return
    dup = df.duplicated(key, keep=False)
    if dup.any():
        examples = df.loc[dup, key].drop_duplicates().head(3).to_dict("records")
        raise ConsistencyError(
            f"{what} key {tuple(key)} not unique; examples: {examples}"
        )
