"""Deterministic uncertainty helpers shared by the S4 builders.

Every function here is a pure numeric rule — the two-sided Wald p-value, the
Benjamini–Hochberg step-up adjustment, the τ²/σ̄² variance split and its regime
label, and the 1/β_se²-weighted mean. No file IO, no mutation, no randomness, and
**no SciPy dependency** (the normal CDF is computed from :func:`math.erf`).
Keeping these in one place makes the uncertainty vocabulary auditable and
guarantees the builders and the validators agree.

The load-bearing rules follow the design document:

- **CI-based significance** (§3.5): a signed effect is significant when its
  confidence interval excludes 0. For a normal (Wald) interval this is equivalent
  to ``|β/β_se| ≥ z_{1-α/2}``; the same z-score yields the two-sided p-value.
- **FDR control** (§5.2): when screening the ~200 positions of a serotype, control
  the false-discovery rate with Benjamini–Hochberg rather than reporting raw
  counts.
- **τ²/σ̄² diagnostic** (§3.1): a region with a high τ² fraction is
  replicate-dominated (candidate metastability); a low fraction is
  sampling-dominated.
- **β_se-weighted summaries** (§3.5): average effects within a domain weighting by
  1/β_se².
"""
from __future__ import annotations

import math

from ..models.schema import (
    REPLICATE_DOMINATED_FRACTION,
    SAMPLING_DOMINATED_FRACTION,
    VARIANCE_REGIME_BALANCED,
    VARIANCE_REGIME_REPLICATE,
    VARIANCE_REGIME_SAMPLING,
    VARIANCE_REGIME_UNDEFINED,
)


def normal_cdf(z: float) -> float:
    """Standard-normal cumulative distribution function Φ(z).

    Computed from :func:`math.erf` so S4 needs no SciPy:
    ``Φ(z) = ½(1 + erf(z / √2))``.

    Parameters
    ----------
    z
        The quantile.

    Returns
    -------
    float
        ``Φ(z)`` in ``[0, 1]``.
    """
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def two_sided_wald_p(beta: float, beta_se: float) -> float:
    """Two-sided Wald p-value for ``H0: effect = 0``.

    Parameters
    ----------
    beta
        The signed effect estimate.
    beta_se
        Its standard error; must be strictly positive for a defined p-value.

    Returns
    -------
    float
        ``2·(1 − Φ(|β / β_se|))``, or ``nan`` if either input is ``nan`` or
        ``beta_se`` is not strictly positive.
    """
    if _is_nan(beta) or _is_nan(beta_se) or beta_se <= 0.0:
        return float("nan")
    z = abs(beta / beta_se)
    return 2.0 * (1.0 - normal_cdf(z))


def wald_z(beta: float, beta_se: float) -> float:
    """Signed Wald z-score ``β / β_se`` (``nan`` if undefined)."""
    if _is_nan(beta) or _is_nan(beta_se) or beta_se <= 0.0:
        return float("nan")
    return beta / beta_se


def ci_excludes_zero(ci_lower: float, ci_upper: float) -> bool:
    """True if the closed interval ``[ci_lower, ci_upper]`` excludes 0.

    A ``nan`` bound (e.g. a mixed mechanism with no signed CI) yields ``False``.
    """
    if _is_nan(ci_lower) or _is_nan(ci_upper):
        return False
    return ci_lower > 0.0 or ci_upper < 0.0


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Benjamini–Hochberg step-up adjusted p-values.

    Parameters
    ----------
    p_values
        Raw p-values. ``nan`` entries are carried through as ``nan`` and excluded
        from the multiple-testing family (they neither receive nor affect an
        adjusted value).

    Returns
    -------
    list[float]
        Adjusted p-values aligned to the input order, each clamped to ``[0, 1]``
        and enforced monotone non-decreasing in the rank order of the raw
        p-values (the standard BH cumulative-minimum from the largest rank down).

    Notes
    -----
    With ``m`` the number of non-``nan`` p-values, the value ranked ``i`` (1-based,
    ascending) is adjusted to ``p·m/i`` and then made monotone so a smaller raw
    p-value never yields a larger adjusted one.
    """
    indexed = [(i, p) for i, p in enumerate(p_values) if not _is_nan(p)]
    out: list[float] = [float("nan")] * len(p_values)
    m = len(indexed)
    if m == 0:
        return out
    # sort by raw p ascending
    indexed.sort(key=lambda t: t[1])
    # step-up: walk from the largest rank down, taking the cumulative minimum
    running_min = float("inf")
    adjusted_by_pos: dict[int, float] = {}
    for rank in range(m, 0, -1):
        orig_index, p = indexed[rank - 1]
        val = p * m / rank
        running_min = min(running_min, val)
        adjusted_by_pos[orig_index] = min(1.0, max(0.0, running_min))
    for idx, val in adjusted_by_pos.items():
        out[idx] = val
    return out


def variance_fractions(
    tau2: float, sigma2_bar: float
) -> tuple[float, float, float, float, str]:
    """Split the unreproduced variance into τ² and σ̄² fractions with a regime.

    Parameters
    ----------
    tau2
        Between-replicate variance component τ² (replicate disagreement).
    sigma2_bar
        Autocorrelation-corrected within-replicate variance σ̄² (sampling noise).

    Returns
    -------
    (total, frac_tau2, frac_sigma2, ratio, regime)
        ``total`` is ``τ² + σ̄²``; ``frac_tau2`` / ``frac_sigma2`` are the shares
        of that total (``nan`` when the total is 0); ``ratio`` is ``τ²/σ̄²``
        (``nan`` when σ̄² is 0); ``regime`` is one of the design's labels —
        replicate-dominated, sampling-dominated, balanced, or undefined (when the
        total is 0).
    """
    if _is_nan(tau2) or _is_nan(sigma2_bar):
        return (float("nan"), float("nan"), float("nan"), float("nan"),
                VARIANCE_REGIME_UNDEFINED)
    total = tau2 + sigma2_bar
    if total <= 0.0:
        return (total, float("nan"), float("nan"), float("nan"),
                VARIANCE_REGIME_UNDEFINED)
    frac_tau2 = tau2 / total
    frac_sigma2 = sigma2_bar / total
    ratio = tau2 / sigma2_bar if sigma2_bar > 0.0 else float("nan")
    regime = classify_variance_regime(frac_tau2)
    return total, frac_tau2, frac_sigma2, ratio, regime


def classify_variance_regime(frac_tau2: float) -> str:
    """Label a τ² fraction as replicate/sampling-dominated or balanced.

    ``frac_tau2 ≥ REPLICATE_DOMINATED_FRACTION`` → replicate-dominated;
    ``frac_tau2 ≤ SAMPLING_DOMINATED_FRACTION`` → sampling-dominated;
    strictly between → balanced; ``nan`` → undefined.
    """
    if _is_nan(frac_tau2):
        return VARIANCE_REGIME_UNDEFINED
    if frac_tau2 >= REPLICATE_DOMINATED_FRACTION:
        return VARIANCE_REGIME_REPLICATE
    if frac_tau2 <= SAMPLING_DOMINATED_FRACTION:
        return VARIANCE_REGIME_SAMPLING
    return VARIANCE_REGIME_BALANCED


def weighted_mean_se(
    values: list[float], ses: list[float]
) -> tuple[float, float, float]:
    """Inverse-variance (1/β_se²) weighted mean of signed effects.

    Parameters
    ----------
    values
        Signed effect estimates.
    ses
        Their standard errors; only strictly-positive, non-``nan`` pairs
        contribute.

    Returns
    -------
    (weighted_mean, weighted_se, unweighted_mean)
        The 1/se²-weighted mean, its standard error ``sqrt(1/Σ w)`` with
        ``w = 1/se²``, and the plain arithmetic mean over the contributing
        values. All three are ``nan`` when no pair contributes.
    """
    num = 0.0
    den = 0.0
    plain: list[float] = []
    for v, se in zip(values, ses, strict=True):
        if _is_nan(v) or _is_nan(se) or se <= 0.0:
            continue
        w = 1.0 / (se * se)
        num += w * v
        den += w
        plain.append(v)
    if den <= 0.0 or not plain:
        return float("nan"), float("nan"), float("nan")
    weighted_mean = num / den
    weighted_se = math.sqrt(1.0 / den)
    unweighted_mean = sum(plain) / len(plain)
    return weighted_mean, weighted_se, unweighted_mean


def _is_nan(value: object) -> bool:
    return isinstance(value, float) and math.isnan(value)
