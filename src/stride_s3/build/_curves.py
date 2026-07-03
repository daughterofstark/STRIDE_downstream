"""Deterministic hierarchy-reduction helpers shared by the S3 builders.

Every function here is a pure rule over a locus's ρ-vs-scale curve — computing
per-step gains, auditing monotonicity, and classifying a scale into its tier. No
file IO, no mutation, no randomness. Keeping them in one place makes the
reduction vocabulary auditable and guarantees the builders agree.

The load-bearing rules are taken verbatim from the design document:

- **Upward closure / monotonicity (I2)** (§3.4): ρ should be non-decreasing as
  the scale coarsens (from residue at index 0 to complex at index 6). A locus
  whose ρ ever drops when coarsening violates this and is flagged.
- **ρ is read, never recomputed** from its variance components (§1.2): these
  helpers only ever compare and difference the ρ values already in the profile.
"""
from __future__ import annotations

import math
from collections.abc import Sequence

from ..models.schema import (
    DOMAIN_SCALE_LEVEL,
    SCALE_LEVEL_TO_INDEX,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)


def scale_tier(scale_level: str) -> str:
    """Classify a scale into the licensed / exploratory tier.

    Domain scale and coarser is the *licensed* claim level at K=3; anything finer
    than domain (residue, secondary_structure, motif) is *exploratory*
    (§5.3–5.4).
    """
    domain_index = SCALE_LEVEL_TO_INDEX[DOMAIN_SCALE_LEVEL]
    idx = SCALE_LEVEL_TO_INDEX.get(scale_level)
    if idx is None:
        return TIER_EXPLORATORY
    return TIER_LICENSED if idx >= domain_index else TIER_EXPLORATORY


def step_gains(rho_by_index: Sequence[float]) -> list[float]:
    """Per-step ρ gain when coarsening: ``rho[i] - rho[i-1]``.

    Parameters
    ----------
    rho_by_index
        ρ ordered by scale index, finest (index 0) → coarsest. Assumed dense
        (one value per scale).

    Returns
    -------
    list[float]
        Same length as the input; the first element is ``nan`` (no finer scale
        to compare against), each subsequent element is the gain over the
        next-finer scale.
    """
    gains: list[float] = []
    for i, rho in enumerate(rho_by_index):
        if i == 0:
            gains.append(float("nan"))
        else:
            gains.append(float(rho) - float(rho_by_index[i - 1]))
    return gains


def audit_monotonicity(
    rho_by_index: Sequence[float],
) -> tuple[bool, int, float, int]:
    """Audit whether ρ is non-decreasing as the scale coarsens.

    Parameters
    ----------
    rho_by_index
        ρ ordered by scale index, finest → coarsest.

    Returns
    -------
    (is_monotone, n_violations, max_decrease, first_violation_index)
        ``is_monotone`` is True when ρ never decreases as the scale coarsens.
        ``n_violations`` counts adjacent steps where ρ decreases. ``max_decrease``
        is the largest single-step drop as a positive magnitude (0.0 if
        monotone). ``first_violation_index`` is the *finer* scale index of the
        first decreasing step, or ``-1`` if monotone.

    Notes
    -----
    ``NaN`` ρ values are skipped when forming adjacent comparisons (a gap in the
    curve neither counts as a violation nor masks a later one).
    """
    n_violations = 0
    max_decrease = 0.0
    first_violation_index = -1
    prev_idx: int | None = None
    prev_val: float | None = None
    for idx, rho in enumerate(rho_by_index):
        val = float(rho)
        if math.isnan(val):
            continue
        if prev_val is not None:
            delta = val - prev_val
            if delta < 0.0:
                n_violations += 1
                drop = -delta
                if drop > max_decrease:
                    max_decrease = drop
                if first_violation_index == -1 and prev_idx is not None:
                    first_violation_index = prev_idx
        prev_idx = idx
        prev_val = val
    is_monotone = n_violations == 0
    return is_monotone, n_violations, max_decrease, first_violation_index


def is_distributed_effect(
    rho_residue: float, rho_domain: float, rho_star: float
) -> bool:
    """True if a locus is reproducible only once aggregated to domain scale.

    A *distributed* effect is one that does not clear the (provisional) gate at
    residue scale but does at domain scale: ``rho_residue < rho_star ≤
    rho_domain`` (§3.4, "reproducible only when aggregated"). ``NaN`` ρ values
    never satisfy the comparison.
    """
    if _is_nan(rho_residue) or _is_nan(rho_domain):
        return False
    return rho_residue < rho_star <= rho_domain


def _is_nan(value: object) -> bool:
    return isinstance(value, float) and math.isnan(value)
