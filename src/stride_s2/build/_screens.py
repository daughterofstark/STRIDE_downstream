"""Deterministic reduction helpers shared by the S2 builders.

Every function here is a pure rule over profile/mechanism facts — re-gating a
locus at a given ρ*, deciding a tier, and evaluating the signed/significant
screen. No file IO, no mutation, no randomness. Keeping them in one place makes
the reduction vocabulary auditable and guarantees the builders agree.

Two rules are load-bearing and taken verbatim from the design document:

- **Re-gating** (§1.5): a locus gates at the *finest* scale (smallest scale
  index) at which ρ ≥ ρ*. Because the profile already carries every scale,
  re-gating at any ρ* is a filter, not a recomputation (§5.3). ρ itself is never
  recomputed from its variance components (§1.2).
- **Signed/significant screen** (§2.2, §3.5): a mechanism yields a signed,
  significant claim when its direction is not ``mixed`` **and** its β confidence
  interval excludes 0 **and** its ρ meets the (provisional) threshold.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Mapping

from ..models.schema import (
    DIRECTION_MIXED,
    DOMAIN_SCALE_LEVEL,
    RESIDUE_SCALE_LEVEL,
    SCALE_LEVEL_TO_INDEX,
    SCALE_UNRESOLVED,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)


def regate_scale(
    rho_by_scale: Mapping[str, float], rho_star: float
) -> tuple[str, int]:
    """Return the finest scale at which ρ ≥ ρ* for one locus.

    Parameters
    ----------
    rho_by_scale
        Mapping ``scale_level -> ρ`` for a single locus (one entry per scale).
    rho_star
        The threshold to gate at.

    Returns
    -------
    (gated_scale_level, gated_scale_index)
        The finest (smallest-index) scale whose ρ ≥ ρ*. If no scale reaches ρ*,
        returns the :data:`SCALE_UNRESOLVED` sentinel with index ``-1``.

    Notes
    -----
    ``NaN`` ρ values never satisfy ``ρ ≥ ρ*`` and are therefore skipped.
    """
    best_index: int | None = None
    best_level = SCALE_UNRESOLVED
    for level, rho in rho_by_scale.items():
        if rho is None or (isinstance(rho, float) and math.isnan(rho)):
            continue
        if rho >= rho_star:
            idx = SCALE_LEVEL_TO_INDEX[level]
            if best_index is None or idx < best_index:
                best_index = idx
                best_level = level
    if best_index is None:
        return SCALE_UNRESOLVED, -1
    return best_level, best_index


def scale_tier(scale_level: str) -> str:
    """Classify a gated scale into the licensed / exploratory tier.

    Domain scale and coarser is the *licensed* claim level at K=3; anything
    finer than domain (residue, secondary_structure, motif) is *exploratory*
    (§5.3–5.4). The unresolved sentinel is exploratory (no licensed claim).
    """
    if scale_level == SCALE_UNRESOLVED:
        return TIER_EXPLORATORY
    domain_index = SCALE_LEVEL_TO_INDEX[DOMAIN_SCALE_LEVEL]
    idx = SCALE_LEVEL_TO_INDEX.get(scale_level)
    if idx is None:
        return TIER_EXPLORATORY
    return TIER_LICENSED if idx >= domain_index else TIER_EXPLORATORY


def is_residue_scale(scale_level: str) -> bool:
    """True if the gated scale is the finest (residue) scale."""
    return scale_level == RESIDUE_SCALE_LEVEL


def is_signed(direction: object) -> bool:
    """True if a mechanism carries a coherent signed direction (not ``mixed``)."""
    if direction is None:
        return False
    return str(direction) != DIRECTION_MIXED


def ci_excludes_zero(
    beta_ci_lower: float | None, beta_ci_upper: float | None
) -> bool:
    """True if a signed β confidence interval strictly excludes 0.

    A ``None`` bound (the mechanism is mixed, so no signed CI) never excludes 0.
    An interval touching 0 at either bound does not exclude it.
    """
    if beta_ci_lower is None or beta_ci_upper is None:
        return False
    if _is_nan(beta_ci_lower) or _is_nan(beta_ci_upper):
        return False
    lo, hi = sorted((float(beta_ci_lower), float(beta_ci_upper)))
    return lo > 0.0 or hi < 0.0


def passes_signed_screen(
    direction: object,
    beta_ci_lower: float | None,
    beta_ci_upper: float | None,
    rho: float,
    rho_star: float,
) -> bool:
    """The full signed/significant screen: signed AND CI≠0 AND ρ ≥ ρ*."""
    if not is_signed(direction):
        return False
    if not ci_excludes_zero(beta_ci_lower, beta_ci_upper):
        return False
    if rho is None or _is_nan(rho):
        return False
    return rho >= rho_star


def normalise_band(band: Iterable[float], decimals: int) -> tuple[float, ...]:
    """Round, de-duplicate, and sort a ρ* band ascending.

    Deterministic: identical input multisets always yield the same tuple.
    """
    seen: set[float] = set()
    out: list[float] = []
    for value in band:
        rounded = round(float(value), decimals)
        if rounded not in seen:
            seen.add(rounded)
            out.append(rounded)
    return tuple(sorted(out))


def _is_nan(value: object) -> bool:
    return isinstance(value, float) and math.isnan(value)
