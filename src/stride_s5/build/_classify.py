"""Deterministic cross-serotype label helpers shared by the S5 builders.

Every function here is a pure classification rule over already-aggregated
serotype-level counts — the conservation-of-reproducibility label, the
direction-concordance label, and the catalytic-machinery membership tests. No
file IO, no mutation, no randomness. Keeping them in one place makes the
cross-serotype vocabulary auditable and guarantees the builders and the
validators agree.

The load-bearing rules follow the design document:

- **conservation of reproducibility** (§3.3): a shared position is reproducible in
  all / a majority / some / none of the serotypes it is present in;
- **direction concordance** (§3.3): for positions signed in ≥ 2 serotypes, the
  serotypes agree, split by a majority, or conflict on the sign;
- **conserved catalytic machinery** (§3.3): the Catalytic Triad residues and the
  catalytic-machinery domains are flagged so their conservation can be inspected
  directly.
"""
from __future__ import annotations

from ..models.schema import (
    CATALYTIC_DOMAINS,
    CATALYTIC_TRIAD_LABELS,
    CONCORDANCE_AGREE,
    CONCORDANCE_CONFLICT,
    CONCORDANCE_MAJORITY,
    CONSERVATION_ALL,
    CONSERVATION_MAJORITY,
    CONSERVATION_NONE,
    CONSERVATION_SOME,
    DIRECTION_DECREASE,
    DIRECTION_INCREASE,
    DIRECTION_NONE,
)


def is_catalytic_triad(canon_label: str) -> bool:
    """True if ``canon_label`` is one of the Catalytic Triad residues."""
    return canon_label in CATALYTIC_TRIAD_LABELS


def is_catalytic_domain(domain: str) -> bool:
    """True if ``domain`` is a catalytic-machinery domain (Triad / Oxyanion Loop)."""
    return domain in CATALYTIC_DOMAINS


def conservation_class(n_reproducible: int, n_present: int) -> str:
    """Label how conserved a position's reproducibility is across serotypes.

    Parameters
    ----------
    n_reproducible
        Serotypes in which the position is present *and* reproducible (ρ ≥ ρ\\*).
    n_present
        Serotypes in which the position is present (the denominator).

    Returns
    -------
    str
        ``reproducible_all`` when reproducible in every serotype it is present in;
        ``reproducible_none`` when reproducible in none; ``reproducible_majority``
        when in a strict majority; ``reproducible_some`` otherwise. When
        ``n_present == 0`` the position is treated as ``reproducible_none``.
    """
    if n_present <= 0 or n_reproducible <= 0:
        return CONSERVATION_NONE
    if n_reproducible >= n_present:
        return CONSERVATION_ALL
    if n_reproducible * 2 > n_present:
        return CONSERVATION_MAJORITY
    return CONSERVATION_SOME


def concordance_class(n_increase: int, n_decrease: int) -> str:
    """Label the direction concordance of a position's signed serotypes.

    Parameters
    ----------
    n_increase
        Serotypes where the position is signed with an ``increase``.
    n_decrease
        Serotypes where the position is signed with a ``decrease``.

    Returns
    -------
    str
        ``agree`` when only one direction occurs; ``conflict`` when the two
        directions tie; ``majority`` when both occur but one is a strict majority.
    """
    if n_increase == 0 or n_decrease == 0:
        return CONCORDANCE_AGREE
    if n_increase == n_decrease:
        return CONCORDANCE_CONFLICT
    return CONCORDANCE_MAJORITY


def majority_direction(n_increase: int, n_decrease: int) -> str:
    """The strict-majority direction, or ``none`` on a tie or when unsigned."""
    if n_increase == 0 and n_decrease == 0:
        return DIRECTION_NONE
    if n_increase > n_decrease:
        return DIRECTION_INCREASE
    if n_decrease > n_increase:
        return DIRECTION_DECREASE
    return DIRECTION_NONE
