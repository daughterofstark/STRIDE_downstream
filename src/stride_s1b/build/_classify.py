"""Deterministic classification helpers shared by the S1B builders.

Every function here maps a structural fact from S1A to a fixed category label.
These are pure lookups and boolean rules — no statistics, no scoring, no
ranking. Keeping them in one place makes the annotation vocabulary auditable and
guarantees the builders agree on the rules.
"""
from __future__ import annotations

from ..models.schema import (
    AVAIL_ALL,
    AVAIL_NONE,
    AVAIL_SOME,
    CONSERVATION_PAN,
    CONSERVATION_PARTIAL,
    CONSERVATION_UNIQUE,
    DOMAIN_STATUS_ASSIGNED,
    DOMAIN_STATUS_UNASSIGNED,
    SS_STATUS_RESOLVED,
    SS_STATUS_UNRESOLVED,
    UNRESOLVED_HIERARCHY_SENTINELS,
)


def is_resolved(level_value: object) -> bool:
    """Return True if a hierarchy level value is a real assignment.

    A level is *unresolved* if it is one of STRIDE's sentinels
    (``unassigned``/``none``/``unknown``/empty). The comparison is
    case-insensitive and whitespace-tolerant; the underlying value is never
    modified.
    """
    if level_value is None:
        return False
    return str(level_value).strip().lower() not in UNRESOLVED_HIERARCHY_SENTINELS


def domain_status(domain_value: str) -> str:
    """Classify a domain value as ``assigned`` or ``unassigned``."""
    return DOMAIN_STATUS_ASSIGNED if is_resolved(domain_value) else DOMAIN_STATUS_UNASSIGNED


def secondary_structure_status(ss_value: str) -> str:
    """Classify a secondary-structure value as ``resolved`` or ``unresolved``."""
    return SS_STATUS_RESOLVED if is_resolved(ss_value) else SS_STATUS_UNRESOLVED


def conservation_class(n_serotypes_present: int, n_serotypes_total: int) -> str:
    """Classify a residue's conservation from its serotype-presence count.

    - ``pan_serotype``   : present in every serotype (and there is >1 serotype,
      or the single-serotype dataset trivially has it everywhere);
    - ``serotype_unique``: present in exactly one serotype (of several);
    - ``partial``        : present in some-but-not-all serotypes.
    """
    if n_serotypes_total <= 0:
        return CONSERVATION_UNIQUE
    if n_serotypes_present >= n_serotypes_total:
        return CONSERVATION_PAN
    if n_serotypes_present <= 1:
        return CONSERVATION_UNIQUE
    return CONSERVATION_PARTIAL


def availability_class(
    n_replicates: int, available: bool, in_all_replicates: bool
) -> str:
    """Classify a residue's replicate availability.

    - ``all_replicates`` : observed in every replicate of its serotype;
    - ``some_replicates``: observed in at least one but not all replicates;
    - ``no_replicates``  : not observed in any replicate.
    """
    if not available or n_replicates <= 0:
        return AVAIL_NONE
    if in_all_replicates:
        return AVAIL_ALL
    return AVAIL_SOME
