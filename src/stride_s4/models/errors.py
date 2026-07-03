"""Exception hierarchy for Stage S4.

Every failure subclasses :class:`S4Error`. Messages name the serotype/region and
the offending value(s) so failures are actionable.
"""
from __future__ import annotations


class S4Error(Exception):
    """Base class for every Stage-S4 failure."""


class InputError(S4Error):
    """The S0 STRIDE table is missing, unreadable, or missing required columns.

    S4 consumes the S0 STRIDE table (the tidy profile); this is raised when that
    input is absent or does not carry the columns S4 depends on.
    """


class ConsistencyError(S4Error):
    """An uncertainty-layer invariant fails.

    Examples: a non-unique output key, a variance fraction that does not lie in
    ``[0, 1]`` or whose τ²/σ̄² split does not sum to 1, a regime label that
    disagrees with the recomputed fraction, a significance row whose
    signed/CI/p-value fields are mutually inconsistent, or a domain effect
    summary whose direction counts do not partition its mechanisms.
    """
