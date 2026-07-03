"""Exception hierarchy for Stage S2.

Every failure subclasses :class:`S2Error`. Messages name the serotype/region and
the offending value(s) so failures are actionable.
"""
from __future__ import annotations


class S2Error(Exception):
    """Base class for every Stage-S2 failure."""


class InputError(S2Error):
    """An input table is missing, unreadable, or missing required columns.

    S2 consumes the S0 STRIDE table and the S1B annotation tables; this is
    raised when those inputs are absent or do not carry the columns S2 depends
    on.
    """


class ConfigError(S2Error):
    """The requested reduction configuration is invalid.

    Examples: an empty ρ* band, or a ρ* value outside the ``[0, 1]`` range that
    the reproducibility coefficient is defined on.
    """


class ConsistencyError(S2Error):
    """A per-serotype reduction invariant fails.

    Examples: a non-unique output key, a census whose per-scale counts do not
    sum to the serotype's locus count, a domain-scale ρ that is not a
    region-level constant across the domain's member loci, or a serotype
    scorecard whose direction counts do not partition the mechanisms.
    """
