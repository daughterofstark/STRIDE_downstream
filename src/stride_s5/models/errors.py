"""Exception hierarchy for Stage S5.

Every failure subclasses :class:`S5Error`. Messages name the serotype/position or
region and the offending value(s) so failures are actionable.
"""
from __future__ import annotations


class S5Error(Exception):
    """Base class for every Stage-S5 failure."""


class InputError(S5Error):
    """An S5 input is missing, unreadable, or missing required columns.

    S5 consumes the S0 STRIDE table (the tidy profile) and the S1A
    ``conservation_table`` (the shared-position index); this is raised when either
    input is absent or does not carry the columns S5 depends on.
    """


class ConsistencyError(S5Error):
    """A cross-serotype-layer invariant fails.

    Examples: a non-unique output key; a conservation row whose
    ``n_serotypes_reproducible`` exceeds ``n_serotypes_present`` or whose
    ``conservation_class`` disagrees with the recomputed reproducible fraction; a
    concordance row whose direction counts do not partition its signed serotypes
    or whose ``concordance_class`` disagrees with those counts; a domain × serotype
    row whose ρ is not region-constant across the domain's member loci; or a
    scorecard whose signed / mixed counts do not partition its mechanisms.
    """
