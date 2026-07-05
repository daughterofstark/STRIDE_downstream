"""Exception hierarchy for Stage S6.

Every failure subclasses :class:`S6Error`. Messages name the serotype / position
and the offending value(s) so failures are actionable.
"""
from __future__ import annotations


class S6Error(Exception):
    """Base class for every Stage-S6 failure."""


class InputError(S6Error):
    """An S6 input is missing, unreadable, or missing required columns.

    S6's required input is the S1A ``replicate_inventory`` (the replicate-structure
    index). The S0 ``replicate_table`` (the Level-1 per-run observations) is
    *optional* — its absence is not an error but the design's anticipated
    "per-run correlation CSVs unavailable" state, in which the per-run analyses are
    recorded as blocked. This is raised when the inventory is absent/unreadable, or
    when either input is present but does not carry the columns S6 depends on.
    """


class ConsistencyError(S6Error):
    """A replicate-layer invariant fails.

    Examples: a non-unique output key; a regime row whose
    ``per_replicate_effects_available`` disagrees with its
    ``n_replicates_with_effects`` or whose ``residue_claims_licensed`` disagrees
    with ``n_replicates``; a spread row whose ``theta_max < theta_min`` or whose
    range / absolute-mean is negative; a concordance row whose Kendall's *W* is
    outside ``[0, 1]`` or whose ``concordance_class`` disagrees with the recomputed
    coefficient; or a blocked-analysis row whose ``available`` flag disagrees with
    its ``status`` (or a leave-one-replicate-out row not marked blocked).
    """
