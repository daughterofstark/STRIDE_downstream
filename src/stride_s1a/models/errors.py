"""Exception hierarchy for Stage S1A.

Every failure subclasses :class:`S1AError`. Messages name the serotype/residue
and the offending value(s) so failures are actionable.
"""
from __future__ import annotations


class S1AError(Exception):
    """Base class for every Stage-S1A failure."""


class InputError(S1AError):
    """An S0 canonical table is missing, unreadable, or missing required columns.

    S1A consumes only the canonical parquet tables; this is raised when those
    inputs are absent or do not carry the columns S1A depends on.
    """


class ConsistencyError(S1AError):
    """A derived-dataset invariant fails.

    Examples: a non-unique canonical residue key, a STRIDE locus that maps to no
    canonical residue, a residue with more than one hierarchy path, a replicate
    row that maps to no canonical residue, or an orphan residue.
    """
