"""Exception hierarchy for Stage S1B.

Every failure subclasses :class:`S1BError`. Messages name the serotype/residue
and the offending value(s) so failures are actionable.
"""
from __future__ import annotations


class S1BError(Exception):
    """Base class for every Stage-S1B failure."""


class InputError(S1BError):
    """An S1A table is missing, unreadable, or missing required columns.

    S1B consumes only the S1A parquet tables; this is raised when those inputs
    are absent or do not carry the columns S1B depends on.
    """


class ConsistencyError(S1BError):
    """A structural annotation invariant fails.

    Examples: a residue with more than one biological annotation, a non-unique
    hierarchy path, an internally inconsistent domain membership, an orphan
    annotation, a serotype referencing a non-existent canonical residue, or a
    referential-integrity break between two generated tables.
    """
