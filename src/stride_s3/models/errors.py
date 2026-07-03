"""Exception hierarchy for Stage S3.

Every failure subclasses :class:`S3Error`. Messages name the serotype/locus and
the offending value(s) so failures are actionable.
"""
from __future__ import annotations


class S3Error(Exception):
    """Base class for every Stage-S3 failure."""


class InputError(S3Error):
    """The S0 STRIDE table is missing, unreadable, or missing required columns.

    S3 consumes the S0 STRIDE table (the profile); this is raised when that input
    is absent or does not carry the columns S3 depends on.
    """


class ConsistencyError(S3Error):
    """A hierarchy-reduction invariant fails.

    Examples: a non-unique output key, a locus that does not carry exactly one
    profile row per scale, a domain-vs-residue gap that does not equal
    ``rho_domain - rho_residue``, a monotonicity flag that disagrees with the
    recomputed violation count, or a chain contrast whose direction counts do not
    partition the chain's mechanisms.
    """
