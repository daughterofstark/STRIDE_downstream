"""Exception hierarchy for the analysis framework.

Every failure subclasses :class:`StrideAnalysisError`, so callers can catch all
framework errors with one ``except``. Messages name the serotype/replicate/file
and the offending value(s) so failures are actionable.
"""
from __future__ import annotations


class StrideAnalysisError(Exception):
    """Base class for every framework failure."""


class DiscoveryError(StrideAnalysisError):
    """A dataset directory cannot be resolved into the expected structure.

    Examples: missing profile/mechanism, missing replicate directories,
    inconsistent replicate counts across serotypes, duplicate serotypes,
    malformed run-directory layout, or an empty dataset.
    """


class SchemaError(StrideAnalysisError):
    """A file violates the frozen schema for its data level.

    Covers: missing/extra columns or JSON keys, wrong dtypes, out-of-range
    values, missing required identifiers, illegal enum values, and the
    null-iff-mixed rule for mechanism beta fields.
    """


class HierarchyError(StrideAnalysisError):
    """A ``region_id`` cannot be parsed against the frozen hierarchy grammar."""


class ConsistencyError(StrideAnalysisError):
    """A cross-file / cross-row invariant fails.

    Examples: mechanism ``rho`` disagreeing with its gated profile row, more
    than one gated row per locus, orphan mechanisms or orphan gated loci, or a
    canonical-table key that is not unique.
    """
