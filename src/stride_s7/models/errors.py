"""Exception hierarchy for Stage S7.

Every failure subclasses :class:`S7Error`. Messages name the offending
figure/table id and input so failures are actionable.
"""
from __future__ import annotations


class S7Error(Exception):
    """Base class for every Stage-S7 failure."""


class InputError(S7Error):
    """An S7 input is missing, unreadable, or missing required columns.

    S7 reads the S2–S6 reduction tables (never the raw STRIDE files). This is
    raised when a required prior-stage table is absent or unreadable, or is present
    but does not carry the columns S7 depends on. A required input's absence is a
    hard failure (fail loudly) — unlike an *empty but present* input, which is
    handled gracefully by emitting an empty-but-valid artifact.
    """


class ConsistencyError(S7Error):
    """A reporting-layer structural invariant fails.

    Examples: a prepared figure-data table missing its declared columns; a
    manuscript table missing a declared column; a non-deterministic or missing
    output filename; an artifact recorded in the manifest but not written to disk;
    or an incomplete provenance header (missing input digests, ρ\\*, or the
    calibrated flag).
    """
