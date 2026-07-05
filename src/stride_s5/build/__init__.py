"""S5 build subpackage: the four reusable cross-serotype builders.

Each builder is pure (no file IO, no mutation of inputs) and independently
importable by later stages. The shared deterministic rules live in the private
helpers: :mod:`._frames` (per-serotype aggregation of the profile — the
anti-pseudoreplication "aggregate first" step) and :mod:`._classify` (the
conservation, concordance, and catalytic-machinery label vocabularies).
"""
from __future__ import annotations

from .cross_serotype_scorecard import build_cross_serotype_scorecard
from .direction_concordance import build_direction_concordance
from .domain_serotype_matrix import build_domain_serotype_matrix
from .position_conservation import build_position_conservation

__all__ = [
    "build_position_conservation",
    "build_direction_concordance",
    "build_domain_serotype_matrix",
    "build_cross_serotype_scorecard",
]
