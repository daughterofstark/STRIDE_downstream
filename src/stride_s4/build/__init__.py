"""S4 build subpackage: the four reusable uncertainty-layer builders.

Each builder is pure (no file IO, no mutation of inputs) and independently
importable by later stages. The shared deterministic rules (the Wald p-value, the
Benjamini–Hochberg adjustment, the τ²/σ̄² variance split, the 1/β_se²-weighted
mean) live in :mod:`._stats`.
"""
from __future__ import annotations

from .domain_effect_summary import build_domain_effect_summary
from .residue_variance import build_residue_variance
from .significance_screen import build_significance_screen
from .variance_budget import build_variance_budget

__all__ = [
    "build_variance_budget",
    "build_residue_variance",
    "build_significance_screen",
    "build_domain_effect_summary",
]
