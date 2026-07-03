"""S2 build subpackage: the five reusable per-serotype reduction builders.

Each builder is pure (no file IO, no mutation of inputs) and independently
importable by later stages. The shared deterministic rules (re-gating, tier
classification, the signed/significant screen) live in :mod:`._screens`.
"""
from __future__ import annotations

from .domain_reproducibility import build_domain_reproducibility
from .residue_landscape import build_residue_landscape
from .resolution_census import build_resolution_census
from .serotype_summary import build_serotype_summary
from .signed_screen import build_signed_screen

__all__ = [
    "build_resolution_census",
    "build_residue_landscape",
    "build_domain_reproducibility",
    "build_signed_screen",
    "build_serotype_summary",
]
