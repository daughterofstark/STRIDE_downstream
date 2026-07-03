"""S1A build subpackage: the four reusable-table builders.

Each builder is pure (no file IO, no mutation of inputs) and independently
importable by later stages.
"""
from __future__ import annotations

from .canonical_residues import build_canonical_residues
from .conservation_table import build_conservation_table
from .domain_table import build_domain_table
from .replicate_inventory import build_replicate_inventory

__all__ = [
    "build_canonical_residues",
    "build_domain_table",
    "build_replicate_inventory",
    "build_conservation_table",
]
