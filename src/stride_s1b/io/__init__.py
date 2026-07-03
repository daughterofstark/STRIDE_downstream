"""S1B IO subpackage: loaders for S1A tables + writers for S1B artifacts."""
from __future__ import annotations

from .loaders import (
    load_canonical_residues,
    load_conservation_table,
    load_domain_table,
    load_replicate_inventory,
)
from .writers import write_annotation_summary, write_tables

__all__ = [
    "load_canonical_residues",
    "load_domain_table",
    "load_replicate_inventory",
    "load_conservation_table",
    "write_tables",
    "write_annotation_summary",
]
