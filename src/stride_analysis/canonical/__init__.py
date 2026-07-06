"""Canonical subpackage: builds the two separate canonical tables.

The replicate table (Level 1) and STRIDE table (Level 2) are deliberately kept
separate — later stages import whichever they need.
"""
from __future__ import annotations

from .replicate_table import (
    assemble_replicate_table,
    build_replicate_rows,
    replicate_canon_label,
)
from .stride_table import assemble_stride_table, build_stride_rows

__all__ = [
    "build_replicate_rows",
    "replicate_canon_label",
    "assemble_replicate_table",
    "build_stride_rows",
    "assemble_stride_table",
]
