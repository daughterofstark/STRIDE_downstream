"""S2 IO subpackage: loaders for the S0/S1B inputs + writers for S2 artifacts."""
from __future__ import annotations

from .loaders import (
    file_digest,
    load_domain_annotation,
    load_residue_annotation,
    load_stride_table,
)
from .writers import write_reduction_summary, write_tables

__all__ = [
    "load_stride_table",
    "load_residue_annotation",
    "load_domain_annotation",
    "file_digest",
    "write_tables",
    "write_reduction_summary",
]
