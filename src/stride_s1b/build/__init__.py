"""S1B build subpackage: the four reusable annotation-table builders.

Each builder is pure (no file IO, no mutation of inputs) and independently
importable by later stages.
"""
from __future__ import annotations

from .domain_annotation import build_domain_annotation
from .hierarchy_annotation import build_hierarchy_annotation
from .residue_annotation import build_residue_annotation
from .serotype_annotation import build_serotype_annotation

__all__ = [
    "build_residue_annotation",
    "build_domain_annotation",
    "build_hierarchy_annotation",
    "build_serotype_annotation",
]
