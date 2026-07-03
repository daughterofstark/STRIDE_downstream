"""S1B validation subpackage: structural annotation checks."""
from __future__ import annotations

from .checks import (
    validate_domain_membership,
    validate_one_annotation_per_residue,
    validate_referential_integrity,
    validate_serotype_references,
    validate_unique_hierarchy_paths,
)

__all__ = [
    "validate_one_annotation_per_residue",
    "validate_unique_hierarchy_paths",
    "validate_domain_membership",
    "validate_serotype_references",
    "validate_referential_integrity",
]
