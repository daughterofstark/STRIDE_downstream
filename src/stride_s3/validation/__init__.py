"""S3 validation subpackage: structural hierarchy-reduction checks."""
from __future__ import annotations

from .checks import (
    validate_chain_contrast_totals,
    validate_gap_consistency,
    validate_monotonicity_audit_consistency,
    validate_scale_curve_completeness,
    validate_unique_keys,
)

__all__ = [
    "validate_unique_keys",
    "validate_scale_curve_completeness",
    "validate_gap_consistency",
    "validate_monotonicity_audit_consistency",
    "validate_chain_contrast_totals",
]
