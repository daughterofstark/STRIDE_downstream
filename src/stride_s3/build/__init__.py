"""S3 build subpackage: the four reusable hierarchy-reduction builders.

Each builder is pure (no file IO, no mutation of inputs) and independently
importable by later stages. The shared deterministic rules (per-step gain, the
monotonicity audit, tier classification, the distributed-effect flag) live in
:mod:`._curves`.
"""
from __future__ import annotations

from .chain_contrast import build_chain_contrast
from .monotonicity_audit import build_monotonicity_audit
from .resolution_gap import build_resolution_gap
from .scale_curve import build_scale_curve

__all__ = [
    "build_scale_curve",
    "build_resolution_gap",
    "build_monotonicity_audit",
    "build_chain_contrast",
]
