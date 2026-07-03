"""stride_s3 — Stage S3: the hierarchy reduction layer.

S3 consumes the S0 STRIDE table (``stride_table.parquet``) — the profile, with
every locus evaluated at all seven hierarchy scales — and reduces it *along the
scale axis* into per-locus and per-chain hierarchy products:

- ``scale_curve.parquet``        — the ρ-vs-scale curve for every locus (F7)
- ``resolution_gap.parquet``     — the domain-vs-residue reproducibility gap Δρ
- ``monotonicity_audit.parquet`` — the upward-closure (I2) audit per locus
- ``chain_contrast.parquet``     — chain-level contrast (e.g. NS2B vs NS3)
- ``hierarchy_summary.json``     — facts, provenance header, validation outcomes

S3's products are **ρ*-independent** descriptions of the profile: it makes no
calibrated pass/fail resolution claim (the gate is uncalibrated, §0.1), never
recomputes ρ from its variance components, performs **no cross-serotype tests**
(that is S5), and produces no figures. It never modifies S0, S1A, S1B, or S2.

Two tiers are labelled explicitly per the design: **licensed** (domain-scale and
coarser, the claim level at K=3) and **exploratory** (residue-scale, outside the
operating range).

Public API
----------
- :func:`run_s3`   — load → build → validate → write artifacts.
- :func:`build_s3` — same, in memory (no writes).
- :class:`S3Tables`— the four returned tables.
- :class:`S3Error` — base of every S3 exception.
"""
from __future__ import annotations

from .build import (
    build_chain_contrast,
    build_monotonicity_audit,
    build_resolution_gap,
    build_scale_curve,
)
from .models import S3Report
from .models.errors import ConsistencyError, InputError, S3Error
from .s3 import S3Tables, build_s3, run_s3

__all__ = [
    # orchestration
    "run_s3",
    "build_s3",
    "S3Tables",
    # builders (reusable by later stages)
    "build_scale_curve",
    "build_resolution_gap",
    "build_monotonicity_audit",
    "build_chain_contrast",
    # report + errors
    "S3Report",
    "S3Error",
    "InputError",
    "ConsistencyError",
]

__version__ = "0.1.0"
