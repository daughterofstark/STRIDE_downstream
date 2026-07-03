"""stride_s4 — Stage S4: the uncertainty layer.

S4 consumes the S0 STRIDE table (``stride_table.parquet``) — the tidy profile,
with every locus at every scale carrying the variance components τ² and σ̄² and
the gated mechanism payload — and derives per-domain and per-residue uncertainty
products:

- ``variance_budget.parquet``        — τ²/σ̄² budget per domain (Tier A)
- ``residue_variance.parquet``       — per-residue decomposition + τ² ranking (Tier B)
- ``significance_screen.parquet``    — CI-exclusion + Wald p + BH-FDR per mechanism
- ``domain_effect_summary.parquet``  — β_se-weighted effect summary per domain (Tier A)
- ``uncertainty_summary.json``       — facts, provenance header, validation outcomes

S4's products are **ρ*-independent** descriptions of the profile and the emitted
mechanisms: it makes no calibrated pass/fail resolution claim (the gate is
uncalibrated, §0.1), never recomputes ρ or the variance components, performs **no
cross-serotype tests** (that is S5 — the serotype is the unit of replication and
the FDR family), and produces no figures. It never modifies S0, S1A, S1B, S2, or
S3.

Two tiers are labelled explicitly per the design: **licensed** (domain-scale,
the claim level at K=3) and **exploratory** (residue-scale, outside the operating
range).

Public API
----------
- :func:`run_s4`   — load → build → validate → write artifacts.
- :func:`build_s4` — same, in memory (no writes).
- :class:`S4Tables`— the four returned tables.
- :class:`S4Error` — base of every S4 exception.
"""
from __future__ import annotations

from .build import (
    build_domain_effect_summary,
    build_residue_variance,
    build_significance_screen,
    build_variance_budget,
)
from .models import S4Report
from .models.errors import ConsistencyError, InputError, S4Error
from .s4 import S4Tables, build_s4, run_s4

__all__ = [
    # orchestration
    "run_s4",
    "build_s4",
    "S4Tables",
    # builders (reusable by later stages)
    "build_variance_budget",
    "build_residue_variance",
    "build_significance_screen",
    "build_domain_effect_summary",
    # report + errors
    "S4Report",
    "S4Error",
    "InputError",
    "ConsistencyError",
]

__version__ = "0.1.0"
