"""stride_s2 — Stage S2: the per-serotype reduction layer.

S2 consumes the S0 STRIDE table (``stride_table.parquet``) and the two S1B
annotation tables (``residue_annotation``, ``domain_annotation``) and derives
per-serotype reproducibility summaries: the achieved-resolution census, the
residue-scale reproducibility landscape, domain-scale reproducibility, the
signed/significant screen, and a per-serotype scorecard. Everything is reported
over a **ρ\\* band**, never at a single threshold, because the gate is
uncalibrated (design §0.1, §5.3). Two tiers are labelled explicitly:
**licensed** (domain-scale, the claim level at K=3) and **exploratory**
(residue-scale, outside the operating range).

S2 performs **no cross-serotype tests** (that is S5), makes no calibrated
pass/fail resolution claims, never recomputes ρ from its variance components,
and produces no figures. It never modifies S0, S1A, or S1B.

Outputs
-------
- ``resolution_census.parquet``      — per (serotype, ρ*, gated scale) locus counts
- ``residue_landscape.parquet``      — per-residue ρ landscape (Tier B, exploratory)
- ``domain_reproducibility.parquet`` — per-domain ρ (Tier A, licensed)
- ``signed_screen.parquet``          — signed/significant screen over the band
- ``serotype_summary.parquet``       — per-serotype scorecard over the band
- ``reduction_summary.json``         — facts, provenance header, validation outcomes

Public API
----------
- :func:`run_s2`   — load → build → validate → write artifacts.
- :func:`build_s2` — same, in memory (no writes).
- :class:`S2Tables`— the five returned tables.
- :class:`S2Error` — base of every S2 exception.
"""
from __future__ import annotations

from .build import (
    build_domain_reproducibility,
    build_residue_landscape,
    build_resolution_census,
    build_serotype_summary,
    build_signed_screen,
)
from .models import S2Report
from .models.errors import ConfigError, ConsistencyError, InputError, S2Error
from .s2 import S2Tables, build_s2, run_s2

__all__ = [
    # orchestration
    "run_s2",
    "build_s2",
    "S2Tables",
    # builders (reusable by later stages)
    "build_resolution_census",
    "build_residue_landscape",
    "build_domain_reproducibility",
    "build_signed_screen",
    "build_serotype_summary",
    # report + errors
    "S2Report",
    "S2Error",
    "InputError",
    "ConfigError",
    "ConsistencyError",
]

__version__ = "0.1.0"
