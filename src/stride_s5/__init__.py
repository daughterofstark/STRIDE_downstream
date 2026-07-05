"""stride_s5 — Stage S5: the cross-serotype layer (n = 4).

S5 consumes the S0 STRIDE table (``stride_table.parquet``) — the tidy profile —
and the S1A ``conservation_table.parquet`` — the shared-position index — and
compares the four dengue serotypes at the level of shared canonical positions and
named regions:

- ``position_conservation.parquet``    — conservation of reproducibility across
  the shared positions (all / majority / some / none), with serotype-divergent
  and Catalytic-Triad flags (Tier B)
- ``direction_concordance.parquet``    — do serotypes agree on the sign of shared
  signed positions? (agree / majority / conflict; Tier B)
- ``domain_serotype_matrix.parquet``   — the tidy-long ρ(domain × serotype) matrix
  over the NS3 domains + NS2B, catalytic domains flagged (Tier A)
- ``cross_serotype_scorecard.parquet`` — per-serotype scorecard: n_loci,
  %residue-gated, %signed, %mixed, ρ median, shared-position counts (Tier B)
- ``conservation_summary.json``        — facts, provenance header, validation
  outcomes

S5 respects the design's cross-serotype guardrails: **serotype is the unit of
replication (n = 4)** — per-serotype values are aggregated first, then compared,
and residues are never treated as independent samples; results are **descriptive**
(counts / effect sizes), not p-values across residues. The gate is **uncalibrated**
— every reproducibility statement is relative to a provisional ρ\\* = 0.5 and
labelled provisional. **Domain scale is the licensed claim level** at K = 3; the
position and scorecard products are residue-scale and labelled exploratory. S5
reads ρ and the variance components (never recomputes them), produces no figures,
and never modifies S0, S1A, S1B, S2, S3, or S4.

Public API
----------
- :func:`run_s5`   — load → build → validate → write artifacts.
- :func:`build_s5` — same, in memory (no writes).
- :class:`S5Tables`— the four returned tables.
- :class:`S5Error` — base of every S5 exception.
"""
from __future__ import annotations

from .build import (
    build_cross_serotype_scorecard,
    build_direction_concordance,
    build_domain_serotype_matrix,
    build_position_conservation,
)
from .models import S5Report
from .models.errors import ConsistencyError, InputError, S5Error
from .s5 import S5Tables, build_s5, run_s5

__all__ = [
    # orchestration
    "run_s5",
    "build_s5",
    "S5Tables",
    # builders (reusable by later stages)
    "build_position_conservation",
    "build_direction_concordance",
    "build_domain_serotype_matrix",
    "build_cross_serotype_scorecard",
    # report + errors
    "S5Report",
    "S5Error",
    "InputError",
    "ConsistencyError",
]

__version__ = "0.1.0"
