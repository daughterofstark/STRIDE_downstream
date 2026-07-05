"""stride_s6 — Stage S6: the replicate layer.

S6 is the **per-replicate-observation** stage. It consumes the S1A
``replicate_inventory.parquet`` (the replicate-structure index — how many runs
contain each ``canon_label`` per serotype, and whether the replicate set is
complete) and, when it exists, the S0 ``replicate_table.parquet`` (the Level-1
per-run observations — one row per ``(serotype, replicate, canon_label)`` carrying
that run's effect ``r``, i.e. the per-run θ). From these it derives the
replicate-axis products of design §3.1:

- ``replicate_regime.parquet``           — per-serotype replicate structure: K,
  completeness, residue-claim licensing (K>=5), per-run-effect availability
- ``replicate_effect_spread.parquet``    — descriptive across-run spread of each
  position's per-run θ (Tier B, exploratory; empty when θ is unavailable)
- ``replicate_concordance.parquet``      — per-serotype Kendall's W + mean pairwise
  Spearman of the per-run effect rankings (Tier B; empty when θ is unavailable)
- ``replicate_blocked_analyses.parquet`` — the ledger of blocked/unavailable
  replicate-level analyses, with the reason and required (absent) input for each
- ``replicate_summary.json``             — facts, provenance header, the blocked
  ledger, and validation outcomes

S6 is the only stage that reads the replicate inputs; it never reads
``stride_table`` and never re-implements the **τ²-based** replicate-disagreement
mapping or τ²/σ̄² regime diagnostic of §3.1 — those read the *aggregate* variance
components and are produced by **S4** (``residue_variance`` / ``variance_budget``).
S6 owns the *per-replicate* axis instead.

Two design-anticipated analyses are **blocked**, documented, and never
approximated: the **per-run** rank concordance / effect spread are blocked
whenever the per-run correlation CSVs are absent (the real dengue upload's state,
§4.1), and **leave-one-replicate-out** stability of ρ / the gated scale is always
blocked (it requires a STRIDE re-run, §4.2). Per-run θ is read, never recomputed,
and is not reconstructible from the K-aggregate (§4.2). Residue-scale per-run
products are exploratory at K = 3; the gate is uncalibrated (``calibrated=false``).
S6 never modifies S0, S1A, S1B, S2, S3, S4, or S5, and produces no figures.

Public API
----------
- :func:`run_s6`   — load → build → validate → write artifacts.
- :func:`build_s6` — same, in memory (no writes).
- :class:`S6Tables`— the four returned tables.
- :class:`S6Error` — base of every S6 exception.
"""
from __future__ import annotations

from .build import (
    build_replicate_blocked_analyses,
    build_replicate_concordance,
    build_replicate_effect_spread,
    build_replicate_regime,
)
from .models import S6Report
from .models.errors import ConsistencyError, InputError, S6Error
from .s6 import S6Tables, build_s6, run_s6

__version__ = "0.1.0"

__all__ = [
    # orchestration
    "run_s6",
    "build_s6",
    "S6Tables",
    # builders (reusable by later stages)
    "build_replicate_regime",
    "build_replicate_effect_spread",
    "build_replicate_concordance",
    "build_replicate_blocked_analyses",
    # report + errors
    "S6Report",
    "S6Error",
    "InputError",
    "ConsistencyError",
    # version
    "__version__",
]
