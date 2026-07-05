r"""Stage S7 — the reporting layer (publication figures & manuscript tables).

S7 is the terminal, reporting-only stage. It consumes the already-reduced outputs
of S2–S6 and assembles the design's publication figures (F1–F8) and manuscript
tables (T1–T5), deterministically and provenance-stamped. It performs **no new
statistics and no biological inference**, never reads the raw STRIDE files, and
never recomputes a quantity a prior stage produced.

Figure → source (all fully supported by existing outputs):

- F1 reproducibility landscape ← S2 ``residue_landscape``
- F2 achieved-resolution census ← S2 ``resolution_census``
- F3 domain × serotype ρ heatmap ← S5 ``domain_serotype_matrix``
- F4 signed-effect forest ← S4 ``significance_screen``
- F5 cross-serotype conservation ← S5 ``position_conservation``
- F6 variance composition ← S4 ``variance_budget``
- F7 ρ-vs-scale (catalytic) ← S3 ``scale_curve``
- F8 coherence vs ρ ← S2 ``domain_reproducibility``

Table → source:

- T1 per-serotype summary ← S2 ``serotype_summary``
- T2 domain ρ + signed effect ← S2 ``domain_reproducibility`` + S4 ``domain_effect_summary``
- T3 catalytic cross-serotype ← S5 ``domain_serotype_matrix`` (catalytic)
- T4 top shared signed positions ← S5 ``position_conservation``
- T5 variance-component budget ← S4 ``variance_budget``

The replicate layer (S6) feeds no design figure/table; its blocked-analysis ledger
is surfaced under the report's ``limitations`` and its inputs digested for
provenance. Figures are emitted as deterministic, dependency-free SVG (plus
prepared-data CSV/Parquet); tables as CSV/Parquet/Markdown. See ``docs/s7.md``.
"""
from __future__ import annotations

from .build import build_all_figures, build_all_tables
from .models import ArtifactRecord, S7Report, ValidationCheck
from .models.errors import ConsistencyError, InputError, S7Error
from .s7 import DEFAULT_STAGE_DIRS, S7Artifacts, build_s7, run_s7

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "build_s7",
    "run_s7",
    "S7Artifacts",
    "DEFAULT_STAGE_DIRS",
    "S7Report",
    "ArtifactRecord",
    "ValidationCheck",
    "S7Error",
    "InputError",
    "ConsistencyError",
    "build_all_figures",
    "build_all_tables",
]
