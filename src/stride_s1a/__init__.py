"""stride_s1a — Stage S1A: the reusable biological data layer.

S1A consumes **only** the S0 canonical tables (``stride_table.parquet`` and
``replicate_table.parquet``) and builds reusable biological annotations and
derived datasets. It performs no statistics, no scoring, no clustering, no
interpretation, and produces no figures.

Outputs
-------
- ``canonical_residues.parquet``  — one canonical residue object per (serotype, canon_label)
- ``domain_table.parquet``        — structural domain summaries
- ``replicate_inventory.parquet`` — per-residue replicate availability
- ``conservation_table.parquet``  — cross-serotype presence map
- ``dataset_summary.json``        — machine-readable facts + validation outcomes

Public API
----------
- :func:`run_s1a`   — load → build → validate → write artifacts.
- :func:`build_s1a` — same, in memory (no writes).
- :class:`S1ATables`— the four returned tables.
- :class:`S1AError` — base of every S1A exception.
"""
from __future__ import annotations

from .build import (
    build_canonical_residues,
    build_conservation_table,
    build_domain_table,
    build_replicate_inventory,
)
from .models import S1AReport
from .models.errors import ConsistencyError, InputError, S1AError
from .s1a import S1ATables, build_s1a, run_s1a

__all__ = [
    # orchestration
    "run_s1a",
    "build_s1a",
    "S1ATables",
    # builders (reusable by later stages)
    "build_canonical_residues",
    "build_domain_table",
    "build_replicate_inventory",
    "build_conservation_table",
    # report + errors
    "S1AReport",
    "S1AError",
    "InputError",
    "ConsistencyError",
]

__version__ = "0.1.0"
