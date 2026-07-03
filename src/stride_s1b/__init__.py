"""stride_s1b — Stage S1B: the biological annotation layer.

S1B consumes **only** the four S1A parquet tables (``canonical_residues``,
``domain_table``, ``replicate_inventory``, ``conservation_table``) and builds
reusable, deterministic biological annotation tables. It never reads raw STRIDE
outputs, ``profile.csv``/``mechanism.json``, or MD trajectories, and it performs
no statistics, no ranking, no clustering, no hypothesis generation, and produces
no figures.

Outputs
-------
- ``residue_annotation.parquet``  — per-residue structural annotation
- ``domain_annotation.parquet``   — per-domain structural annotation
- ``hierarchy_annotation.parquet``— per-residue hierarchy resolution annotation
- ``serotype_annotation.parquet`` — per-serotype structural composition
- ``annotation_summary.json``     — machine-readable facts + validation outcomes

Public API
----------
- :func:`run_s1b`   — load → build → validate → write artifacts.
- :func:`build_s1b` — same, in memory (no writes).
- :class:`S1BTables`— the four returned tables.
- :class:`S1BError` — base of every S1B exception.
"""
from __future__ import annotations

from .build import (
    build_domain_annotation,
    build_hierarchy_annotation,
    build_residue_annotation,
    build_serotype_annotation,
)
from .models import S1BReport
from .models.errors import ConsistencyError, InputError, S1BError
from .s1b import S1BTables, build_s1b, run_s1b

__all__ = [
    # orchestration
    "run_s1b",
    "build_s1b",
    "S1BTables",
    # builders (reusable by later stages)
    "build_residue_annotation",
    "build_domain_annotation",
    "build_hierarchy_annotation",
    "build_serotype_annotation",
    # report + errors
    "S1BReport",
    "S1BError",
    "InputError",
    "ConsistencyError",
]

__version__ = "0.1.0"
