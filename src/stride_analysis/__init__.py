"""stride_analysis — downstream analysis framework for STRIDE outputs.

This package is an *analysis framework*: real datasets are user-supplied inputs,
not part of the source. It currently implements **Stage S0** (ingestion,
validation, and canonical-table construction) as reusable subpackages:

- :mod:`stride_analysis.io`          — dataset discovery + loaders
- :mod:`stride_analysis.validation`  — schema + consistency checks
- :mod:`stride_analysis.canonical`   — canonical-table builders
- :mod:`stride_analysis.models`      — schemas, pydantic models, errors

Later stages import these subpackages directly; there is no monolithic pipeline.

Top-level convenience API
-------------------------
- :func:`run_s0`      — discover → validate → build → write artifacts.
- :func:`build_tables`— same, in memory (no writes).
"""
from __future__ import annotations

from .io import discover_dataset
from .models import Dataset, MechanismFile, Report, SerotypeDataset
from .models.errors import (
    ConsistencyError,
    DiscoveryError,
    HierarchyError,
    SchemaError,
    StrideAnalysisError,
)
from .s0 import build_tables, run_s0

__all__ = [
    # orchestration
    "run_s0",
    "build_tables",
    "discover_dataset",
    # models
    "Dataset",
    "SerotypeDataset",
    "MechanismFile",
    "Report",
    # errors
    "StrideAnalysisError",
    "DiscoveryError",
    "SchemaError",
    "HierarchyError",
    "ConsistencyError",
]

__version__ = "0.2.0"
