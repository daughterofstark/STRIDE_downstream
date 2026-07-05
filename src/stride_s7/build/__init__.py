"""S7 build subpackage: deterministic figure-data prep and manuscript-table assembly."""
from __future__ import annotations

from .figures import (
    FIGURE_BUILDERS,
    build_all_figures,
    prepare_f1,
    prepare_f2,
    prepare_f3,
    prepare_f4,
    prepare_f5,
    prepare_f6,
    prepare_f7,
    prepare_f8,
)
from .tables import (
    build_all_tables,
    build_t1,
    build_t2,
    build_t3,
    build_t4,
    build_t5,
)

__all__ = [
    "FIGURE_BUILDERS",
    "build_all_figures",
    "prepare_f1",
    "prepare_f2",
    "prepare_f3",
    "prepare_f4",
    "prepare_f5",
    "prepare_f6",
    "prepare_f7",
    "prepare_f8",
    "build_all_tables",
    "build_t1",
    "build_t2",
    "build_t3",
    "build_t4",
    "build_t5",
]
