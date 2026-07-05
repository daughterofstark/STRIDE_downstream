"""Shared pytest fixtures for S7 tests (synthetic inputs only)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from tests.s7.fixtures import make_all_inputs, write_empty_inputs, write_inputs


@pytest.fixture
def inputs() -> dict[str, pd.DataFrame]:
    """The full synthetic S2–S6 input frames keyed by filename."""
    return make_all_inputs()


@pytest.fixture
def stage_dirs(tmp_path: Path) -> dict[str, Path]:
    """Full synthetic inputs written to disk; returns the per-stage directories."""
    return write_inputs(tmp_path / "in")


@pytest.fixture
def empty_stage_dirs(tmp_path: Path) -> dict[str, Path]:
    """Present-but-empty inputs written to disk; returns the per-stage directories."""
    return write_empty_inputs(tmp_path / "empty_in")
