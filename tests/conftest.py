"""Shared pytest fixtures.

All fixtures build synthetic datasets on disk via the framework's own synthetic
generator. No test depends on real data.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from stride_analysis._synthetic import (
    make_correlations_df,
    make_mechanism_dict,
    make_profile_df,
    write_dataset,
)


@pytest.fixture
def serotypes() -> list[str]:
    return ["DENV1", "DENV2"]


@pytest.fixture
def run_dirs() -> list[str]:
    return ["1st_run", "2nd_run", "3rd_run"]


@pytest.fixture
def dataset_root(
    tmp_path: Path, serotypes: list[str], run_dirs: list[str]
) -> Path:
    """A complete, valid dataset (both levels) on disk."""
    return write_dataset(tmp_path / "data", serotypes, run_dirs)


@pytest.fixture
def correlations_df() -> pd.DataFrame:
    return make_correlations_df("DENV1", 1)


@pytest.fixture
def profile_df() -> pd.DataFrame:
    return make_profile_df("DENV1")


@pytest.fixture
def mechanism_dict() -> dict:
    return make_mechanism_dict("DENV1")
