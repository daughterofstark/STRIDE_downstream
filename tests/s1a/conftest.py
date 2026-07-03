"""Shared pytest fixtures for S1A tests."""
from __future__ import annotations

import pandas as pd
import pytest

from tests.s1a.fixtures import make_replicate_table, make_stride_table


@pytest.fixture
def stride_table() -> pd.DataFrame:
    return make_stride_table()


@pytest.fixture
def replicate_table() -> pd.DataFrame:
    return make_replicate_table()
