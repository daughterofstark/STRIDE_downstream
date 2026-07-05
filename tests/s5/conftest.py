"""Shared pytest fixtures for S5 tests."""
from __future__ import annotations

import pandas as pd
import pytest

from tests.s5.fixtures import make_conservation_table, make_stride_table


@pytest.fixture
def stride_table() -> pd.DataFrame:
    return make_stride_table()


@pytest.fixture
def conservation_table() -> pd.DataFrame:
    return make_conservation_table()
