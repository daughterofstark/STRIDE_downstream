"""Shared pytest fixtures for S6 tests."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_s6.build import per_run_effects
from tests.s6.fixtures import (
    make_empty_inventory,
    make_replicate_inventory,
    make_replicate_table,
)


@pytest.fixture
def replicate_inventory() -> pd.DataFrame:
    return make_replicate_inventory()


@pytest.fixture
def replicate_table() -> pd.DataFrame:
    return make_replicate_table()


@pytest.fixture
def effects(replicate_table: pd.DataFrame) -> pd.DataFrame:
    return per_run_effects(replicate_table)


@pytest.fixture
def empty_inventory() -> pd.DataFrame:
    return make_empty_inventory()
