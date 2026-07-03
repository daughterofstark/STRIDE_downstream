"""Shared pytest fixtures for S1B tests."""
from __future__ import annotations

import pandas as pd
import pytest

from tests.s1b.fixtures import (
    make_canonical_residues,
    make_conservation_table,
    make_domain_table,
    make_replicate_inventory,
)


@pytest.fixture
def canonical_residues() -> pd.DataFrame:
    return make_canonical_residues()


@pytest.fixture
def domain_table() -> pd.DataFrame:
    return make_domain_table()


@pytest.fixture
def conservation_table() -> pd.DataFrame:
    return make_conservation_table()


@pytest.fixture
def replicate_inventory() -> pd.DataFrame:
    return make_replicate_inventory()
