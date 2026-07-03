"""Shared pytest fixtures for S2 tests."""
from __future__ import annotations

import pandas as pd
import pytest

from tests.s2.fixtures import (
    make_domain_annotation,
    make_residue_annotation,
    make_stride_table,
)


@pytest.fixture
def stride_table() -> pd.DataFrame:
    return make_stride_table()


@pytest.fixture
def residue_annotation() -> pd.DataFrame:
    return make_residue_annotation()


@pytest.fixture
def domain_annotation() -> pd.DataFrame:
    return make_domain_annotation()
