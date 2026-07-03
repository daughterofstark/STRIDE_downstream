"""Shared pytest fixtures for S3 tests."""
from __future__ import annotations

import pandas as pd
import pytest

from tests.s3.fixtures import make_stride_table


@pytest.fixture
def stride_table() -> pd.DataFrame:
    return make_stride_table()
