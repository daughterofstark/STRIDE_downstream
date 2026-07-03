"""Loader tests: malformed files at the IO boundary."""
from __future__ import annotations

from pathlib import Path

import pytest

from stride_analysis.io import load_correlations, load_mechanism, load_profile
from stride_analysis.models.errors import SchemaError


def test_unparseable_json_raises(tmp_path: Path) -> None:
    p = tmp_path / "X_mechanism.json"
    p.write_text("{ not json ")
    with pytest.raises(SchemaError, match="could not parse mechanism JSON"):
        load_mechanism(p)


def test_json_missing_fields_raises(tmp_path: Path) -> None:
    p = tmp_path / "X_mechanism.json"
    p.write_text('{"schema_version": "m5"}')
    with pytest.raises(SchemaError, match="failed schema validation"):
        load_mechanism(p)


def test_empty_profile_raises(tmp_path: Path) -> None:
    p = tmp_path / "X_profile.csv"
    p.write_text("protein,locus\n")
    with pytest.raises(SchemaError, match="empty"):
        load_profile(p)


def test_empty_correlations_raises(tmp_path: Path) -> None:
    p = tmp_path / "X_correlations_v5.csv"
    p.write_text("file_resid,label\n")
    with pytest.raises(SchemaError, match="empty"):
        load_correlations(p)


def test_binary_csv_raises(tmp_path: Path) -> None:
    p = tmp_path / "X_correlations_v5.csv"
    p.write_bytes(b"\x00\x01\x02\xff\xfe")
    with pytest.raises(SchemaError):
        load_correlations(p)
