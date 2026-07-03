"""Level-1 replicate schema validation tests."""
from __future__ import annotations

import pandas as pd
import pytest

from stride_analysis.models.errors import SchemaError
from stride_analysis.validation import validate_correlations_schema


def test_valid_correlations_passes(correlations_df: pd.DataFrame) -> None:
    extra = validate_correlations_schema(correlations_df, "DENV1", "1st_run")
    assert extra == []  # synthetic table uses only known columns


def test_missing_required_column_raises(correlations_df: pd.DataFrame) -> None:
    df = correlations_df.drop(columns=["r"])
    with pytest.raises(SchemaError, match="missing required columns.*r"):
        validate_correlations_schema(df, "DENV1", "1st_run")


def test_missing_label_raises(correlations_df: pd.DataFrame) -> None:
    df = correlations_df.copy()
    df.loc[0, "label"] = None
    with pytest.raises(SchemaError, match="label.*missing"):
        validate_correlations_schema(df, "DENV1", "1st_run")


def test_blank_label_raises(correlations_df: pd.DataFrame) -> None:
    df = correlations_df.copy()
    df.loc[0, "label"] = "   "
    with pytest.raises(SchemaError, match="label.*blank"):
        validate_correlations_schema(df, "DENV1", "1st_run")


def test_nonfinite_effect_raises(correlations_df: pd.DataFrame) -> None:
    df = correlations_df.copy()
    df.loc[0, "r"] = float("inf")
    with pytest.raises(SchemaError, match="effect column.*non-finite"):
        validate_correlations_schema(df, "DENV1", "1st_run")


def test_non_numeric_effect_raises(correlations_df: pd.DataFrame) -> None:
    df = correlations_df.copy()
    df["abs_r"] = df["abs_r"].astype(object)
    df.loc[0, "abs_r"] = "oops"
    with pytest.raises(SchemaError, match="non-numeric|non-finite"):
        validate_correlations_schema(df, "DENV1", "1st_run")


def test_duplicate_residue_label_raises(correlations_df: pd.DataFrame) -> None:
    df = pd.concat([correlations_df, correlations_df.iloc[[0]]], ignore_index=True)
    with pytest.raises(SchemaError, match="duplicate residue label"):
        validate_correlations_schema(df, "DENV1", "1st_run")


def test_unknown_extra_columns_are_reported_not_rejected(
    correlations_df: pd.DataFrame,
) -> None:
    df = correlations_df.copy()
    df["some_future_stride_column"] = 1.0
    extra = validate_correlations_schema(df, "DENV1", "1st_run")
    assert extra == ["some_future_stride_column"]


def test_bad_int_column_raises(correlations_df: pd.DataFrame) -> None:
    df = correlations_df.copy()
    df["file_resid"] = df["file_resid"].astype(object)
    df.loc[0, "file_resid"] = 3.5
    with pytest.raises(SchemaError, match="integer-valued"):
        validate_correlations_schema(df, "DENV1", "1st_run")
