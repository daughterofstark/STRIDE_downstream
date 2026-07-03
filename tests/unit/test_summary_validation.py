"""Level-2 summary schema, consistency, and hierarchy tests."""
from __future__ import annotations

import copy

import pandas as pd
import pytest

from stride_analysis.models import MechanismFile
from stride_analysis.models.errors import (
    ConsistencyError,
    HierarchyError,
    SchemaError,
)
from stride_analysis.validation import (
    check_profile_mechanism_consistency,
    parse_residue_path,
    split_region_id,
    validate_path_depth,
    validate_profile_schema,
)


# ---------------------------------------------------------------------------
# profile schema
# ---------------------------------------------------------------------------
def test_valid_profile_passes(profile_df: pd.DataFrame) -> None:
    validate_profile_schema(profile_df, "DENV1")


def test_missing_column_raises(profile_df: pd.DataFrame) -> None:
    with pytest.raises(SchemaError, match="missing required columns.*rho"):
        validate_profile_schema(profile_df.drop(columns=["rho"]), "DENV1")


def test_extra_column_raises(profile_df: pd.DataFrame) -> None:
    df = profile_df.copy()
    df["surprise"] = 1
    with pytest.raises(SchemaError, match="extra columns"):
        validate_profile_schema(df, "DENV1")


def test_rho_out_of_range_raises(profile_df: pd.DataFrame) -> None:
    df = profile_df.copy()
    df.loc[0, "rho"] = 2.0
    with pytest.raises(SchemaError, match=r"outside \[0, 1\]"):
        validate_profile_schema(df, "DENV1")


def test_negative_variance_raises(profile_df: pd.DataFrame) -> None:
    df = profile_df.copy()
    df.loc[0, "tau2"] = -1.0
    with pytest.raises(SchemaError, match="tau2.*negative"):
        validate_profile_schema(df, "DENV1")


def test_scale_mismatch_raises(profile_df: pd.DataFrame) -> None:
    df = profile_df.copy()
    idx = df[df["scale_level"] == "residue"].index[0]
    df.loc[idx, "scale_level"] = "domain"
    with pytest.raises(SchemaError, match="implies"):
        validate_profile_schema(df, "DENV1")


def test_duplicate_locus_scale_raises(profile_df: pd.DataFrame) -> None:
    df = pd.concat([profile_df, profile_df.iloc[[0]]], ignore_index=True)
    with pytest.raises(SchemaError, match="duplicate \\(locus, scale_level\\)"):
        validate_profile_schema(df, "DENV1")


def test_protein_mismatch_raises(profile_df: pd.DataFrame) -> None:
    with pytest.raises(SchemaError, match="protein column"):
        validate_profile_schema(profile_df, "WRONG")


def test_residue_region_id_must_equal_locus(profile_df: pd.DataFrame) -> None:
    df = profile_df.copy()
    idx = df[df["scale_level"] == "residue"].index[0]
    df.loc[idx, "region_id"] = "CPLX/prot/NS3/Catalytic Triad/none/unknown/NS3:999"
    with pytest.raises(SchemaError, match="region_id != locus"):
        validate_profile_schema(df, "DENV1")


def test_multiple_gated_per_locus_raises(profile_df: pd.DataFrame) -> None:
    df = profile_df.copy()
    a_locus = df["locus"].iloc[0]
    dom = df[(df["locus"] == a_locus) & (df["scale_level"] == "domain")].index[0]
    df.loc[dom, "gated"] = True
    with pytest.raises(SchemaError, match="more than one|>1 gated") as exc:
        validate_profile_schema(df, "DENV1")
    # the message must name the offending locus, not leak a literal f-string
    assert a_locus in str(exc.value)
    assert "{k:" not in str(exc.value)


# ---------------------------------------------------------------------------
# mechanism (pydantic)
# ---------------------------------------------------------------------------
def test_mechanism_mixed_with_beta_raises(mechanism_dict: dict) -> None:
    bad = copy.deepcopy(mechanism_dict)
    bad["mechanisms"][1]["beta_signed"] = 0.5  # B is mixed
    with pytest.raises(Exception, match="mixed"):
        MechanismFile.model_validate(bad)


def test_mechanism_bad_direction_raises(mechanism_dict: dict) -> None:
    bad = copy.deepcopy(mechanism_dict)
    bad["mechanisms"][0]["direction"] = "sideways"
    with pytest.raises(Exception, match="direction"):
        MechanismFile.model_validate(bad)


def test_mechanism_summary_mismatch_raises(mechanism_dict: dict) -> None:
    bad = copy.deepcopy(mechanism_dict)
    bad["summary"]["n_mechanisms"] = 99
    with pytest.raises(Exception, match="n_mechanisms"):
        MechanismFile.model_validate(bad)


# ---------------------------------------------------------------------------
# consistency
# ---------------------------------------------------------------------------
def test_consistency_passes(profile_df: pd.DataFrame, mechanism_dict: dict) -> None:
    mech = MechanismFile.model_validate(mechanism_dict)
    check_profile_mechanism_consistency(profile_df, mech, "DENV1")


def test_rho_mismatch_raises(profile_df: pd.DataFrame, mechanism_dict: dict) -> None:
    bad = copy.deepcopy(mechanism_dict)
    bad["mechanisms"][0]["rho"] = 0.123
    mech = MechanismFile.model_validate(bad)
    with pytest.raises(ConsistencyError, match="gated rho.*mechanism rho"):
        check_profile_mechanism_consistency(profile_df, mech, "DENV1")


def test_orphan_mechanism_raises(profile_df: pd.DataFrame, mechanism_dict: dict) -> None:
    bad = copy.deepcopy(mechanism_dict)
    bad["mechanisms"][0]["loci"] = ["CPLX/prot/NS3/x/none/unknown/NS3:999"]
    mech = MechanismFile.model_validate(bad)
    with pytest.raises(ConsistencyError, match="orphan mechanism"):
        check_profile_mechanism_consistency(profile_df, mech, "DENV1")


def test_orphan_profile_row_raises(profile_df: pd.DataFrame, mechanism_dict: dict) -> None:
    bad = copy.deepcopy(mechanism_dict)
    bad["mechanisms"] = bad["mechanisms"][:1]
    bad["summary"]["n_mechanisms"] = 1
    bad["summary"]["n_gate_uncertain"] = 1
    mech = MechanismFile.model_validate(bad)
    with pytest.raises(ConsistencyError, match="orphan profile rows"):
        check_profile_mechanism_consistency(profile_df, mech, "DENV1")


# ---------------------------------------------------------------------------
# hierarchy
# ---------------------------------------------------------------------------
def test_parse_residue_path() -> None:
    parsed = parse_residue_path("CPLX/prot/NS3/Catalytic Triad/none/unknown/NS3:51")
    assert parsed["h_chain"] == "NS3"
    assert parsed["h_domain"] == "Catalytic Triad"


def test_split_double_slash_raises() -> None:
    with pytest.raises(HierarchyError, match="empty segment"):
        split_region_id("a//b")


def test_validate_path_depth_mismatch_raises() -> None:
    with pytest.raises(HierarchyError, match="expected 4"):
        validate_path_depth("CPLX/prot/NS3", "domain")


def test_parse_wrong_depth_raises() -> None:
    with pytest.raises(HierarchyError, match="expected 7"):
        parse_residue_path("CPLX/prot/NS3")
