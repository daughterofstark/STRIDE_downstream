"""Canonical table construction tests (both tables kept separate)."""
from __future__ import annotations

from pathlib import Path

from stride_analysis._synthetic import (
    make_correlations_df,
    make_mechanism_dict,
    make_profile_df,
    make_real_correlations_df,
)
from stride_analysis.canonical import (
    assemble_replicate_table,
    assemble_stride_table,
    build_replicate_rows,
    build_stride_rows,
    replicate_canon_label,
)
from stride_analysis.models import MechanismFile, ReplicateInput
from stride_analysis.models.schema import (
    MECHANISM_MERGE_COLUMNS,
    STRIDE_TABLE_COLUMNS,
)


def _rep_input(sero: str, run: str, idx: int) -> ReplicateInput:
    return ReplicateInput(
        serotype=sero, run_dir=run, replicate_index=idx,
        correlations_path=Path(f"/x/{run}/{sero}/analysis_output/c.csv"),
    )


# -- canonical join key construction (Bug 2 regression) ---------------------
def test_replicate_canon_label_rebuilds_from_chain_and_resid() -> None:
    """Real STRIDE style: key is ``{chain}:{canon_resid}``, not the label."""
    import pandas as pd

    df = pd.DataFrame(
        {
            "canon_resid": [51, -24, 112],
            "name": ["HIS", "ALA", "GLY"],
            # resname-style labels that are NOT the canonical join key
            "label": ["HIS-51", "ALA-24", "GLY112"],
            "chain": ["NS3", "NS2B", "NS3"],
        }
    )
    assert replicate_canon_label(df).tolist() == ["NS3:51", "NS2B:-24", "NS3:112"]


def test_replicate_canon_label_falls_back_to_label_without_chain() -> None:
    """Synthetic style: no ``chain`` column, ``label`` already canonical."""
    import pandas as pd

    df = pd.DataFrame(
        {
            "canon_resid": [51, 200],
            "name": ["HIS", "LYS"],
            "label": ["NS3:51", "NS3:200"],
        }
    )
    assert replicate_canon_label(df).tolist() == ["NS3:51", "NS3:200"]


def test_build_replicate_rows_uses_canonical_key_for_real_style() -> None:
    """The real-style ``label`` (``HIS-51``) must not leak into ``canon_label``."""
    df = make_real_correlations_df("DENV1", 1)
    out = build_replicate_rows(df, _rep_input("DENV1", "1st_run", 1))
    assert set(out["canon_label"]) == {"NS3:51", "NS3:200"}
    assert not set(out["canon_label"]) & set(out["label"])  # disjoint from label


# -- replicate table --------------------------------------------------------
def test_replicate_rows_have_identity_and_provenance() -> None:
    df = make_correlations_df("DENV1", 1)
    out = build_replicate_rows(df, _rep_input("DENV1", "1st_run", 1))
    assert list(out.columns[:4]) == [
        "serotype", "replicate", "replicate_index", "canon_label",
    ]
    assert (out["serotype"] == "DENV1").all()
    assert (out["replicate_index"] == 1).all()
    assert "source_path" in out.columns


def test_replicate_table_key_unique() -> None:
    frames = [
        build_replicate_rows(make_correlations_df("DENV1", i), _rep_input("DENV1", f"{i}_run", i))
        for i in (1, 2, 3)
    ]
    tab = assemble_replicate_table(frames)
    assert not tab.duplicated(["serotype", "replicate", "canon_label"]).any()
    assert len(tab) == 6  # 3 replicates x 2 residues


def test_replicate_table_unions_columns() -> None:
    # one replicate has an extra column; union-align, no crash
    a = build_replicate_rows(make_correlations_df("DENV1", 1), _rep_input("DENV1", "1st_run", 1))
    b_df = make_correlations_df("DENV1", 2)
    b_df["future_col"] = 7.0
    b = build_replicate_rows(b_df, _rep_input("DENV1", "2nd_run", 2))
    tab = assemble_replicate_table([a, b])
    assert "future_col" in tab.columns
    assert tab["future_col"].isna().sum() == 2  # only b's rows have it


# -- stride table -----------------------------------------------------------
def test_stride_rows_shape_and_columns() -> None:
    prof = make_profile_df("DENV1")
    mech = MechanismFile.model_validate(make_mechanism_dict("DENV1"))
    out = build_stride_rows(prof, mech, "DENV1", "p.csv", "m.json")
    assert list(out.columns) == list(STRIDE_TABLE_COLUMNS)
    assert len(out) == 14  # 2 loci x 7 scales


def test_stride_mech_payload_only_on_gated_rows() -> None:
    prof = make_profile_df("DENV1")
    mech = MechanismFile.model_validate(make_mechanism_dict("DENV1"))
    out = build_stride_rows(prof, mech, "DENV1", "p.csv", "m.json")
    gated = out[out["is_gated_scale"]]
    nongated = out[~out["is_gated_scale"]]
    assert len(gated) == 2
    assert gated["mech_direction"].notna().all()
    for col in MECHANISM_MERGE_COLUMNS:
        assert nongated[col].isna().all()


def test_stride_table_key_unique() -> None:
    frames = []
    for sero in ("DENV1", "DENV2"):
        prof = make_profile_df(sero)
        mech = MechanismFile.model_validate(make_mechanism_dict(sero))
        frames.append(build_stride_rows(prof, mech, sero, "p", "m"))
    tab = assemble_stride_table(frames)
    assert not tab.duplicated(["serotype", "canon_label", "scale_level"]).any()
    assert len(tab) == 28


def test_tables_are_separate_no_shared_columns_beyond_keys() -> None:
    # the two canonical tables must not be collapsed: verify their schemas differ
    rep = build_replicate_rows(make_correlations_df("DENV1", 1), _rep_input("DENV1", "1st_run", 1))
    prof = make_profile_df("DENV1")
    mech = MechanismFile.model_validate(make_mechanism_dict("DENV1"))
    stride = build_stride_rows(prof, mech, "DENV1", "p", "m")
    # replicate table has 'r'/'abs_r' effect columns; stride table has 'rho'/'beta'
    assert "r" in rep.columns and "rho" not in rep.columns
    assert "rho" in stride.columns and "abs_r" not in stride.columns
