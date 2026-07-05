"""Loader, orchestration, and CLI tests for S5."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from stride_s5 import build_s5, run_s5
from stride_s5.__main__ import main
from stride_s5.io import (
    file_digest,
    load_conservation_table,
    load_stride_table,
)
from stride_s5.models.errors import InputError
from stride_s5.models.schema import (
    CONSERVATION_TABLE_REQUIRED,
    OUT_CONSERVATION_SUMMARY,
    OUT_CROSS_SEROTYPE_SCORECARD,
    OUT_DIRECTION_CONCORDANCE,
    OUT_DOMAIN_SEROTYPE_MATRIX,
    OUT_POSITION_CONSERVATION,
    STRIDE_TABLE_REQUIRED,
)
from tests.s5.fixtures import write_inputs


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------
def test_load_missing_stride_raises(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="stride_table not found"):
        load_stride_table(tmp_path / "nope.parquet")


def test_load_missing_conservation_raises(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="conservation_table not found"):
        load_conservation_table(tmp_path / "nope.parquet")


def test_load_missing_columns_raises(tmp_path: Path) -> None:
    p = tmp_path / "stride_table.parquet"
    pd.DataFrame({"serotype": ["DENV1"]}).to_parquet(p, index=False)
    with pytest.raises(InputError, match="missing required column"):
        load_stride_table(p)


def test_load_conservation_missing_columns_raises(tmp_path: Path) -> None:
    p = tmp_path / "conservation_table.parquet"
    pd.DataFrame({"canon_label": ["NS3:51"]}).to_parquet(p, index=False)
    with pytest.raises(InputError, match="missing required column"):
        load_conservation_table(p)


def test_load_unreadable_raises(tmp_path: Path) -> None:
    p = tmp_path / "stride_table.parquet"
    p.write_bytes(b"not a parquet file")
    with pytest.raises(InputError, match="could not read stride_table"):
        load_stride_table(p)


def test_load_roundtrip(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    st = load_stride_table(paths["stride_table"])
    ct = load_conservation_table(paths["conservation_table"])
    for col in STRIDE_TABLE_REQUIRED:
        assert col in st.columns
    for col in CONSERVATION_TABLE_REQUIRED:
        assert col in ct.columns


def test_file_digest_stable_and_missing(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    d1 = file_digest(paths["stride_table"])
    d2 = file_digest(paths["stride_table"])
    assert d1 == d2 and len(d1) == 64
    assert file_digest(tmp_path / "absent.parquet") == ""


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def test_build_s5_end_to_end(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    tables, report = build_s5(
        paths["stride_table"], paths["conservation_table"]
    )
    assert report.all_passed
    assert report.serotypes == ["DENV1", "DENV2", "DENV3", "DENV4"]
    assert len(tables.position_conservation) == 6
    assert len(tables.direction_concordance) == 4
    assert len(tables.domain_serotype_matrix) == 15
    assert len(tables.cross_serotype_scorecard) == 4
    # provenance: uncalibrated, provisional ρ*, both input digests
    assert report.provenance["calibrated"] is False
    assert report.provenance["provisional_rho_star"] == 0.5
    assert report.provenance["inputs"]["stride_table"]["sha256"]
    assert report.provenance["inputs"]["conservation_table"]["sha256"]
    # facts capture the cross-serotype summary
    assert report.facts["n_serotype_divergent"] == 2
    assert report.facts["n_catalytic_triad_positions"] == 3
    assert report.facts["conservation_class_counts"]["reproducible_all"] == 3
    assert report.facts["concordance_class_counts"]["conflict"] == 1


def test_run_s5_writes_all_artifacts(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    out = tmp_path / "out"
    run_s5(paths["stride_table"], paths["conservation_table"], out)
    for name in (
        OUT_POSITION_CONSERVATION,
        OUT_DIRECTION_CONCORDANCE,
        OUT_DOMAIN_SEROTYPE_MATRIX,
        OUT_CROSS_SEROTYPE_SCORECARD,
        OUT_CONSERVATION_SUMMARY,
    ):
        assert (out / name).exists(), f"missing {name}"
    payload = json.loads((out / OUT_CONSERVATION_SUMMARY).read_text())
    assert payload["stage"] == "S5"
    assert payload["all_checks_passed"] is True
    assert payload["provenance"]["calibrated"] is False


def test_run_s5_deterministic(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    t1, _ = run_s5(
        paths["stride_table"], paths["conservation_table"], tmp_path / "o1"
    )
    t2, _ = run_s5(
        paths["stride_table"], paths["conservation_table"], tmp_path / "o2"
    )
    pd.testing.assert_frame_equal(
        t1.position_conservation, t2.position_conservation
    )
    pd.testing.assert_frame_equal(
        t1.direction_concordance, t2.direction_concordance
    )
    pd.testing.assert_frame_equal(
        t1.domain_serotype_matrix, t2.domain_serotype_matrix
    )
    pd.testing.assert_frame_equal(
        t1.cross_serotype_scorecard, t2.cross_serotype_scorecard
    )


def test_run_s5_written_parquet_matches(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    out = tmp_path / "out"
    tables, _ = run_s5(paths["stride_table"], paths["conservation_table"], out)
    reread = pd.read_parquet(out / OUT_POSITION_CONSERVATION)
    pd.testing.assert_frame_equal(reread, tables.position_conservation)


def test_build_s5_empty_dataset(tmp_path: Path) -> None:
    empty_st = pd.DataFrame(columns=list(STRIDE_TABLE_REQUIRED))
    empty_ct = pd.DataFrame(columns=list(CONSERVATION_TABLE_REQUIRED))
    paths = write_inputs(
        tmp_path, stride_table=empty_st, conservation_table=empty_ct
    )
    tables, report = build_s5(
        paths["stride_table"], paths["conservation_table"]
    )
    assert report.all_passed
    assert tables.position_conservation.empty
    assert tables.cross_serotype_scorecard.empty
    assert report.facts == {}


def test_rho_star_override_changes_reproducibility(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    # at ρ*=0.65, NS3:135 (min ρ 0.60) drops out of "all" reproducible
    tables, _ = build_s5(
        paths["stride_table"], paths["conservation_table"], rho_star=0.65
    )
    ns135 = tables.position_conservation[
        tables.position_conservation["canon_label"] == "NS3:135"
    ].iloc[0]
    assert ns135["conservation_class"] != "reproducible_all"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def test_cli_success(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    write_inputs(tmp_path / "in")
    rc = main(
        [
            "--input-dir",
            str(tmp_path / "in"),
            "--conservation-input-dir",
            str(tmp_path / "in"),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 0
    assert "S5 OK" in capsys.readouterr().out
    assert (tmp_path / "out" / OUT_POSITION_CONSERVATION).exists()


def test_cli_explicit_paths(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    paths = write_inputs(tmp_path / "in")
    rc = main(
        [
            "--stride-table",
            str(paths["stride_table"]),
            "--conservation-table",
            str(paths["conservation_table"]),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 0
    assert (tmp_path / "out" / OUT_DIRECTION_CONCORDANCE).exists()


def test_cli_failure_returns_1(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc = main(
        [
            "--input-dir",
            str(tmp_path / "missing"),
            "--conservation-input-dir",
            str(tmp_path / "missing"),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 1
    assert "S5 FAILED" in capsys.readouterr().err
