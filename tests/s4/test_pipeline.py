"""Loader, orchestration, and CLI tests for S4."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from stride_s4 import build_s4, run_s4
from stride_s4.__main__ import main
from stride_s4.io import file_digest, load_stride_table
from stride_s4.models.errors import InputError
from stride_s4.models.schema import (
    OUT_DOMAIN_EFFECT_SUMMARY,
    OUT_RESIDUE_VARIANCE,
    OUT_SIGNIFICANCE_SCREEN,
    OUT_UNCERTAINTY_SUMMARY,
    OUT_VARIANCE_BUDGET,
    STRIDE_TABLE_REQUIRED,
)
from tests.s4.fixtures import write_inputs


# ---------------------------------------------------------------------------
# loader
# ---------------------------------------------------------------------------
def test_load_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="stride_table not found"):
        load_stride_table(tmp_path / "nope.parquet")


def test_load_missing_columns_raises(tmp_path: Path) -> None:
    p = tmp_path / "stride_table.parquet"
    pd.DataFrame({"serotype": ["DENVA"]}).to_parquet(p, index=False)
    with pytest.raises(InputError, match="missing required column"):
        load_stride_table(p)


def test_load_unreadable_raises(tmp_path: Path) -> None:
    p = tmp_path / "stride_table.parquet"
    p.write_bytes(b"not a parquet file")
    with pytest.raises(InputError, match="could not read stride_table"):
        load_stride_table(p)


def test_load_roundtrip(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    df = load_stride_table(paths["stride_table"])
    assert not df.empty
    for col in STRIDE_TABLE_REQUIRED:
        assert col in df.columns


def test_file_digest_stable_and_missing(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    d1 = file_digest(paths["stride_table"])
    d2 = file_digest(paths["stride_table"])
    assert d1 == d2 and len(d1) == 64
    assert file_digest(tmp_path / "absent.parquet") == ""


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def test_build_s4_end_to_end(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    tables, report = build_s4(paths["stride_table"])
    assert report.all_passed
    assert report.serotypes == ["DENVA", "DENVB"]
    assert len(tables.variance_budget) == 6      # 2 serotypes × 3 domains
    assert len(tables.residue_variance) == 8     # 2 × 4 loci
    assert len(tables.significance_screen) == 8  # 2 × 4 gated mechanisms
    assert len(tables.domain_effect_summary) == 6
    # provenance stamped uncalibrated with alpha + digest
    assert report.provenance["calibrated"] is False
    assert report.provenance["fdr_alpha"] == 0.05
    assert report.provenance["inputs"]["stride_table"]["sha256"]
    # facts capture signed / significant counts
    assert report.facts["n_signed_mechanisms"] == 6  # 3 signed × 2 serotypes
    assert report.facts["n_ci_excludes_zero"] == 4   # NS3:51 + NS3:200 × 2


def test_run_s4_writes_all_artifacts(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    out = tmp_path / "out"
    run_s4(paths["stride_table"], out)
    for name in (
        OUT_VARIANCE_BUDGET,
        OUT_RESIDUE_VARIANCE,
        OUT_SIGNIFICANCE_SCREEN,
        OUT_DOMAIN_EFFECT_SUMMARY,
        OUT_UNCERTAINTY_SUMMARY,
    ):
        assert (out / name).exists(), f"missing {name}"
    payload = json.loads((out / OUT_UNCERTAINTY_SUMMARY).read_text())
    assert payload["stage"] == "S4"
    assert payload["all_checks_passed"] is True
    assert payload["provenance"]["calibrated"] is False


def test_run_s4_deterministic(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    t1, _ = run_s4(paths["stride_table"], tmp_path / "o1")
    t2, _ = run_s4(paths["stride_table"], tmp_path / "o2")
    pd.testing.assert_frame_equal(t1.variance_budget, t2.variance_budget)
    pd.testing.assert_frame_equal(t1.residue_variance, t2.residue_variance)
    pd.testing.assert_frame_equal(
        t1.significance_screen, t2.significance_screen
    )
    pd.testing.assert_frame_equal(
        t1.domain_effect_summary, t2.domain_effect_summary
    )


def test_build_s4_empty_dataset(tmp_path: Path) -> None:
    empty = pd.DataFrame(columns=list(STRIDE_TABLE_REQUIRED))
    paths = write_inputs(tmp_path, stride_table=empty)
    tables, report = build_s4(paths["stride_table"])
    assert report.all_passed
    assert tables.variance_budget.empty
    assert tables.significance_screen.empty
    assert report.facts == {}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def test_cli_success(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    write_inputs(tmp_path / "in")
    rc = main(
        [
            "--input-dir",
            str(tmp_path / "in"),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 0
    assert "S4 OK" in capsys.readouterr().out
    assert (tmp_path / "out" / OUT_VARIANCE_BUDGET).exists()


def test_cli_explicit_stride_table(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    paths = write_inputs(tmp_path / "in")
    rc = main(
        [
            "--stride-table",
            str(paths["stride_table"]),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 0
    assert (tmp_path / "out" / OUT_SIGNIFICANCE_SCREEN).exists()


def test_cli_failure_returns_1(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc = main(
        [
            "--input-dir",
            str(tmp_path / "missing"),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 1
    assert "S4 FAILED" in capsys.readouterr().err
