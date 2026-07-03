"""Loader, orchestration, and CLI tests for S3."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from stride_s3 import build_s3, run_s3
from stride_s3.__main__ import main
from stride_s3.io import file_digest, load_stride_table
from stride_s3.models.errors import InputError
from stride_s3.models.schema import (
    OUT_CHAIN_CONTRAST,
    OUT_HIERARCHY_SUMMARY,
    OUT_MONOTONICITY_AUDIT,
    OUT_RESOLUTION_GAP,
    OUT_SCALE_CURVE,
    STRIDE_TABLE_REQUIRED,
)
from tests.s3.fixtures import write_inputs


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
def test_build_s3_end_to_end(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    tables, report = build_s3(paths["stride_table"])
    assert report.all_passed
    assert report.serotypes == ["DENVA", "DENVB"]
    # 2 serotypes × 4 loci × 7 scales
    assert len(tables.scale_curve) == 56
    assert len(tables.resolution_gap) == 8
    assert len(tables.monotonicity_audit) == 8
    assert len(tables.chain_contrast) == 4
    # provenance stamped uncalibrated with a digest
    assert report.provenance["calibrated"] is False
    assert report.provenance["inputs"]["stride_table"]["sha256"]
    # facts capture the distributed + non-monotone counts
    assert report.facts["n_non_monotone_loci"] == 2  # NS3:99 in each serotype
    assert report.facts["n_distributed_effects"] >= 1


def test_run_s3_writes_all_artifacts(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    out = tmp_path / "out"
    run_s3(paths["stride_table"], out)
    for name in (
        OUT_SCALE_CURVE,
        OUT_RESOLUTION_GAP,
        OUT_MONOTONICITY_AUDIT,
        OUT_CHAIN_CONTRAST,
        OUT_HIERARCHY_SUMMARY,
    ):
        assert (out / name).exists(), f"missing {name}"
    payload = json.loads((out / OUT_HIERARCHY_SUMMARY).read_text())
    assert payload["stage"] == "S3"
    assert payload["all_checks_passed"] is True
    assert payload["provenance"]["calibrated"] is False


def test_run_s3_deterministic(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    t1, _ = run_s3(paths["stride_table"], tmp_path / "o1")
    t2, _ = run_s3(paths["stride_table"], tmp_path / "o2")
    pd.testing.assert_frame_equal(t1.scale_curve, t2.scale_curve)
    pd.testing.assert_frame_equal(t1.resolution_gap, t2.resolution_gap)
    pd.testing.assert_frame_equal(t1.monotonicity_audit, t2.monotonicity_audit)
    pd.testing.assert_frame_equal(t1.chain_contrast, t2.chain_contrast)


def test_build_s3_empty_dataset(tmp_path: Path) -> None:
    empty = pd.DataFrame(columns=list(STRIDE_TABLE_REQUIRED))
    paths = write_inputs(tmp_path, stride_table=empty)
    tables, report = build_s3(paths["stride_table"])
    assert report.all_passed
    assert tables.scale_curve.empty
    assert tables.chain_contrast.empty
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
    assert "S3 OK" in capsys.readouterr().out
    assert (tmp_path / "out" / OUT_SCALE_CURVE).exists()


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
    assert (tmp_path / "out" / OUT_CHAIN_CONTRAST).exists()


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
    assert "S3 FAILED" in capsys.readouterr().err
