"""Integration tests: end-to-end S0 orchestration on synthetic datasets."""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from stride_analysis import build_tables, run_s0
from stride_analysis._synthetic import write_dataset
from stride_analysis.models.schema import (
    OUT_REPLICATE_PARQUET,
    OUT_SCHEMA_REPORT,
    OUT_STRIDE_PARQUET,
    OUT_VALIDATION_REPORT,
    STRIDE_TABLE_COLUMNS,
)


def test_build_tables_end_to_end(dataset_root: Path) -> None:
    rep, stride, report = build_tables(dataset_root)
    assert report.all_passed
    # 2 serotypes x 3 replicates x 2 residues
    assert len(rep) == 12
    # 2 serotypes x 2 loci x 7 scales
    assert len(stride) == 28
    assert not rep.duplicated(["serotype", "replicate", "canon_label"]).any()
    assert not stride.duplicated(["serotype", "canon_label", "scale_level"]).any()


def test_run_s0_writes_all_artifacts(dataset_root: Path, tmp_path: Path) -> None:
    out = tmp_path / "outputs"
    rep, stride, report = run_s0(dataset_root, out)
    for name in (
        OUT_REPLICATE_PARQUET,
        "replicate_table.csv",
        OUT_STRIDE_PARQUET,
        "stride_table.csv",
        OUT_SCHEMA_REPORT,
        OUT_VALIDATION_REPORT,
    ):
        assert (out / name).exists(), f"missing {name}"
    rt = pd.read_parquet(out / OUT_STRIDE_PARQUET)
    assert list(rt.columns) == list(STRIDE_TABLE_COLUMNS)


def test_run_s0_deterministic(dataset_root: Path, tmp_path: Path) -> None:
    r1, s1, _ = run_s0(dataset_root, tmp_path / "o1")
    r2, s2, _ = run_s0(dataset_root, tmp_path / "o2")
    pd.testing.assert_frame_equal(r1, r2)
    pd.testing.assert_frame_equal(s1, s2)


def test_summaries_only_dataset(tmp_path: Path) -> None:
    root = write_dataset(tmp_path / "d", ["DENV1"], ["1st_run"], with_summaries=True)
    shutil.rmtree(root / "1st_run")
    rep, stride, report = build_tables(root, require_replicates=False)
    assert rep.empty
    assert len(stride) == 14
    assert report.all_passed


def test_replicates_only_dataset(tmp_path: Path) -> None:
    root = write_dataset(
        tmp_path / "d", ["DENV1"], ["1st_run", "2nd_run"], with_summaries=False
    )
    rep, stride, report = build_tables(root, require_summaries=False)
    assert len(rep) == 4  # 2 replicates x 2 residues
    assert stride.empty
    assert report.all_passed


def test_cross_level_alignment_reported(dataset_root: Path) -> None:
    _, _, report = build_tables(dataset_root)
    align = [c for c in report.checks if "alignment" in c.name]
    assert align and all(c.passed for c in align)
    assert "shared" in align[0].detail


def test_example_dataset_ships_and_runs(tmp_path: Path) -> None:
    # the committed example must load out of the box
    example = Path("examples/small_synthetic_dataset")
    if not example.is_dir():
        pytest.skip("example dataset not present in this checkout")
    rep, stride, report = build_tables(example)
    assert report.all_passed
    assert len(rep) == 24  # 4 serotypes x 3 replicates x 2 residues
    assert len(stride) == 56
