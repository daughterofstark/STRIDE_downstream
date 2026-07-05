"""Tests for the structural validation layer."""
from __future__ import annotations

from pathlib import Path

import pytest

from stride_s7.models import ArtifactRecord, S7Report
from stride_s7.models.errors import ConsistencyError
from stride_s7.s7 import build_s7
from stride_s7.validation import (
    validate_columns,
    validate_completeness,
    validate_filenames,
    validate_on_disk,
    validate_provenance,
)


def test_completeness_passes_on_full_report(stage_dirs: dict[str, Path]) -> None:
    _artifacts, report = build_s7(stage_dirs)
    fresh = S7Report(figures=report.figures, tables=report.tables)
    validate_completeness(fresh)
    assert fresh.all_passed


def test_completeness_flags_missing_figure(stage_dirs: dict[str, Path]) -> None:
    _artifacts, report = build_s7(stage_dirs)
    broken = S7Report(figures=report.figures[:-1], tables=report.tables)
    with pytest.raises(ConsistencyError, match="missing figure"):
        validate_completeness(broken)


def test_completeness_flags_duplicate_table(stage_dirs: dict[str, Path]) -> None:
    _artifacts, report = build_s7(stage_dirs)
    broken = S7Report(
        figures=report.figures, tables=[*report.tables, report.tables[0]]
    )
    with pytest.raises(ConsistencyError, match="duplicate table"):
        validate_completeness(broken)


def test_columns_flag_wrong_schema(
    stage_dirs: dict[str, Path],
) -> None:
    artifacts, report = build_s7(stage_dirs)
    fig = dict(artifacts.figure_data)
    fig["F1"] = fig["F1"].drop(columns=["rho_residue"])
    with pytest.raises(ConsistencyError, match="figure F1"):
        validate_columns(fig, artifacts.tables, S7Report())


def test_columns_flag_wrong_table_schema(stage_dirs: dict[str, Path]) -> None:
    artifacts, _report = build_s7(stage_dirs)
    tabs = dict(artifacts.tables)
    tabs["T5"] = tabs["T5"].drop(columns=["variance_regime"])
    with pytest.raises(ConsistencyError, match="table T5"):
        validate_columns(artifacts.figure_data, tabs, S7Report())


def test_filenames_flag_bad_slug() -> None:
    report = S7Report(
        figures=[
            ArtifactRecord(
                artifact_id="F1",
                kind="figure",
                title="x",
                slug="WRONG_SLUG",
                sources=[],
                files=["WRONG_SLUG.svg", "WRONG_SLUG.csv", "WRONG_SLUG.parquet"],
                n_rows=0,
            )
        ]
    )
    with pytest.raises(ConsistencyError, match="figure F1"):
        validate_filenames(report)


def test_on_disk_flags_missing_file(
    stage_dirs: dict[str, Path], tmp_path: Path
) -> None:
    _artifacts, report = build_s7(stage_dirs)
    empty_dir = tmp_path / "nothing"
    empty_dir.mkdir()
    with pytest.raises(ConsistencyError, match="not on disk"):
        validate_on_disk(report, empty_dir)


def test_provenance_requires_calibrated_false() -> None:
    report = S7Report()
    report.provenance = {
        "calibrated": True,
        "provisional_rho_star": 0.5,
        "n_replicates_note": "x",
        "inputs": {"a": {"sha256": "abc"}},
    }
    with pytest.raises(ConsistencyError, match="calibrated"):
        validate_provenance(report)


def test_provenance_requires_input_digests() -> None:
    report = S7Report()
    report.provenance = {
        "calibrated": False,
        "provisional_rho_star": 0.5,
        "n_replicates_note": "x",
        "inputs": {"a": {"path": "p"}},  # no sha256
    }
    with pytest.raises(ConsistencyError, match="sha256"):
        validate_provenance(report)


def test_provenance_missing_key() -> None:
    report = S7Report()
    report.provenance = {"calibrated": False}
    with pytest.raises(ConsistencyError, match="provenance missing"):
        validate_provenance(report)
