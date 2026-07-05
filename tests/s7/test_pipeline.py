"""End-to-end tests for S7 orchestration and the CLI."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from stride_s7.__main__ import main
from stride_s7.models.errors import InputError
from stride_s7.models.schema import (
    FIGURE_IDS,
    FIGURE_SLUGS,
    OUT_MANIFEST,
    OUT_SUMMARY,
    TABLE_IDS,
    TABLE_SLUGS,
)
from stride_s7.s7 import build_s7, run_s7


def test_build_s7_passes_all_checks(stage_dirs: dict[str, Path]) -> None:
    _artifacts, report = build_s7(stage_dirs)
    assert report.all_passed
    assert len(report.figures) == len(FIGURE_IDS)
    assert len(report.tables) == len(TABLE_IDS)
    assert report.serotypes == ["DENV1", "DENV2", "DENV3", "DENV4"]


def test_run_s7_writes_every_artifact(
    stage_dirs: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "out"
    _artifacts, report = run_s7(stage_dirs, out)
    assert report.all_passed
    for fid in FIGURE_IDS:
        slug = FIGURE_SLUGS[fid]
        for suffix in (".svg", ".csv", ".parquet"):
            assert (out / f"{slug}{suffix}").is_file()
    for tid in TABLE_IDS:
        slug = TABLE_SLUGS[tid]
        for suffix in (".csv", ".parquet", ".md"):
            assert (out / f"{slug}{suffix}").is_file()
    assert (out / OUT_SUMMARY).is_file()
    assert (out / OUT_MANIFEST).is_file()


def test_summary_json_provenance_and_limitations(
    stage_dirs: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "out"
    run_s7(stage_dirs, out)
    payload = json.loads((out / OUT_SUMMARY).read_text())
    assert payload["stage"] == "S7"
    assert payload["all_checks_passed"] is True
    prov = payload["provenance"]
    assert prov["calibrated"] is False
    assert prov["provisional_rho_star"] == 0.5
    assert len(prov["inputs"]) == 12
    for meta in prov["inputs"].values():
        assert len(meta["sha256"]) == 64
    # S6 replicate layer surfaced as limitations, not fabricated into a figure
    assert any(item["source"] == "S6" for item in payload["limitations"])


def test_run_s7_is_byte_identical(
    stage_dirs: dict[str, Path], tmp_path: Path
) -> None:
    out_a = tmp_path / "a"
    out_b = tmp_path / "b"
    run_s7(stage_dirs, out_a)
    run_s7(stage_dirs, out_b)
    names = sorted(p.name for p in out_a.iterdir())
    assert names == sorted(p.name for p in out_b.iterdir())
    for name in names:
        assert (out_a / name).read_bytes() == (out_b / name).read_bytes(), name


def test_empty_inputs_are_handled(
    empty_stage_dirs: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "empty_out"
    _artifacts, report = run_s7(empty_stage_dirs, out)
    assert report.all_passed
    # all artifacts still written, figures still valid SVG
    for fid in FIGURE_IDS:
        svg = (out / f"{FIGURE_SLUGS[fid]}.svg").read_text()
        assert svg.startswith("<svg")
    assert all(rec.n_rows == 0 for rec in report.figures)


def test_missing_input_raises(stage_dirs: dict[str, Path]) -> None:
    # delete one required S4 table → hard failure
    (stage_dirs["s4"] / "variance_budget.parquet").unlink()
    with pytest.raises(InputError, match="variance_budget"):
        build_s7(stage_dirs)


def test_malformed_input_raises(stage_dirs: dict[str, Path]) -> None:
    # a present file missing a required column
    import pandas as pd

    path = stage_dirs["s5"] / "position_conservation.parquet"
    df = pd.read_parquet(path).drop(columns=["conservation_class"])
    df.to_parquet(path, index=False)
    with pytest.raises(InputError, match="conservation_class"):
        build_s7(stage_dirs)


def test_cli_success(
    stage_dirs: dict[str, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "cli_out"
    rc = main(
        [
            "--s2-input-dir", str(stage_dirs["s2"]),
            "--s3-input-dir", str(stage_dirs["s3"]),
            "--s4-input-dir", str(stage_dirs["s4"]),
            "--s5-input-dir", str(stage_dirs["s5"]),
            "--s6-input-dir", str(stage_dirs["s6"]),
            "--output-dir", str(out),
        ]
    )
    assert rc == 0
    captured = capsys.readouterr()
    assert "S7 OK" in captured.out
    assert (out / OUT_SUMMARY).is_file()


def test_cli_failure_on_missing_dir(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc = main(
        [
            "--s2-input-dir", str(tmp_path / "nope"),
            "--s3-input-dir", str(tmp_path / "nope"),
            "--s4-input-dir", str(tmp_path / "nope"),
            "--s5-input-dir", str(tmp_path / "nope"),
            "--s6-input-dir", str(tmp_path / "nope"),
            "--output-dir", str(tmp_path / "out"),
        ]
    )
    assert rc == 1
    captured = capsys.readouterr()
    assert "S7 FAILED" in captured.err
