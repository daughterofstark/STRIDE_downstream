"""End-to-end tests for the S6 orchestration and CLI."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from stride_s6 import __version__, build_s6, run_s6
from stride_s6.__main__ import main
from stride_s6.models.errors import InputError
from tests.s6.fixtures import (
    make_empty_inventory,
    make_replicate_inventory,
    make_replicate_table,
    write_inputs,
)


def test_version() -> None:
    assert __version__ == "0.1.0"


def test_build_s6_with_effects(tmp_path: Path) -> None:
    inv, tbl = write_inputs(
        tmp_path, make_replicate_inventory(), make_replicate_table()
    )
    tables, report = build_s6(inv, tbl)
    assert report.all_passed
    assert report.per_replicate_effects_available is True
    assert len(tables.replicate_regime) == 4
    assert len(tables.replicate_concordance) == 4
    assert not tables.replicate_effect_spread.empty
    # DENV1's perfect concordance survives the full pipeline
    denv1 = tables.replicate_concordance.query("serotype == 'DENV1'").iloc[0]
    assert denv1["kendalls_w"] == 1.0


def test_build_s6_blocked_state(tmp_path: Path) -> None:
    # inventory present, no replicate table → the design's blocked state
    inv, tbl = write_inputs(tmp_path, make_empty_inventory(), None)
    assert tbl is None
    tables, report = build_s6(inv, None)
    assert report.all_passed
    assert report.per_replicate_effects_available is False
    assert tables.replicate_concordance.empty
    assert tables.replicate_effect_spread.empty
    # the ledger still records all three analyses, all blocked
    assert len(tables.replicate_blocked_analyses) == 3
    assert not tables.replicate_blocked_analyses["available"].any()


def test_run_s6_writes_all_artifacts(tmp_path: Path) -> None:
    inv, tbl = write_inputs(
        tmp_path / "in", make_replicate_inventory(), make_replicate_table()
    )
    out = tmp_path / "out"
    run_s6(inv, out, replicate_table_path=tbl)
    for name in (
        "replicate_regime.parquet",
        "replicate_effect_spread.parquet",
        "replicate_concordance.parquet",
        "replicate_blocked_analyses.parquet",
        "replicate_summary.json",
    ):
        assert (out / name).is_file()
    summary = json.loads((out / "replicate_summary.json").read_text())
    assert summary["stage"] == "S6"
    assert summary["all_checks_passed"] is True
    assert summary["per_replicate_effects_available"] is True
    assert len(summary["blocked_analyses"]) == 3


def test_run_s6_deterministic(tmp_path: Path) -> None:
    inv, tbl = write_inputs(
        tmp_path / "in", make_replicate_inventory(), make_replicate_table()
    )
    out1, out2 = tmp_path / "o1", tmp_path / "o2"
    run_s6(inv, out1, replicate_table_path=tbl)
    run_s6(inv, out2, replicate_table_path=tbl)
    for name in (
        "replicate_regime.parquet",
        "replicate_effect_spread.parquet",
        "replicate_concordance.parquet",
        "replicate_blocked_analyses.parquet",
    ):
        a = pd.read_parquet(out1 / name)
        b = pd.read_parquet(out2 / name)
        pd.testing.assert_frame_equal(a, b)


def test_missing_inventory_raises(tmp_path: Path) -> None:
    with pytest.raises(InputError):
        build_s6(tmp_path / "does_not_exist.parquet", None)


def test_cli_success_with_effects(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    write_inputs(
        tmp_path, make_replicate_inventory(), make_replicate_table()
    )
    rc = main(
        [
            "--input-dir",
            str(tmp_path),
            "--inventory-input-dir",
            str(tmp_path),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 0
    captured = capsys.readouterr()
    assert "S6 OK" in captured.out
    assert "per-run effects available: True" in captured.out


def test_cli_blocked_path_with_flag(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    write_inputs(tmp_path, make_replicate_inventory(), make_replicate_table())
    rc = main(
        [
            "--input-dir",
            str(tmp_path),
            "--inventory-input-dir",
            str(tmp_path),
            "--output-dir",
            str(tmp_path / "out"),
            "--no-replicate-table",
        ]
    )
    assert rc == 0
    captured = capsys.readouterr()
    assert "per-run effects available: False" in captured.out


def test_cli_failure_missing_inventory(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc = main(
        [
            "--inventory-input-dir",
            str(tmp_path / "nope"),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert rc == 1
    captured = capsys.readouterr()
    assert "S6 FAILED" in captured.err
