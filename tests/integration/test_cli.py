"""CLI tests."""
from __future__ import annotations

import shutil
from pathlib import Path

from stride_analysis.__main__ import main
from stride_analysis._synthetic import write_dataset


def test_cli_success(dataset_root: Path, tmp_path: Path, capsys) -> None:
    rc = main(["--data-root", str(dataset_root), "--output-dir", str(tmp_path / "o")])
    assert rc == 0
    out = capsys.readouterr().out
    assert "S0 OK" in out
    assert (tmp_path / "o" / "stride_table.parquet").exists()


def test_cli_failure_returns_1(tmp_path: Path, capsys) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    rc = main(["--data-root", str(empty), "--output-dir", str(tmp_path / "o")])
    assert rc == 1
    assert "S0 FAILED" in capsys.readouterr().err


def test_cli_no_require_replicates(tmp_path: Path, capsys) -> None:
    root = write_dataset(tmp_path / "d", ["DENV1"], ["1st_run"], with_summaries=True)
    shutil.rmtree(root / "1st_run")
    rc = main([
        "--data-root", str(root),
        "--output-dir", str(tmp_path / "o"),
        "--no-require-replicates",
    ])
    assert rc == 0
