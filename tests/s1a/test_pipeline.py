"""Loader, orchestration, and CLI tests for S1A."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from stride_s1a import build_s1a, run_s1a
from stride_s1a.__main__ import main
from stride_s1a.io import load_replicate_table, load_stride_table
from stride_s1a.models.errors import InputError
from stride_s1a.models.schema import (
    OUT_CANONICAL_RESIDUES,
    OUT_CONSERVATION_TABLE,
    OUT_DATASET_SUMMARY,
    OUT_DOMAIN_TABLE,
    OUT_REPLICATE_INVENTORY,
)
from tests.s1a.fixtures import write_tables


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------
def test_load_stride_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="stride_table not found"):
        load_stride_table(tmp_path / "nope.parquet")


def test_load_stride_missing_columns_raises(tmp_path: Path) -> None:
    p = tmp_path / "stride_table.parquet"
    pd.DataFrame({"serotype": ["DENVA"]}).to_parquet(p, index=False)
    with pytest.raises(InputError, match="missing required column"):
        load_stride_table(p)


def test_load_replicate_missing_file_is_empty(tmp_path: Path) -> None:
    # a missing replicate table is a valid summaries-only state
    df = load_replicate_table(tmp_path / "absent.parquet")
    assert df.empty


def test_load_replicate_missing_columns_raises(tmp_path: Path) -> None:
    p = tmp_path / "replicate_table.parquet"
    pd.DataFrame({"serotype": ["DENVA"]}).to_parquet(p, index=False)
    with pytest.raises(InputError, match="missing required column"):
        load_replicate_table(p)


def test_load_stride_unreadable_raises(tmp_path: Path) -> None:
    p = tmp_path / "stride_table.parquet"
    p.write_bytes(b"not a parquet file")
    with pytest.raises(InputError, match="could not read stride_table"):
        load_stride_table(p)


def test_load_replicate_unreadable_raises(tmp_path: Path) -> None:
    p = tmp_path / "replicate_table.parquet"
    p.write_bytes(b"not a parquet file")
    with pytest.raises(InputError, match="could not read replicate_table"):
        load_replicate_table(p)


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def test_build_s1a_end_to_end(tmp_path: Path) -> None:
    sp, rp = write_tables(tmp_path)
    tables, report = build_s1a(sp, rp)
    assert report.all_passed
    assert report.n_canonical_residues == 6
    assert report.n_union == 3
    assert report.n_conserved_all == 1
    assert len(tables.canonical_residues) == 6
    assert len(tables.conservation_table) == 3


def test_run_s1a_writes_all_artifacts(tmp_path: Path) -> None:
    sp, rp = write_tables(tmp_path / "in")
    out = tmp_path / "out"
    tables, report = run_s1a(sp, rp, out)
    for name in (
        OUT_CANONICAL_RESIDUES,
        OUT_DOMAIN_TABLE,
        OUT_REPLICATE_INVENTORY,
        OUT_CONSERVATION_TABLE,
        OUT_DATASET_SUMMARY,
    ):
        assert (out / name).exists(), f"missing {name}"


def test_run_s1a_deterministic(tmp_path: Path) -> None:
    sp, rp = write_tables(tmp_path / "in")
    t1, _ = run_s1a(sp, rp, tmp_path / "o1")
    t2, _ = run_s1a(sp, rp, tmp_path / "o2")
    pd.testing.assert_frame_equal(t1.canonical_residues, t2.canonical_residues)
    pd.testing.assert_frame_equal(t1.conservation_table, t2.conservation_table)
    pd.testing.assert_frame_equal(t1.domain_table, t2.domain_table)
    pd.testing.assert_frame_equal(t1.replicate_inventory, t2.replicate_inventory)


def test_build_s1a_summaries_only(tmp_path: Path) -> None:
    # no replicate table at all -> inventory all-unavailable, still valid
    from tests.s1a.fixtures import make_stride_table

    sp = tmp_path / "stride_table.parquet"
    make_stride_table().to_parquet(sp, index=False)
    tables, report = build_s1a(sp, tmp_path / "absent_replicate.parquet")
    assert report.all_passed
    assert (tables.replicate_inventory["n_replicates"] == 0).all()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def test_cli_success(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    write_tables(tmp_path / "in")
    rc = main([
        "--input-dir", str(tmp_path / "in"),
        "--output-dir", str(tmp_path / "out"),
    ])
    assert rc == 0
    assert "S1A OK" in capsys.readouterr().out
    assert (tmp_path / "out" / OUT_CANONICAL_RESIDUES).exists()


def test_cli_failure_returns_1(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    rc = main([
        "--stride-table", str(tmp_path / "missing.parquet"),
        "--output-dir", str(tmp_path / "out"),
    ])
    assert rc == 1
    assert "S1A FAILED" in capsys.readouterr().err
