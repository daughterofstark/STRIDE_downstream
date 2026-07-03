"""Loader, orchestration, and CLI tests for S1B."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from stride_s1b import build_s1b, run_s1b
from stride_s1b.__main__ import main
from stride_s1b.io import (
    load_canonical_residues,
    load_conservation_table,
    load_domain_table,
    load_replicate_inventory,
)
from stride_s1b.models.errors import InputError
from stride_s1b.models.schema import (
    OUT_ANNOTATION_SUMMARY,
    OUT_DOMAIN_ANNOTATION,
    OUT_HIERARCHY_ANNOTATION,
    OUT_RESIDUE_ANNOTATION,
    OUT_SEROTYPE_ANNOTATION,
)
from tests.s1b.fixtures import make_empty_replicate_inventory, write_s1a_tables


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------
def test_load_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="canonical_residues not found"):
        load_canonical_residues(tmp_path / "nope.parquet")


def test_load_missing_columns_raises(tmp_path: Path) -> None:
    p = tmp_path / "canonical_residues.parquet"
    pd.DataFrame({"serotype": ["DENVA"]}).to_parquet(p, index=False)
    with pytest.raises(InputError, match="missing required column"):
        load_canonical_residues(p)


def test_load_unreadable_raises(tmp_path: Path) -> None:
    p = tmp_path / "domain_table.parquet"
    p.write_bytes(b"not a parquet file")
    with pytest.raises(InputError, match="could not read domain_table"):
        load_domain_table(p)


def test_all_loaders_roundtrip(tmp_path: Path) -> None:
    paths = write_s1a_tables(tmp_path)
    assert not load_canonical_residues(paths["canonical_residues"]).empty
    assert not load_domain_table(paths["domain_table"]).empty
    assert not load_replicate_inventory(paths["replicate_inventory"]).empty
    assert not load_conservation_table(paths["conservation_table"]).empty


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def test_build_s1b_end_to_end(tmp_path: Path) -> None:
    paths = write_s1a_tables(tmp_path)
    tables, report = build_s1b(
        paths["canonical_residues"],
        paths["domain_table"],
        paths["replicate_inventory"],
        paths["conservation_table"],
    )
    assert report.all_passed
    # 3 + 2 + 1 = 6 residues across DENVA/DENVB/DENVC
    assert report.n_residue_annotations == 6
    assert report.n_serotype_annotations == 3
    assert len(tables.residue_annotation) == 6
    assert len(tables.hierarchy_annotation) == 6


def test_run_s1b_writes_all_artifacts(tmp_path: Path) -> None:
    paths = write_s1a_tables(tmp_path / "in")
    out = tmp_path / "out"
    run_s1b(
        paths["canonical_residues"],
        paths["domain_table"],
        paths["replicate_inventory"],
        paths["conservation_table"],
        out,
    )
    for name in (
        OUT_RESIDUE_ANNOTATION,
        OUT_DOMAIN_ANNOTATION,
        OUT_HIERARCHY_ANNOTATION,
        OUT_SEROTYPE_ANNOTATION,
        OUT_ANNOTATION_SUMMARY,
    ):
        assert (out / name).exists(), f"missing {name}"


def test_run_s1b_deterministic(tmp_path: Path) -> None:
    paths = write_s1a_tables(tmp_path / "in")
    args = (
        paths["canonical_residues"],
        paths["domain_table"],
        paths["replicate_inventory"],
        paths["conservation_table"],
    )
    t1, _ = run_s1b(*args, tmp_path / "o1")
    t2, _ = run_s1b(*args, tmp_path / "o2")
    pd.testing.assert_frame_equal(t1.residue_annotation, t2.residue_annotation)
    pd.testing.assert_frame_equal(t1.domain_annotation, t2.domain_annotation)
    pd.testing.assert_frame_equal(t1.hierarchy_annotation, t2.hierarchy_annotation)
    pd.testing.assert_frame_equal(t1.serotype_annotation, t2.serotype_annotation)


def test_build_s1b_summaries_only(tmp_path: Path) -> None:
    # replicate inventory with zero replicates everywhere -> still valid
    paths = write_s1a_tables(
        tmp_path, replicate_inventory=make_empty_replicate_inventory()
    )
    tables, report = build_s1b(
        paths["canonical_residues"],
        paths["domain_table"],
        paths["replicate_inventory"],
        paths["conservation_table"],
    )
    assert report.all_passed
    assert (tables.residue_annotation["n_replicates"] == 0).all()


def test_build_s1b_empty_dataset(tmp_path: Path) -> None:
    # all four S1A tables empty (well-formed headers) -> empty annotations, valid
    from stride_s1b.models.schema import (
        CANONICAL_RESIDUES_REQUIRED,
        CONSERVATION_TABLE_REQUIRED,
        DOMAIN_TABLE_REQUIRED,
        REPLICATE_INVENTORY_REQUIRED,
    )

    empty_cr = pd.DataFrame(columns=list(CANONICAL_RESIDUES_REQUIRED))
    empty_dt = pd.DataFrame(columns=list(DOMAIN_TABLE_REQUIRED))
    empty_ri = pd.DataFrame(columns=list(REPLICATE_INVENTORY_REQUIRED))
    empty_ct = pd.DataFrame(columns=list(CONSERVATION_TABLE_REQUIRED))
    paths = write_s1a_tables(
        tmp_path,
        canonical_residues=empty_cr,
        domain_table=empty_dt,
        replicate_inventory=empty_ri,
        conservation_table=empty_ct,
    )
    tables, report = build_s1b(
        paths["canonical_residues"],
        paths["domain_table"],
        paths["replicate_inventory"],
        paths["conservation_table"],
    )
    assert report.all_passed
    assert report.n_residue_annotations == 0
    assert tables.residue_annotation.empty
    assert report.facts == {}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def test_cli_success(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    write_s1a_tables(tmp_path / "in")
    rc = main([
        "--input-dir", str(tmp_path / "in"),
        "--output-dir", str(tmp_path / "out"),
    ])
    assert rc == 0
    assert "S1B OK" in capsys.readouterr().out
    assert (tmp_path / "out" / OUT_RESIDUE_ANNOTATION).exists()


def test_cli_failure_returns_1(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc = main([
        "--input-dir", str(tmp_path / "missing"),
        "--output-dir", str(tmp_path / "out"),
    ])
    assert rc == 1
    assert "S1B FAILED" in capsys.readouterr().err
