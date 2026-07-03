"""Loader, orchestration, and CLI tests for S2."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from stride_s2 import build_s2, run_s2
from stride_s2.__main__ import main
from stride_s2.io import (
    file_digest,
    load_domain_annotation,
    load_residue_annotation,
    load_stride_table,
)
from stride_s2.models.errors import ConfigError, InputError
from stride_s2.models.schema import (
    DEFAULT_RHO_STAR_BAND,
    OUT_DOMAIN_REPRODUCIBILITY,
    OUT_REDUCTION_SUMMARY,
    OUT_RESIDUE_LANDSCAPE,
    OUT_RESOLUTION_CENSUS,
    OUT_SEROTYPE_SUMMARY,
    OUT_SIGNED_SCREEN,
)
from tests.s2.fixtures import write_inputs


# ---------------------------------------------------------------------------
# loaders
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
    p = tmp_path / "residue_annotation.parquet"
    p.write_bytes(b"not a parquet file")
    with pytest.raises(InputError, match="could not read residue_annotation"):
        load_residue_annotation(p)


def test_all_loaders_roundtrip(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    assert not load_stride_table(paths["stride_table"]).empty
    assert not load_residue_annotation(paths["residue_annotation"]).empty
    assert not load_domain_annotation(paths["domain_annotation"]).empty


def test_file_digest_stable_and_missing(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    d1 = file_digest(paths["stride_table"])
    d2 = file_digest(paths["stride_table"])
    assert d1 == d2 and len(d1) == 64
    assert file_digest(tmp_path / "absent.parquet") == ""


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------
def test_build_s2_end_to_end(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    tables, report = build_s2(
        paths["stride_table"],
        paths["residue_annotation"],
        paths["domain_annotation"],
    )
    assert report.all_passed
    assert report.serotypes == ["DENVA", "DENVB"]
    assert report.rho_star_band == list(DEFAULT_RHO_STAR_BAND)
    # 2 serotypes × 3 loci residue-landscape rows
    assert len(tables.residue_landscape) == 6
    # 2 serotypes × 3 domains
    assert len(tables.domain_reproducibility) == 6
    # provenance header present and stamped uncalibrated
    assert report.provenance["calibrated"] is False
    assert report.provenance["inputs"]["stride_table"]["sha256"]


def test_build_s2_custom_band(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    tables, report = build_s2(
        paths["stride_table"],
        paths["residue_annotation"],
        paths["domain_annotation"],
        rho_star_band=[0.9, 0.5, 0.5],  # unsorted + duplicate
    )
    assert report.rho_star_band == [0.5, 0.9]
    # census has both ρ* values per serotype
    assert set(tables.resolution_census["rho_star"].unique()) == {0.5, 0.9}


def test_build_s2_empty_band_raises(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    with pytest.raises(ConfigError, match="empty"):
        build_s2(
            paths["stride_table"],
            paths["residue_annotation"],
            paths["domain_annotation"],
            rho_star_band=[],
        )


def test_build_s2_out_of_range_band_raises(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    with pytest.raises(ConfigError, match=r"outside \[0, 1\]"):
        build_s2(
            paths["stride_table"],
            paths["residue_annotation"],
            paths["domain_annotation"],
            rho_star_band=[0.5, 1.5],
        )


def test_run_s2_writes_all_artifacts(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    out = tmp_path / "out"
    run_s2(
        paths["stride_table"],
        paths["residue_annotation"],
        paths["domain_annotation"],
        out,
    )
    for name in (
        OUT_RESOLUTION_CENSUS,
        OUT_RESIDUE_LANDSCAPE,
        OUT_DOMAIN_REPRODUCIBILITY,
        OUT_SIGNED_SCREEN,
        OUT_SEROTYPE_SUMMARY,
        OUT_REDUCTION_SUMMARY,
    ):
        assert (out / name).exists(), f"missing {name}"
    # the summary JSON carries the provenance header and passed checks
    payload = json.loads((out / OUT_REDUCTION_SUMMARY).read_text())
    assert payload["stage"] == "S2"
    assert payload["all_checks_passed"] is True
    assert payload["provenance"]["calibrated"] is False


def test_run_s2_deterministic(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path / "in")
    args = (
        paths["stride_table"],
        paths["residue_annotation"],
        paths["domain_annotation"],
    )
    t1, _ = run_s2(*args, tmp_path / "o1")
    t2, _ = run_s2(*args, tmp_path / "o2")
    pd.testing.assert_frame_equal(
        t1.resolution_census, t2.resolution_census
    )
    pd.testing.assert_frame_equal(t1.residue_landscape, t2.residue_landscape)
    pd.testing.assert_frame_equal(
        t1.domain_reproducibility, t2.domain_reproducibility
    )
    pd.testing.assert_frame_equal(t1.signed_screen, t2.signed_screen)
    pd.testing.assert_frame_equal(t1.serotype_summary, t2.serotype_summary)


def test_build_s2_empty_dataset(tmp_path: Path) -> None:
    from stride_s2.models.schema import (
        DOMAIN_ANNOTATION_REQUIRED,
        RESIDUE_ANNOTATION_REQUIRED,
        STRIDE_TABLE_REQUIRED,
    )

    empty_st = pd.DataFrame(columns=list(STRIDE_TABLE_REQUIRED))
    empty_ra = pd.DataFrame(columns=list(RESIDUE_ANNOTATION_REQUIRED))
    empty_da = pd.DataFrame(columns=list(DOMAIN_ANNOTATION_REQUIRED))
    paths = write_inputs(
        tmp_path,
        stride_table=empty_st,
        residue_annotation=empty_ra,
        domain_annotation=empty_da,
    )
    tables, report = build_s2(
        paths["stride_table"],
        paths["residue_annotation"],
        paths["domain_annotation"],
    )
    assert report.all_passed
    assert tables.resolution_census.empty
    assert tables.residue_landscape.empty
    assert report.facts == {}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def test_cli_success(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    write_inputs(tmp_path / "in")
    rc = main([
        "--stride-input-dir", str(tmp_path / "in"),
        "--annotation-input-dir", str(tmp_path / "in"),
        "--output-dir", str(tmp_path / "out"),
    ])
    assert rc == 0
    assert "S2 OK" in capsys.readouterr().out
    assert (tmp_path / "out" / OUT_RESOLUTION_CENSUS).exists()


def test_cli_custom_rho_star(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    write_inputs(tmp_path / "in")
    rc = main([
        "--stride-input-dir", str(tmp_path / "in"),
        "--annotation-input-dir", str(tmp_path / "in"),
        "--output-dir", str(tmp_path / "out"),
        "--rho-star", "0.5", "0.8",
    ])
    assert rc == 0
    census = pd.read_parquet(tmp_path / "out" / OUT_RESOLUTION_CENSUS)
    assert set(census["rho_star"].unique()) == {0.5, 0.8}


def test_cli_failure_returns_1(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc = main([
        "--stride-input-dir", str(tmp_path / "missing"),
        "--annotation-input-dir", str(tmp_path / "missing"),
        "--output-dir", str(tmp_path / "out"),
    ])
    assert rc == 1
    assert "S2 FAILED" in capsys.readouterr().err
