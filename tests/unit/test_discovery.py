"""Discovery tests, including the required new failure modes:
missing replicate, missing profile, missing mechanism, inconsistent replicate
counts, malformed layouts, duplicate serotypes, partial datasets.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from stride_analysis import discover_dataset
from stride_analysis._synthetic import (
    make_mechanism_dict,
    make_profile_df,
    write_dataset,
    write_real_dataset,
)
from stride_analysis.io import discover_replicates
from stride_analysis.models.errors import DiscoveryError


# -- native container layouts (Bug 1 regression) ----------------------------
def test_discovers_native_per_run_container(tmp_path: Path) -> None:
    """The native STRIDE ``per_run/`` container is discovered without renaming."""
    root = write_real_dataset(
        tmp_path / "stride_outs",
        ["DENV1", "DENV2"],
        ["1st_run", "2nd_run", "3rd_run"],
        container="per_run",
    )
    ds = discover_dataset(root)
    assert ds.names == ("DENV1", "DENV2")
    for s in ds.serotypes:
        assert s.n_replicates == 3
        assert s.summary is not None


def test_discovers_replicates_container(tmp_path: Path) -> None:
    """The flat ``replicates/`` container still works (unchanged behaviour)."""
    root = write_real_dataset(
        tmp_path / "d",
        ["DENV1", "DENV2"],
        ["1st_run", "2nd_run", "3rd_run"],
        container="replicates",
    )
    ds = discover_dataset(root)
    assert ds.names == ("DENV1", "DENV2")
    reps = discover_replicates(root)
    ordered = sorted(reps["DENV1"], key=lambda r: r.replicate_index)
    assert [r.run_dir for r in ordered] == ["1st_run", "2nd_run", "3rd_run"]


def test_discovers_run_dirs_directly_under_root(tmp_path: Path) -> None:
    """With no container, run dirs directly under the root still discover."""
    root = write_real_dataset(
        tmp_path / "d",
        ["DENV1", "DENV2"],
        ["1st_run", "2nd_run"],
        container=None,
    )
    ds = discover_dataset(root)
    assert ds.names == ("DENV1", "DENV2")
    for s in ds.serotypes:
        assert s.n_replicates == 2


def test_discovers_full_dataset(dataset_root: Path) -> None:
    ds = discover_dataset(dataset_root)
    assert ds.names == ("DENV1", "DENV2")
    for s in ds.serotypes:
        assert s.n_replicates == 3
        assert s.summary is not None


def test_replicate_index_follows_run_order(dataset_root: Path) -> None:
    reps = discover_replicates(dataset_root)
    d1 = sorted(reps["DENV1"], key=lambda r: r.replicate_index)
    assert [r.run_dir for r in d1] == ["1st_run", "2nd_run", "3rd_run"]
    assert [r.replicate_index for r in d1] == [1, 2, 3]


def test_missing_root_raises(tmp_path: Path) -> None:
    with pytest.raises(DiscoveryError, match="does not exist"):
        discover_dataset(tmp_path / "nope")


def test_empty_root_raises(tmp_path: Path) -> None:
    d = tmp_path / "empty"
    d.mkdir()
    with pytest.raises(DiscoveryError, match="no STRIDE data found"):
        discover_dataset(d)


# -- missing replicate ------------------------------------------------------
def test_missing_replicate_for_serotype_raises(dataset_root: Path) -> None:
    # remove DENV2's 2nd_run replicate entirely -> unequal counts
    shutil.rmtree(dataset_root / "2nd_run" / "DENV2")
    with pytest.raises(DiscoveryError, match="inconsistent replicate counts"):
        discover_dataset(dataset_root)


def test_serotype_with_no_replicates_when_required(tmp_path: Path) -> None:
    # DENV1 has replicates + summary; DENV2 has only a summary
    root = write_dataset(tmp_path / "d", ["DENV1"], ["1st_run", "2nd_run"])
    # add a summary-only DENV2
    make_profile_df("DENV2").to_csv(
        root / "summaries" / "DENV2_profile.csv", index=False
    )
    (root / "summaries" / "DENV2_mechanism.json").write_text(
        json.dumps(make_mechanism_dict("DENV2"))
    )
    with pytest.raises(DiscoveryError, match="missing replicate data"):
        discover_dataset(root, require_replicates=True)


# -- missing profile / mechanism --------------------------------------------
def test_missing_profile_raises(dataset_root: Path) -> None:
    (dataset_root / "summaries" / "DENV2_profile.csv").unlink()
    with pytest.raises(DiscoveryError, match="mechanisms without profile"):
        discover_dataset(dataset_root)


def test_missing_mechanism_raises(dataset_root: Path) -> None:
    (dataset_root / "summaries" / "DENV2_mechanism.json").unlink()
    with pytest.raises(DiscoveryError, match="profiles without mechanism"):
        discover_dataset(dataset_root)


def test_summary_required_but_absent_raises(tmp_path: Path) -> None:
    root = write_dataset(
        tmp_path / "d", ["DENV1"], ["1st_run", "2nd_run"], with_summaries=False
    )
    with pytest.raises(DiscoveryError, match="missing Level-2 summaries"):
        discover_dataset(root, require_summaries=True)


# -- inconsistent replicate counts ------------------------------------------
def test_inconsistent_replicate_counts_raises(tmp_path: Path) -> None:
    # DENV1 gets 3 runs, DENV2 only 2 (remove its 3rd)
    root = write_dataset(tmp_path / "d", ["DENV1", "DENV2"], ["1st_run", "2nd_run", "3rd_run"])
    shutil.rmtree(root / "3rd_run" / "DENV2")
    with pytest.raises(DiscoveryError, match="inconsistent replicate counts"):
        discover_dataset(root)


def test_unequal_counts_allowed_when_flag_off(tmp_path: Path) -> None:
    root = write_dataset(tmp_path / "d", ["DENV1", "DENV2"], ["1st_run", "2nd_run", "3rd_run"])
    shutil.rmtree(root / "3rd_run" / "DENV2")
    ds = discover_dataset(root, enforce_equal_replicate_counts=False)
    counts = {s.serotype: s.n_replicates for s in ds.serotypes}
    assert counts == {"DENV1": 3, "DENV2": 2}


# -- malformed layout -------------------------------------------------------
def test_missing_analysis_output_dir_raises(dataset_root: Path) -> None:
    # rename a serotype's analysis_output so the correlations file is unreachable
    ao = dataset_root / "1st_run" / "DENV1" / "analysis_output"
    # keep the run-dir a run-dir (DENV2 still has analysis_output) but break DENV1
    for f in ao.iterdir():
        f.unlink()
    ao.rmdir()
    # DENV1 now has 2 replicates, DENV2 has 3 -> inconsistent
    with pytest.raises(DiscoveryError, match="inconsistent replicate counts"):
        discover_dataset(dataset_root)


def test_analysis_output_without_correlations_raises(dataset_root: Path) -> None:
    ao = dataset_root / "1st_run" / "DENV1" / "analysis_output"
    for f in ao.iterdir():
        f.unlink()  # empty analysis_output -> no correlations file
    with pytest.raises(DiscoveryError, match="no '.*correlations"):
        discover_dataset(dataset_root)


# -- duplicate serotype -----------------------------------------------------
def test_duplicate_profile_raises(dataset_root: Path) -> None:
    # place a second DENV1 profile at the root (another search dir)
    shutil.copy(
        dataset_root / "summaries" / "DENV1_profile.csv",
        dataset_root / "DENV1_profile.csv",
    )
    with pytest.raises(DiscoveryError, match="duplicate profile"):
        discover_dataset(dataset_root)


# -- partial datasets -------------------------------------------------------
def test_partial_dataset_summaries_only(tmp_path: Path) -> None:
    root = write_dataset(
        tmp_path / "d", ["DENV1", "DENV2"], ["1st_run"], with_summaries=True
    )
    # remove all replicates -> summaries-only
    for run in ("1st_run",):
        shutil.rmtree(root / run)
    ds = discover_dataset(root, require_replicates=False)
    assert ds.names == ("DENV1", "DENV2")
    assert all(s.n_replicates == 0 for s in ds.serotypes)


def test_partial_dataset_replicates_only(tmp_path: Path) -> None:
    root = write_dataset(
        tmp_path / "d", ["DENV1"], ["1st_run", "2nd_run"], with_summaries=False
    )
    ds = discover_dataset(root, require_summaries=False)
    assert ds.names == ("DENV1",)
    assert ds.serotypes[0].summary is None
    assert ds.serotypes[0].n_replicates == 2
