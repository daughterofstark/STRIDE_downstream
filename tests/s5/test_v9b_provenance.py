"""Adversarial provenance-layout checks for V9B."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from tests.s5 import v9b_provenance as provenance


def _repo() -> Path:
    return Path(__file__).resolve().parents[2]


def _copy_sources(destination: Path) -> None:
    repo = _repo()
    for relative in provenance.EXPECTED_SOURCE_HASHES:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / relative, target)


def test_full_descendant_and_detached_checkout(tmp_path: Path) -> None:
    full = provenance.verify_provenance(_repo())
    assert full["git_history_status"] == "frozen_anchor_ancestor_of_head"
    assert full["fully_verified"] is True

    detached = tmp_path / "detached"
    subprocess.run(["git", "clone", "-q", str(_repo()), str(detached)], check=True)
    subprocess.run(
        ["git", "checkout", "-q", "--detach", "HEAD"], cwd=detached, check=True
    )
    status = provenance.verify_provenance(detached)
    assert status["git_history_status"] == "frozen_anchor_ancestor_of_detached_head"
    assert status["fully_verified"] is True


def test_exact_anchor_checkout(tmp_path: Path) -> None:
    checkout = tmp_path / "anchor"
    subprocess.run(["git", "clone", "-q", str(_repo()), str(checkout)], check=True)
    subprocess.run(
        ["git", "checkout", "-q", "--detach", provenance.DOWNSTREAM_COMMIT],
        cwd=checkout,
        check=True,
    )
    status = provenance.verify_provenance(checkout)
    assert status["git_history_status"] == "exact_frozen_anchor"
    assert status["fully_verified"] is True


def test_shallow_checkout_reports_unavailable_ancestry(tmp_path: Path) -> None:
    shallow = tmp_path / "shallow"
    subprocess.run(
        ["git", "clone", "-q", "--depth", "1", f"file://{_repo()}", str(shallow)],
        check=True,
    )
    status = provenance.verify_provenance(shallow)
    assert status["git_history_status"] == "anchor_unavailable_in_shallow_checkout"
    assert status["source_integrity_status"] == "verified_against_frozen_hashes"
    assert status["fully_verified"] is False


def test_source_export_without_git_runs_with_source_hash_status(tmp_path: Path) -> None:
    export = tmp_path / "export"
    _copy_sources(export)
    status = provenance.verify_provenance(export)
    assert status["git_history_status"] == "git_metadata_unavailable"
    assert status["source_integrity_status"] == "verified_against_frozen_hashes"
    assert status["fully_verified"] is False


def test_unrelated_git_history_is_not_verified(tmp_path: Path) -> None:
    unrelated = tmp_path / "unrelated"
    _copy_sources(unrelated)
    subprocess.run(["git", "init", "-q", str(unrelated)], check=True)
    subprocess.run(["git", "config", "user.email", "audit@example.invalid"], cwd=unrelated, check=True)
    subprocess.run(["git", "config", "user.name", "Audit"], cwd=unrelated, check=True)
    subprocess.run(["git", "add", "."], cwd=unrelated, check=True)
    subprocess.run(["git", "commit", "-qm", "fabricated"], cwd=unrelated, check=True)
    status = provenance.provenance_status(unrelated)
    assert status["git_history_status"] == "incompatible_git_history"
    assert status["fully_verified"] is False
    with pytest.raises(RuntimeError, match="incompatible"):
        provenance.verify_provenance(unrelated)


def test_modified_helper_source_is_rejected(tmp_path: Path) -> None:
    export = tmp_path / "modified"
    _copy_sources(export)
    helper = export / "src/stride_s5/build/position_conservation.py"
    helper.write_text(helper.read_text() + "\n# modified\n")
    with pytest.raises(RuntimeError, match="hash mismatch"):
        provenance.verify_provenance(export)
