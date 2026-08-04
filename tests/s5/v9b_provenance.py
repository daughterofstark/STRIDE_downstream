"""Source-integrity and Git-history provenance for V9B validation."""
from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

DOWNSTREAM_COMMIT = "4388dbcefaaa1672b63cc7bbce357e1b4ec40b71"
EXPECTED_SOURCE_HASHES = {
    "src/stride_s5/build/__init__.py": "d0779f5ed1b0bd4dab8a366ca08309f0ae53f2afc2ecef4e0919bacbd0979362",
    "src/stride_s5/build/_classify.py": "b9e0f5cc386e4bfda23ac180c1a16ec07bcb062e81e8c24e68a603fcf95e341c",
    "src/stride_s5/build/_frames.py": "be9c3c5e9f52fbe7dea5cd66c7e8d1e249db7b170529f39ba6eee6997a6edbbd",
    "src/stride_s5/build/position_conservation.py": "4c1d6efa5701995a41d5761a0928f052568f9f3f9c208cecf5888d42107457ad",
    "src/stride_s5/models/errors.py": "9abe6d84806d8a9bed44f716bd4e6e05165aa852d4079a43664ca25821210165",
    "src/stride_s5/models/schema.py": "1f7f656b3b8a08313d017c0ab0d8548d4ab5386dbff5b4a4d8ffb2d4e9834e75",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *args], cwd=repo, text=True, capture_output=True, check=False
    )


def _git_value(repo: Path, *args: str) -> str | None:
    result = _git(repo, *args)
    return result.stdout.strip() if result.returncode == 0 else None


def verify_source_hashes(repo_root: str | Path) -> dict[str, str]:
    repo = Path(repo_root).resolve()
    observed: dict[str, str] = {}
    missing: list[str] = []
    for relative, expected in EXPECTED_SOURCE_HASHES.items():
        path = repo / relative
        if not path.is_file():
            missing.append(relative)
            continue
        digest = _sha256(path)
        observed[relative] = digest
        if digest != expected:
            raise RuntimeError(
                f"V9B production helper hash mismatch for {relative}: "
                f"expected {expected}, observed {digest}"
            )
    if missing:
        raise RuntimeError(f"missing V9B production helper sources: {missing}")
    return observed


def provenance_status(repo_root: str | Path) -> dict[str, object]:
    """Return explicit source-integrity and Git-history evidence."""
    repo = Path(repo_root).resolve()
    hashes = verify_source_hashes(repo)
    head = _git_value(repo, "rev-parse", "HEAD")
    if head is None:
        return {
            "required_downstream_anchor": DOWNSTREAM_COMMIT,
            "observed_head": None,
            "head_detached": None,
            "git_history_status": "git_metadata_unavailable",
            "source_integrity_status": "verified_against_frozen_hashes",
            "fully_verified": False,
            "production_source_hashes": hashes,
        }

    detached = _git(repo, "symbolic-ref", "-q", "HEAD").returncode != 0
    anchor_exists = _git(
        repo, "cat-file", "-e", f"{DOWNSTREAM_COMMIT}^{{commit}}"
    ).returncode == 0
    if anchor_exists:
        ancestor = _git(
            repo, "merge-base", "--is-ancestor", DOWNSTREAM_COMMIT, head
        ).returncode == 0
        if not ancestor:
            history = "incompatible_git_history"
            fully_verified = False
        else:
            history = "exact_frozen_anchor" if head == DOWNSTREAM_COMMIT else (
                "frozen_anchor_ancestor_of_detached_head" if detached
                else "frozen_anchor_ancestor_of_head"
            )
            fully_verified = True
    else:
        shallow = _git_value(repo, "rev-parse", "--is-shallow-repository")
        if shallow == "true":
            history = "anchor_unavailable_in_shallow_checkout"
        else:
            history = "incompatible_git_history"
        fully_verified = False
    return {
        "required_downstream_anchor": DOWNSTREAM_COMMIT,
        "observed_head": head,
        "head_detached": detached,
        "git_history_status": history,
        "source_integrity_status": "verified_against_frozen_hashes",
        "fully_verified": fully_verified,
        "production_source_hashes": hashes,
    }


def verify_provenance(repo_root: str | Path) -> dict[str, object]:
    """Require compatible Git history when Git metadata is available."""
    status = provenance_status(repo_root)
    if status["git_history_status"] == "incompatible_git_history":
        raise RuntimeError(
            f"V9B Git history is incompatible with {DOWNSTREAM_COMMIT}"
        )
    return status
