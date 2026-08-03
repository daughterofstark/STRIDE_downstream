"""V9B deterministic validation of presence-aware S5 recurrence semantics.

The eight-row fixture is the exact truth table registered in Phase 3A. It is
converted to minimal S0/S1A-shaped in-memory frames and evaluated by the
production ``build_position_conservation`` helper. Explicit all-four, absence,
and present-below-cutoff fields are added only to the validation result.

Run this module with ``--write-artifacts`` to write the four deterministic V9B
artifacts and their manifest under ``validation_artifacts/``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from stride_s5.build import build_position_conservation


DOWNSTREAM_COMMIT = "4388dbcefaaa1672b63cc7bbce357e1b4ec40b71"
DATASETS = ("D1", "D2", "D3", "D4")
N_TOTAL = 4
RHO_STAR = 0.5
ROW_ORDER = ("P4A", "P4B", "P3A", "P3B", "P2A", "P2B", "P1A", "P1B")
ARTIFACT_DIR = Path(__file__).resolve().parents[2] / "validation_artifacts"
FIXTURE_NAME = "v9b_partial_presence_fixture.csv"
EXPECTED_NAME = "v9b_partial_presence_expected.csv"
RESULTS_NAME = "v9b_partial_presence_results.json"
MANIFEST_NAME = "v9b_manifest.json"

FIXTURE_ROWS = (
    ("P4A", "D1;D2;D3;D4", "D1;D2;D3;D4", "D1;D2;D3;D4"),
    ("P4B", "D1;D2;D3;D4", "D1;D2", "D1;D2"),
    ("P3A", "D1;D2;D3", "D1;D2;D3", "D1;D2;D3"),
    ("P3B", "D1;D2;D3", "D1", "D1"),
    ("P2A", "D1;D2", "D1;D2", "D1;D2"),
    ("P2B", "D1;D2", "", ""),
    ("P1A", "D1", "D1", "D1"),
    ("P1B", "D1", "", ""),
)

EXPECTED_ROWS = (
    ("P4A", 4, 4, 4, 1.0, "reproducible_all", 4, False, True, True, 0, 0),
    ("P4B", 4, 4, 2, 0.5, "reproducible_some", 2, True, True, False, 0, 2),
    ("P3A", 4, 3, 3, 1.0, "reproducible_all", 3, True, False, False, 1, 0),
    ("P3B", 4, 3, 1, 0.333333, "reproducible_some", 1, True, False, False, 1, 2),
    ("P2A", 4, 2, 2, 1.0, "reproducible_all", 2, True, False, False, 2, 0),
    ("P2B", 4, 2, 0, 0.0, "reproducible_none", 0, False, False, False, 2, 2),
    ("P1A", 4, 1, 1, 1.0, "reproducible_all", 1, True, False, False, 3, 0),
    ("P1B", 4, 1, 0, 0.0, "reproducible_none", 0, False, False, False, 3, 1),
)

FIXTURE_COLUMNS = (
    "coordinate_key",
    "datasets_present",
    "datasets_crossing_provisional_0_5",
    "datasets_signed_crossing",
)
SEMANTIC_COLUMNS = (
    "coordinate_key",
    "n_total",
    "n_present",
    "n_reproducible",
    "frac_reproducible",
    "conservation_class",
    "n_signed_reproducible",
    "is_serotype_divergent",
    "present_in_all_four",
    "crossing_all_four",
    "n_absent",
    "n_present_below_cutoff",
)


def fixture_frame() -> pd.DataFrame:
    """Return the exact registered fixture in its fixed row order."""
    return pd.DataFrame(FIXTURE_ROWS, columns=FIXTURE_COLUMNS)


def expected_frame() -> pd.DataFrame:
    """Return the exact registered semantic truth table."""
    return pd.DataFrame(EXPECTED_ROWS, columns=SEMANTIC_COLUMNS)


def _dataset_set(value: object) -> set[str]:
    if value is None or pd.isna(value) or str(value) == "":
        return set()
    return {part for part in str(value).split(";") if part}


def build_input_frames(fixture: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build minimal S0/S1A-shaped frames from the neutral fixture."""
    stride_rows: list[dict[str, object]] = []
    presence_rows: list[dict[str, object]] = []
    for row in fixture.itertuples(index=False):
        key = str(row.coordinate_key)
        present = _dataset_set(row.datasets_present)
        crossing = _dataset_set(row.datasets_crossing_provisional_0_5)
        signed = _dataset_set(row.datasets_signed_crossing)
        assert signed <= crossing <= present <= set(DATASETS)
        presence_rows.append({
            "canon_label": key,
            "n_serotypes": len(present),
            "serotypes_present": sorted(present),
            "serotypes_absent": sorted(set(DATASETS) - present),
            "in_all_serotypes": present == set(DATASETS),
            "in_any_serotype": bool(present),
            "chain": "C1",
            "domain": "G1",
        })
        for dataset in DATASETS:
            if dataset not in present:
                continue
            crosses = dataset in crossing
            is_signed = dataset in signed
            stride_rows.append(_stride_row(
                dataset=dataset,
                key=key,
                scale_level="residue",
                scale_index=0,
                rho=0.75 if crosses else 0.25,
                gated=crosses,
                direction="increase" if is_signed else "mixed",
            ))
            if not crosses:
                stride_rows.append(_stride_row(
                    dataset=dataset,
                    key=key,
                    scale_level="domain",
                    scale_index=3,
                    rho=0.75,
                    gated=True,
                    direction="mixed",
                ))
    stride = pd.DataFrame.from_records(stride_rows).sort_values(
        ["serotype", "canon_label", "scale_index"], kind="mergesort"
    ).reset_index(drop=True)
    presence = pd.DataFrame.from_records(presence_rows)
    presence["_order"] = presence["canon_label"].map(
        {key: index for index, key in enumerate(ROW_ORDER)}
    )
    presence = presence.sort_values("_order").drop(columns="_order").reset_index(
        drop=True
    )
    return stride, presence


def _stride_row(
    *, dataset: str, key: str, scale_level: str, scale_index: int,
    rho: float, gated: bool, direction: str,
) -> dict[str, object]:
    return {
        "serotype": dataset,
        "canon_label": key,
        "scale_level": scale_level,
        "scale_index": scale_index,
        "region_id": f"root/C1/G1/{key}/{scale_level}",
        "rho": rho,
        "gated": gated,
        "beta": 1.0,
        "beta_se": 0.1,
        "tau2": 0.2,
        "sigma2_bar": 0.3,
        "h_chain": "C1",
        "h_domain": "G1",
        "is_gated_scale": gated,
        "mech_direction": direction,
    }


def build_results(fixture: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the production S5 builder and assemble explicit validation fields."""
    stride, presence = build_input_frames(fixture)
    production = build_position_conservation(stride, presence, rho_star=RHO_STAR)
    indexed = production.set_index("canon_label")
    rows = []
    for key in ROW_ORDER:
        row = indexed.loc[key]
        n_total = int(row["n_serotypes_total"])
        n_present = int(row["n_serotypes_present"])
        n_reproducible = int(row["n_serotypes_reproducible"])
        rows.append({
            "coordinate_key": key,
            "n_total": n_total,
            "n_present": n_present,
            "n_reproducible": n_reproducible,
            "frac_reproducible": float(row["frac_reproducible"]),
            "conservation_class": str(row["conservation_class"]),
            "n_signed_reproducible": int(
                row["n_serotypes_signed_reproducible"]
            ),
            "is_serotype_divergent": bool(row["is_serotype_divergent"]),
            "present_in_all_four": n_present == N_TOTAL,
            "crossing_all_four": (
                n_total == N_TOTAL
                and n_present == N_TOTAL
                and n_reproducible == N_TOTAL
            ),
            "n_absent": n_total - n_present,
            "n_present_below_cutoff": n_present - n_reproducible,
        })
    return pd.DataFrame.from_records(rows, columns=SEMANTIC_COLUMNS), production


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False
        )
        + "\n"
    ).encode("ascii")


def _csv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False, lineterminator="\n").encode("ascii")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def artifact_payloads() -> dict[str, bytes]:
    """Render all deterministic V9B artifacts in memory."""
    fixture = fixture_frame()
    expected = expected_frame()
    observed, production = build_results(fixture)
    result_document = {
        "schema_version": "v9b.results.1",
        "purpose": "deterministic_presence_and_recurrence_semantics_validation",
        "rho_star": RHO_STAR,
        "rho_star_status": "provisional_operational_reference",
        "datasets": list(DATASETS),
        "row_order": list(ROW_ORDER),
        "production_helper": (
            "stride_s5.build.build_position_conservation"
        ),
        "production_output_columns": list(production.columns),
        "rows": observed.to_dict(orient="records"),
        "semantic_findings": {
            "reproducible_all_means": (
                "crossing in every dataset where present"
            ),
            "all_four_requires": (
                "n_total == 4 and n_present == 4 and n_reproducible == 4"
            ),
            "legacy_divergence_rule": (
                "0 < n_signed_reproducible < n_total"
            ),
            "p3a_divergence_due_to_absence_is_true": True,
            "absence_and_present_below_cutoff_are_distinct": True,
        },
        "interpretation_boundaries": [
            "does_not_calibrate_rho_star",
            "does_not_test_biology_or_homology",
            "does_not_validate_significance_or_fdr",
            "does_not_modify_frozen_s5_outputs",
        ],
    }
    payloads = {
        FIXTURE_NAME: _csv_bytes(fixture),
        EXPECTED_NAME: _csv_bytes(expected),
        RESULTS_NAME: _json_bytes(result_document),
    }
    repo = Path(__file__).resolve().parents[2]
    observed_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repo, text=True
    ).strip()
    if observed_commit != DOWNSTREAM_COMMIT:
        raise RuntimeError(
            f"V9B requires downstream {DOWNSTREAM_COMMIT}, "
            f"observed {observed_commit}"
        )
    manifest = {
        "schema_version": "v9b.manifest.1",
        "purpose": "deterministic_presence_and_recurrence_semantics_validation",
        "source": {
            "downstream_commit_required": DOWNSTREAM_COMMIT,
            "downstream_commit_observed": observed_commit,
            "test_module_sha256": _sha256(Path(__file__).read_bytes()),
        },
        "environment": {
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "packages": {
                name: version(name)
                for name in ("pandas", "pyarrow", "pydantic", "pytest")
            },
        },
        "design": {
            "datasets": list(DATASETS),
            "n_total": N_TOTAL,
            "coordinate_keys": list(ROW_ORDER),
            "rho_star": RHO_STAR,
            "stochastic": False,
        },
        "artifact_sha256": {
            name: _sha256(data) for name, data in sorted(payloads.items())
        },
    }
    payloads[MANIFEST_NAME] = _json_bytes(manifest)
    return payloads


def write_artifacts(output_dir: Path = ARTIFACT_DIR) -> dict[str, str]:
    """Write the four deterministic artifacts and return their SHA-256 values."""
    output_dir.mkdir(parents=True, exist_ok=True)
    payloads = artifact_payloads()
    for name, data in payloads.items():
        (output_dir / name).write_bytes(data)
    return {name: _sha256(data) for name, data in sorted(payloads.items())}


def test_all_eight_rows_match_registered_truth_table_exactly() -> None:
    observed, _ = build_results(fixture_frame())
    pd.testing.assert_frame_equal(observed, expected_frame(), check_exact=True)


def test_all_four_and_legacy_all_semantics() -> None:
    observed, _ = build_results(fixture_frame())
    all_four = observed.loc[observed["crossing_all_four"], "coordinate_key"].tolist()
    assert all_four == ["P4A"]
    legacy_all = observed.loc[
        observed["conservation_class"] == "reproducible_all", "coordinate_key"
    ].tolist()
    assert legacy_all == ["P4A", "P3A", "P2A", "P1A"]
    for key in ("P3A", "P2A", "P1A"):
        row = observed.set_index("coordinate_key").loc[key]
        assert not bool(row["present_in_all_four"])
        assert not bool(row["crossing_all_four"])


def test_absence_and_present_below_cutoff_remain_distinct() -> None:
    observed, _ = build_results(fixture_frame())
    rows = observed.set_index("coordinate_key")
    assert int(rows.loc["P3A", "n_absent"]) == 1
    assert int(rows.loc["P3A", "n_present_below_cutoff"]) == 0
    assert int(rows.loc["P4B", "n_absent"]) == 0
    assert int(rows.loc["P4B", "n_present_below_cutoff"]) == 2


def test_legacy_divergence_rule_exposes_p3a_absence_case() -> None:
    observed, _ = build_results(fixture_frame())
    p3a = observed.set_index("coordinate_key").loc["P3A"]
    assert int(p3a["n_signed_reproducible"]) == 3
    assert int(p3a["n_present"]) == 3
    assert int(p3a["n_total"]) == 4
    assert bool(p3a["is_serotype_divergent"])
    assert int(p3a["n_absent"]) == 1
    assert int(p3a["n_present_below_cutoff"]) == 0


def test_row_order_and_artifact_serialization_are_deterministic(
    tmp_path: Path,
) -> None:
    first = artifact_payloads()
    second = artifact_payloads()
    assert first == second
    hashes = write_artifacts(tmp_path)
    assert list(build_results(fixture_frame())[0]["coordinate_key"]) == list(
        ROW_ORDER
    )
    for name, data in first.items():
        assert (tmp_path / name).read_bytes() == data
        assert hashes[name] == _sha256(data)


def test_committed_artifacts_match_current_execution() -> None:
    payloads = artifact_payloads()
    assert set(payloads) == {
        FIXTURE_NAME, EXPECTED_NAME, RESULTS_NAME, MANIFEST_NAME
    }
    for name, data in payloads.items():
        assert (ARTIFACT_DIR / name).read_bytes() == data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-artifacts", action="store_true")
    args = parser.parse_args()
    if not args.write_artifacts:
        parser.error("--write-artifacts is required outside pytest")
    print(json.dumps(write_artifacts(), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
