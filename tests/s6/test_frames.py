"""Unit tests for the frame-extraction helpers (``build/_frames``)."""
from __future__ import annotations

import pandas as pd

from stride_s6.build._frames import (
    chain_from_canon,
    per_replicate_effects_available,
    per_run_effects,
    runs_with_effects,
    serotype_effect_matrix,
)


def test_chain_from_canon() -> None:
    assert chain_from_canon("NS3:51") == "NS3"
    assert chain_from_canon("NS2B:-1") == "NS2B"
    assert chain_from_canon("NS3") == "NS3"


def test_per_run_effects_none_returns_empty() -> None:
    out = per_run_effects(None)
    assert out.empty
    assert list(out.columns) == ["serotype", "canon_label", "replicate_index", "r"]


def test_per_run_effects_drops_nan_and_dedups() -> None:
    df = pd.DataFrame(
        {
            "serotype": ["DENV1", "DENV1", "DENV1"],
            "replicate_index": [1, 1, 2],
            "canon_label": ["NS3:51", "NS3:51", "NS3:51"],
            "r": [0.4, float("nan"), 0.5],
        }
    )
    out = per_run_effects(df)
    # the NaN row is dropped; the duplicate (serotype, run, position) collapses
    assert len(out) == 2
    assert set(out["replicate_index"]) == {1, 2}


def test_per_run_effects_does_not_mutate_input(replicate_table: pd.DataFrame) -> None:
    before = replicate_table.copy()
    per_run_effects(replicate_table)
    pd.testing.assert_frame_equal(replicate_table, before)


def test_runs_with_effects(effects: pd.DataFrame) -> None:
    assert runs_with_effects(effects, "DENV1") == [1, 2, 3]
    assert runs_with_effects(effects, "NOPE") == []


def test_per_replicate_effects_available_true_false(effects: pd.DataFrame) -> None:
    assert per_replicate_effects_available(effects) is True
    assert per_replicate_effects_available(per_run_effects(None)) is False


def test_serotype_effect_matrix_complete_cases(effects: pd.DataFrame) -> None:
    positions, runs, matrix = serotype_effect_matrix(effects, "DENV1")
    assert runs == [1, 2, 3]
    assert len(positions) == 6  # all six positions complete for DENV1
    assert all(len(row) == 3 for row in matrix)


def test_serotype_effect_matrix_drops_incomplete(effects: pd.DataFrame) -> None:
    # DENV4 is missing NS3:250 in run 3 → only five complete-case positions
    positions, runs, matrix = serotype_effect_matrix(effects, "DENV4")
    assert runs == [1, 2, 3]
    assert len(positions) == 5
    assert "NS3:250" not in positions
