"""Unit tests for the deterministic S3 hierarchy-reduction helpers."""
from __future__ import annotations

import math

from stride_s3.build._curves import (
    audit_monotonicity,
    is_distributed_effect,
    scale_tier,
    step_gains,
)
from stride_s3.models.schema import TIER_EXPLORATORY, TIER_LICENSED


# ---------------------------------------------------------------------------
# scale_tier
# ---------------------------------------------------------------------------
def test_scale_tier_residue_and_finer_are_exploratory() -> None:
    assert scale_tier("residue") == TIER_EXPLORATORY
    assert scale_tier("secondary_structure") == TIER_EXPLORATORY
    assert scale_tier("motif") == TIER_EXPLORATORY


def test_scale_tier_domain_and_coarser_are_licensed() -> None:
    assert scale_tier("domain") == TIER_LICENSED
    assert scale_tier("chain") == TIER_LICENSED
    assert scale_tier("protein") == TIER_LICENSED
    assert scale_tier("complex") == TIER_LICENSED


def test_scale_tier_unknown_is_exploratory() -> None:
    assert scale_tier("nonsense") == TIER_EXPLORATORY


# ---------------------------------------------------------------------------
# step_gains
# ---------------------------------------------------------------------------
def test_step_gains_first_is_nan() -> None:
    gains = step_gains([0.3, 0.35, 0.7])
    assert math.isnan(gains[0])
    assert gains[1] == 0.35 - 0.30
    assert round(gains[2], 6) == round(0.70 - 0.35, 6)


def test_step_gains_length_matches_input() -> None:
    curve = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    assert len(step_gains(curve)) == len(curve)


def test_step_gains_can_be_negative() -> None:
    gains = step_gains([0.5, 0.4])
    assert round(gains[1], 6) == -0.1


# ---------------------------------------------------------------------------
# audit_monotonicity
# ---------------------------------------------------------------------------
def test_monotone_curve_has_no_violations() -> None:
    is_mono, n_viol, max_dec, first = audit_monotonicity(
        [0.3, 0.35, 0.4, 0.7, 0.75, 0.8, 0.85]
    )
    assert is_mono
    assert n_viol == 0
    assert max_dec == 0.0
    assert first == -1


def test_single_dip_is_flagged() -> None:
    # dip at index 2 (0.45 -> 0.40): violation recorded at the finer index (1)
    is_mono, n_viol, max_dec, first = audit_monotonicity(
        [0.42, 0.45, 0.40, 0.72, 0.78, 0.82, 0.90]
    )
    assert not is_mono
    assert n_viol == 1
    assert round(max_dec, 6) == round(0.05, 6)
    assert first == 1


def test_multiple_violations_counted_and_max_tracked() -> None:
    is_mono, n_viol, max_dec, first = audit_monotonicity(
        [0.9, 0.5, 0.8, 0.2, 0.3, 0.4, 0.5]
    )
    assert not is_mono
    assert n_viol == 2  # 0.9->0.5 and 0.8->0.2
    assert round(max_dec, 6) == round(0.6, 6)  # 0.8 -> 0.2
    assert first == 0


def test_monotonicity_skips_nan() -> None:
    is_mono, n_viol, _max, _first = audit_monotonicity(
        [0.3, float("nan"), 0.5, 0.7, 0.75, 0.8, 0.85]
    )
    assert is_mono
    assert n_viol == 0


# ---------------------------------------------------------------------------
# is_distributed_effect
# ---------------------------------------------------------------------------
def test_distributed_when_residue_below_and_domain_above_threshold() -> None:
    assert is_distributed_effect(0.30, 0.70, 0.5)


def test_not_distributed_when_residue_already_reproducible() -> None:
    assert not is_distributed_effect(0.80, 0.88, 0.5)


def test_not_distributed_when_domain_below_threshold() -> None:
    assert not is_distributed_effect(0.30, 0.45, 0.5)


def test_distributed_boundary_domain_equal_threshold() -> None:
    # domain exactly at ρ* counts (>=), residue strictly below (<)
    assert is_distributed_effect(0.49, 0.50, 0.5)


def test_distributed_nan_is_false() -> None:
    assert not is_distributed_effect(float("nan"), 0.7, 0.5)
    assert not is_distributed_effect(0.3, float("nan"), 0.5)
