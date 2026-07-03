"""Unit tests for the deterministic S2 screen/regate helpers."""
from __future__ import annotations

import math

from stride_s2.build._screens import (
    ci_excludes_zero,
    is_residue_scale,
    is_signed,
    normalise_band,
    passes_signed_screen,
    regate_scale,
    scale_tier,
)
from stride_s2.models.schema import (
    SCALE_UNRESOLVED,
    TIER_EXPLORATORY,
    TIER_LICENSED,
)


# ---------------------------------------------------------------------------
# regate_scale
# ---------------------------------------------------------------------------
def test_regate_picks_finest_scale_meeting_threshold() -> None:
    rho = {
        "residue": 0.80,
        "secondary_structure": 0.82,
        "motif": 0.85,
        "domain": 0.88,
        "chain": 0.90,
        "protein": 0.92,
        "complex": 0.95,
    }
    assert regate_scale(rho, 0.5) == ("residue", 0)
    # raising ρ* past residue/SS pushes the gate to motif
    assert regate_scale(rho, 0.85) == ("motif", 2)


def test_regate_gates_coarse_when_fine_scales_below_threshold() -> None:
    rho = {
        "residue": 0.30,
        "secondary_structure": 0.35,
        "motif": 0.40,
        "domain": 0.70,
        "chain": 0.75,
        "protein": 0.80,
        "complex": 0.85,
    }
    assert regate_scale(rho, 0.5) == ("domain", 3)


def test_regate_unresolved_when_no_scale_meets_threshold() -> None:
    rho = {"residue": 0.3, "domain": 0.7, "complex": 0.85}
    assert regate_scale(rho, 0.95) == (SCALE_UNRESOLVED, -1)


def test_regate_skips_nan() -> None:
    rho = {"residue": float("nan"), "domain": 0.7}
    assert regate_scale(rho, 0.5) == ("domain", 3)


# ---------------------------------------------------------------------------
# scale_tier
# ---------------------------------------------------------------------------
def test_scale_tier_residue_is_exploratory() -> None:
    assert scale_tier("residue") == TIER_EXPLORATORY
    assert scale_tier("secondary_structure") == TIER_EXPLORATORY
    assert scale_tier("motif") == TIER_EXPLORATORY


def test_scale_tier_domain_and_coarser_is_licensed() -> None:
    assert scale_tier("domain") == TIER_LICENSED
    assert scale_tier("chain") == TIER_LICENSED
    assert scale_tier("complex") == TIER_LICENSED


def test_scale_tier_unresolved_is_exploratory() -> None:
    assert scale_tier(SCALE_UNRESOLVED) == TIER_EXPLORATORY


def test_is_residue_scale() -> None:
    assert is_residue_scale("residue")
    assert not is_residue_scale("domain")


# ---------------------------------------------------------------------------
# signed / CI screens
# ---------------------------------------------------------------------------
def test_is_signed() -> None:
    assert is_signed("increase")
    assert is_signed("decrease")
    assert not is_signed("mixed")
    assert not is_signed(None)


def test_ci_excludes_zero_positive_interval() -> None:
    assert ci_excludes_zero(0.02, 0.22)


def test_ci_excludes_zero_negative_interval() -> None:
    assert ci_excludes_zero(-0.30, -0.05)


def test_ci_touching_zero_does_not_exclude() -> None:
    assert not ci_excludes_zero(-0.15, 0.0)
    assert not ci_excludes_zero(0.0, 0.22)


def test_ci_spanning_zero_does_not_exclude() -> None:
    assert not ci_excludes_zero(-0.1, 0.1)


def test_ci_none_or_nan_does_not_exclude() -> None:
    assert not ci_excludes_zero(None, 0.2)
    assert not ci_excludes_zero(0.1, None)
    assert not ci_excludes_zero(float("nan"), 0.2)


def test_ci_reversed_bounds_are_tolerated() -> None:
    # a CI provided high-then-low still excludes zero if both are same-sign
    assert ci_excludes_zero(0.22, 0.02)


# ---------------------------------------------------------------------------
# passes_signed_screen
# ---------------------------------------------------------------------------
def test_passes_screen_all_conditions_met() -> None:
    assert passes_signed_screen("increase", 0.02, 0.22, rho=0.8, rho_star=0.5)


def test_passes_screen_fails_when_mixed() -> None:
    assert not passes_signed_screen("mixed", None, None, rho=0.8, rho_star=0.5)


def test_passes_screen_fails_when_ci_touches_zero() -> None:
    assert not passes_signed_screen(
        "decrease", -0.15, 0.0, rho=0.8, rho_star=0.5
    )


def test_passes_screen_fails_when_rho_below_threshold() -> None:
    assert not passes_signed_screen(
        "increase", 0.02, 0.22, rho=0.7, rho_star=0.85
    )


def test_passes_screen_fails_when_rho_nan() -> None:
    assert not passes_signed_screen(
        "increase", 0.02, 0.22, rho=float("nan"), rho_star=0.5
    )


# ---------------------------------------------------------------------------
# normalise_band
# ---------------------------------------------------------------------------
def test_normalise_band_rounds_dedupes_and_sorts() -> None:
    band = normalise_band([0.9, 0.5, 0.5, 0.70001], decimals=4)
    assert band == (0.5, 0.7, 0.9)


def test_normalise_band_empty() -> None:
    assert normalise_band([], decimals=4) == ()


def test_normalise_band_is_deterministic() -> None:
    a = normalise_band([0.6, 0.5, 0.7], decimals=4)
    b = normalise_band([0.7, 0.6, 0.5], decimals=4)
    assert a == b


def test_normalise_band_floats_are_finite() -> None:
    for value in normalise_band([0.5, 0.55, 0.6], decimals=4):
        assert math.isfinite(value)
