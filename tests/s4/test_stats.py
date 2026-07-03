"""Unit tests for the deterministic S4 statistics helpers."""
from __future__ import annotations

import math

from stride_s4.build._stats import (
    benjamini_hochberg,
    ci_excludes_zero,
    classify_variance_regime,
    normal_cdf,
    two_sided_wald_p,
    variance_fractions,
    wald_z,
    weighted_mean_se,
)
from stride_s4.models.schema import (
    VARIANCE_REGIME_BALANCED,
    VARIANCE_REGIME_REPLICATE,
    VARIANCE_REGIME_SAMPLING,
    VARIANCE_REGIME_UNDEFINED,
)


# ---------------------------------------------------------------------------
# normal_cdf / Wald
# ---------------------------------------------------------------------------
def test_normal_cdf_known_points() -> None:
    assert abs(normal_cdf(0.0) - 0.5) < 1e-12
    assert abs(normal_cdf(1.96) - 0.975) < 1e-3
    assert abs(normal_cdf(-1.96) - 0.025) < 1e-3


def test_wald_p_known_z() -> None:
    # z = 0.12 / 0.05 = 2.4 → two-sided p ≈ 0.016395
    p = two_sided_wald_p(0.12, 0.05)
    assert abs(p - 0.016395) < 1e-5


def test_wald_p_at_196_is_005() -> None:
    # a z of ~1.96 gives p ~ 0.05
    p = two_sided_wald_p(1.96, 1.0)
    assert abs(p - 0.05) < 1e-3


def test_wald_z_and_p_nan_when_se_nonpositive() -> None:
    assert math.isnan(wald_z(0.1, 0.0))
    assert math.isnan(two_sided_wald_p(0.1, 0.0))
    assert math.isnan(two_sided_wald_p(float("nan"), 0.05))


# ---------------------------------------------------------------------------
# ci_excludes_zero
# ---------------------------------------------------------------------------
def test_ci_excludes_zero_positive_interval() -> None:
    assert ci_excludes_zero(0.02, 0.22)


def test_ci_excludes_zero_negative_interval() -> None:
    assert ci_excludes_zero(-0.30, -0.02)


def test_ci_touching_zero_does_not_exclude() -> None:
    assert not ci_excludes_zero(-0.01, 0.21)
    assert not ci_excludes_zero(0.0, 0.2)


def test_ci_nan_bound_is_false() -> None:
    assert not ci_excludes_zero(float("nan"), 0.2)
    assert not ci_excludes_zero(0.1, float("nan"))


# ---------------------------------------------------------------------------
# benjamini_hochberg
# ---------------------------------------------------------------------------
def test_bh_known_small_case() -> None:
    # classic worked example
    raw = [0.01, 0.02, 0.03, 0.04, 0.05]
    adj = benjamini_hochberg(raw)
    # BH: p*m/rank then monotone from the top
    #  rank5: 0.05*5/5 = 0.05
    #  rank4: 0.04*5/4 = 0.05
    #  rank3: 0.03*5/3 = 0.05
    #  rank2: 0.02*5/2 = 0.05
    #  rank1: 0.01*5/1 = 0.05
    for a in adj:
        assert abs(a - 0.05) < 1e-9


def test_bh_monotone_and_clamped() -> None:
    raw = [0.9, 0.001, 0.5]
    adj = benjamini_hochberg(raw)
    # smallest raw (0.001, rank1): 0.001*3/1 = 0.003
    # 0.5 (rank2): 0.5*3/2 = 0.75
    # 0.9 (rank3): 0.9*3/3 = 0.9
    assert abs(adj[1] - 0.003) < 1e-9
    assert abs(adj[2] - 0.75) < 1e-9
    assert abs(adj[0] - 0.9) < 1e-9
    assert all(0.0 <= a <= 1.0 for a in adj)


def test_bh_preserves_input_order() -> None:
    raw = [0.04, 0.01]
    adj = benjamini_hochberg(raw)
    # rank1 = 0.01 (index 1): 0.01*2/1 = 0.02
    # rank2 = 0.04 (index 0): 0.04*2/2 = 0.04
    assert abs(adj[0] - 0.04) < 1e-9
    assert abs(adj[1] - 0.02) < 1e-9


def test_bh_skips_nan() -> None:
    raw = [0.01, float("nan"), 0.02]
    adj = benjamini_hochberg(raw)
    assert math.isnan(adj[1])
    # family size m = 2: rank1 0.01→0.02, rank2 0.02→0.02
    assert abs(adj[0] - 0.02) < 1e-9
    assert abs(adj[2] - 0.02) < 1e-9


def test_bh_empty_and_all_nan() -> None:
    assert benjamini_hochberg([]) == []
    out = benjamini_hochberg([float("nan"), float("nan")])
    assert all(math.isnan(v) for v in out)


# ---------------------------------------------------------------------------
# variance_fractions / regime
# ---------------------------------------------------------------------------
def test_variance_fractions_basic() -> None:
    total, ft, fs, ratio, regime = variance_fractions(0.2, 0.3)
    assert total == 0.5
    assert abs(ft - 0.4) < 1e-12
    assert abs(fs - 0.6) < 1e-12
    assert abs(ratio - (0.2 / 0.3)) < 1e-12
    assert regime == VARIANCE_REGIME_SAMPLING


def test_variance_fractions_replicate_dominated() -> None:
    _t, ft, _fs, _r, regime = variance_fractions(0.42, 0.18)
    assert abs(ft - 0.7) < 1e-12
    assert regime == VARIANCE_REGIME_REPLICATE


def test_variance_fractions_balanced() -> None:
    _t, ft, _fs, _r, regime = variance_fractions(0.25, 0.25)
    assert abs(ft - 0.5) < 1e-12
    assert regime == VARIANCE_REGIME_BALANCED


def test_variance_fractions_zero_total_undefined() -> None:
    total, ft, fs, ratio, regime = variance_fractions(0.0, 0.0)
    assert total == 0.0
    assert math.isnan(ft) and math.isnan(fs) and math.isnan(ratio)
    assert regime == VARIANCE_REGIME_UNDEFINED


def test_variance_fractions_zero_sigma_ratio_nan() -> None:
    _t, ft, _fs, ratio, regime = variance_fractions(0.5, 0.0)
    assert abs(ft - 1.0) < 1e-12
    assert math.isnan(ratio)  # division by zero σ̄²
    assert regime == VARIANCE_REGIME_REPLICATE


def test_classify_regime_boundaries() -> None:
    assert classify_variance_regime(0.6) == VARIANCE_REGIME_REPLICATE
    assert classify_variance_regime(0.4) == VARIANCE_REGIME_SAMPLING
    assert classify_variance_regime(0.5) == VARIANCE_REGIME_BALANCED
    assert classify_variance_regime(float("nan")) == VARIANCE_REGIME_UNDEFINED


# ---------------------------------------------------------------------------
# weighted_mean_se
# ---------------------------------------------------------------------------
def test_weighted_mean_equal_ses_is_plain_mean() -> None:
    wm, wse, um = weighted_mean_se([0.1, 0.3], [0.2, 0.2])
    assert abs(wm - 0.2) < 1e-12
    assert abs(um - 0.2) < 1e-12
    # weighted se = sqrt(1 / (1/0.04 + 1/0.04)) = sqrt(1/50)
    assert abs(wse - math.sqrt(1.0 / 50.0)) < 1e-12


def test_weighted_mean_favours_precise_estimate() -> None:
    # second estimate has a much smaller SE → weighted mean pulled toward it
    wm, _wse, um = weighted_mean_se([0.0, 1.0], [1.0, 0.1])
    assert wm > 0.9
    assert abs(um - 0.5) < 1e-12


def test_weighted_mean_skips_bad_pairs() -> None:
    wm, _wse, um = weighted_mean_se(
        [0.2, float("nan"), 0.4], [0.1, 0.1, 0.0]
    )
    # only the first pair is valid
    assert abs(wm - 0.2) < 1e-12
    assert abs(um - 0.2) < 1e-12


def test_weighted_mean_all_bad_is_nan() -> None:
    wm, wse, um = weighted_mean_se([float("nan")], [0.0])
    assert math.isnan(wm) and math.isnan(wse) and math.isnan(um)
