"""Unit tests for the deterministic rank-concordance math (``build/_concordance``)."""
from __future__ import annotations

import math

from stride_s6.build._concordance import (
    average_ranks,
    kendalls_w,
    mean_pairwise_spearman,
    spearman,
)


def test_average_ranks_no_ties() -> None:
    assert average_ranks([0.3, 0.1, 0.2]) == [3.0, 1.0, 2.0]


def test_average_ranks_with_ties() -> None:
    # two tied smallest share rank (1+2)/2 = 1.5; the third gets rank 3
    assert average_ranks([0.2, 0.9, 0.2]) == [1.5, 3.0, 1.5]


def test_average_ranks_all_tied() -> None:
    assert average_ranks([5.0, 5.0, 5.0]) == [2.0, 2.0, 2.0]


def test_kendalls_w_identical_columns_is_one() -> None:
    # three runs, five positions, identical orderings → perfect concordance
    matrix = [[float(i), float(i), float(i)] for i in range(5)]
    assert kendalls_w(matrix) == 1.0


def test_kendalls_w_two_reversed_raters_is_zero() -> None:
    # m=2, perfectly reversed ranks → row sums constant → S=0 → W=0
    matrix = [[1.0, 5.0], [2.0, 4.0], [3.0, 3.0 + 1e-9], [4.0, 2.0], [5.0, 1.0]]
    assert kendalls_w(matrix) == 0.0


def test_kendalls_w_bounds() -> None:
    matrix = [[0.1, 0.2, 0.9], [0.5, 0.4, 0.1], [0.9, 0.7, 0.3], [0.2, 0.6, 0.5]]
    w = kendalls_w(matrix)
    assert 0.0 <= w <= 1.0


def test_kendalls_w_nan_when_too_few_positions_or_runs() -> None:
    assert math.isnan(kendalls_w([[1.0, 2.0]]))  # one position
    assert math.isnan(kendalls_w([[1.0], [2.0], [3.0]]))  # one run


def test_kendalls_w_constant_runs_returns_one() -> None:
    # every run internally constant → denominator collapses → trivial agreement
    matrix = [[2.0, 2.0], [2.0, 2.0], [2.0, 2.0]]
    assert kendalls_w(matrix) == 1.0


def test_spearman_perfect_and_reversed() -> None:
    assert spearman([1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]) == 1.0
    assert spearman([1.0, 2.0, 3.0, 4.0], [40.0, 30.0, 20.0, 10.0]) == -1.0


def test_spearman_nan_on_constant_vector() -> None:
    assert math.isnan(spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]))


def test_mean_pairwise_spearman_identical_runs() -> None:
    matrix = [[float(i), float(i), float(i)] for i in range(6)]
    assert mean_pairwise_spearman(matrix) == 1.0


def test_mean_pairwise_spearman_averages_pairs() -> None:
    # run A & B identical (+1), run C reversed vs A and B (-1 each) →
    # mean of {+1, -1, -1} = -1/3
    a = [1.0, 2.0, 3.0, 4.0]
    b = [1.0, 2.0, 3.0, 4.0]
    c = [4.0, 3.0, 2.0, 1.0]
    matrix = [[a[i], b[i], c[i]] for i in range(4)]
    assert math.isclose(mean_pairwise_spearman(matrix), -1.0 / 3.0, abs_tol=1e-9)


def test_determinism_repeated_calls() -> None:
    matrix = [[0.1, 0.4, 0.2], [0.5, 0.3, 0.9], [0.9, 0.8, 0.1], [0.2, 0.6, 0.7]]
    first = kendalls_w(matrix)
    for _ in range(5):
        assert kendalls_w(matrix) == first
