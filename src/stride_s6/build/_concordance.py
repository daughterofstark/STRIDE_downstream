r"""Deterministic rank-concordance statistics for the per-run effects.

Pure functions over plain Python lists — no SciPy, no randomness, fully
deterministic — implementing the design's §3.1 per-run rank concordance:
**Kendall's coefficient of concordance** *W* across the runs and the **mean
pairwise Spearman** correlation of the per-run effect rankings.

Kendall's *W* for *m* runs (raters) ranking *N* positions, with average-rank tie
handling and the standard tie correction, is

.. math::

    W = \frac{12 S}{m^2 (N^3 - N) - m \sum_j T_j},

where :math:`S = \sum_i (R_i - \bar R)^2`, :math:`R_i` is the sum over runs of the
average ranks of position *i*, and :math:`T_j = \sum (t^3 - t)` sums over the tie
groups of run *j* (each of size *t*). *W* lies in ``[0, 1]``; it is undefined
(returned as ``nan``) for ``N < 2`` or ``m < 2``, and is ``1.0`` when the
denominator collapses because every run is internally constant.
"""
from __future__ import annotations

import math


def average_ranks(values: list[float]) -> list[float]:
    """Average (fractional) ranks of ``values``, ties sharing the mean rank.

    Ranks are 1-based and assigned in ascending order; deterministic under a
    stable sort. For ``[0.2, 0.9, 0.2]`` returns ``[1.5, 3.0, 1.5]``.
    """
    n = len(values)
    order = sorted(range(n), key=lambda i: (values[i], i))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0  # mean of 1-based positions i+1..j+1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def _tie_correction(ranks: list[float]) -> float:
    """``Σ (t³ − t)`` over the tie groups implied by shared ranks."""
    counts: dict[float, int] = {}
    for r in ranks:
        counts[r] = counts.get(r, 0) + 1
    return float(sum(t**3 - t for t in counts.values() if t > 1))


def kendalls_w(matrix: list[list[float]]) -> float:
    """Kendall's *W* for a position × run θ matrix (rows = positions, cols = runs).

    Returns ``nan`` when there are fewer than two positions or fewer than two runs.
    """
    n = len(matrix)
    if n < 2:
        return math.nan
    m = len(matrix[0])
    if m < 2:
        return math.nan

    # rank positions within each run (column)
    col_ranks: list[list[float]] = []
    tie_sum = 0.0
    for j in range(m):
        column = [matrix[i][j] for i in range(n)]
        r = average_ranks(column)
        col_ranks.append(r)
        tie_sum += _tie_correction(r)

    # R_i = sum over runs of position i's rank
    row_sums = [sum(col_ranks[j][i] for j in range(m)) for i in range(n)]
    mean_rs = sum(row_sums) / n
    s = sum((rs - mean_rs) ** 2 for rs in row_sums)

    denom = (m**2) * (n**3 - n) - m * tie_sum
    if denom <= 0.0:
        # every run internally constant → perfect (trivial) agreement
        return 1.0
    return 12.0 * s / denom


def spearman(x: list[float], y: list[float]) -> float:
    """Spearman rank correlation of two equal-length vectors.

    Pearson correlation of the average-rank vectors. Returns ``nan`` when a vector
    has zero rank variance (all tied) or fewer than two elements.
    """
    n = len(x)
    if n < 2 or len(y) != n:
        return math.nan
    rx = average_ranks(x)
    ry = average_ranks(y)
    mx = sum(rx) / n
    my = sum(ry) / n
    cov = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    vx = sum((rx[i] - mx) ** 2 for i in range(n))
    vy = sum((ry[i] - my) ** 2 for i in range(n))
    if vx <= 0.0 or vy <= 0.0:
        return math.nan
    return cov / math.sqrt(vx * vy)


def mean_pairwise_spearman(matrix: list[list[float]]) -> float:
    """Mean Spearman correlation over all run pairs of a position × run matrix.

    Returns ``nan`` when there are fewer than two positions, fewer than two runs,
    or when no run pair yields a defined correlation.
    """
    n = len(matrix)
    if n < 2:
        return math.nan
    m = len(matrix[0])
    if m < 2:
        return math.nan
    columns = [[matrix[i][j] for i in range(n)] for j in range(m)]
    vals: list[float] = []
    for a in range(m):
        for b in range(a + 1, m):
            rho = spearman(columns[a], columns[b])
            if not math.isnan(rho):
                vals.append(rho)
    if not vals:
        return math.nan
    return sum(vals) / len(vals)
