"""Structural validation of the S3 hierarchy-reduction layer.

Structural / arithmetic checks only — no statistical assertions and no
biological claims (the gate is uncalibrated, §0.1). Each check appends a
:class:`~stride_s3.models.ValidationCheck` to the report and raises
:class:`~stride_s3.models.errors.ConsistencyError` on failure, so later stages
can trust the reduction layer.

Checks:

- every output table's key is unique;
- the scale curve is *complete* — every locus carries exactly one row per scale
  (0..6), and the residue/domain ρ echoed by the gap table match the curve;
- the resolution gap is *arithmetically consistent* — ``delta`` equals
  ``rho_domain - rho_residue`` on every row;
- the monotonicity audit is *self-consistent* — ``is_monotone`` iff
  ``n_violations == 0``, and a monotone row has zero ``max_decrease``;
- the chain contrast *partitions* each chain's mechanisms
  (increase + decrease + mixed = n_mechanisms; increase + decrease = n_signed).
"""
from __future__ import annotations

import math

import pandas as pd

from ..models import S3Report
from ..models.errors import ConsistencyError
from ..models.schema import N_SCALES, RHO_DECIMALS

_TOL = 10.0 ** (-(RHO_DECIMALS - 1))


def validate_unique_keys(
    scale_curve: pd.DataFrame,
    resolution_gap: pd.DataFrame,
    monotonicity_audit: pd.DataFrame,
    chain_contrast: pd.DataFrame,
    report: S3Report,
) -> None:
    """Every output table's declared key is unique."""
    _assert_unique(
        scale_curve,
        ["serotype", "canon_label", "scale_index"],
        "scale_curve",
    )
    _assert_unique(
        resolution_gap, ["serotype", "canon_label"], "resolution_gap"
    )
    _assert_unique(
        monotonicity_audit, ["serotype", "canon_label"], "monotonicity_audit"
    )
    _assert_unique(chain_contrast, ["serotype", "chain"], "chain_contrast")
    report.add(
        "every output table key is unique",
        "global",
        True,
        "4 tables checked",
    )


def validate_scale_curve_completeness(
    scale_curve: pd.DataFrame,
    resolution_gap: pd.DataFrame,
    report: S3Report,
) -> None:
    """Every locus has exactly N_SCALES rows; gap ρ echoes the curve."""
    if scale_curve.empty:
        report.add("scale curve complete", "global", True, "empty")
        return

    counts = scale_curve.groupby(["serotype", "canon_label"]).size()
    bad = counts[counts != N_SCALES]
    if not bad.empty:
        example = bad.index[0]
        raise ConsistencyError(
            f"locus {example} has {int(bad.iloc[0])} scale rows; "
            f"expected exactly {N_SCALES}"
        )
    # every locus must carry the full 0..6 index set exactly once
    idx_sets = scale_curve.groupby(["serotype", "canon_label"])[
        "scale_index"
    ].apply(lambda s: tuple(sorted(int(x) for x in s)))
    expected = tuple(range(N_SCALES))
    for key, got in idx_sets.items():
        if got != expected:
            raise ConsistencyError(
                f"locus {key} scale indices {got} != expected {expected}"
            )

    # gap table's residue/domain ρ must match the curve
    if not resolution_gap.empty:
        curve_res = _rho_at(scale_curve, 0)
        curve_dom = _rho_at(scale_curve, 3)
        for row in resolution_gap.itertuples(index=False):
            key = (row.serotype, row.canon_label)
            _assert_close(
                row.rho_residue,
                curve_res.get(key, float("nan")),
                f"resolution_gap rho_residue for {key}",
            )
            _assert_close(
                row.rho_domain,
                curve_dom.get(key, float("nan")),
                f"resolution_gap rho_domain for {key}",
            )
    report.add(
        "scale curve complete (one row per scale) and gap ρ matches curve",
        "global",
        True,
        f"{len(counts)} loci × {N_SCALES} scales",
    )


def validate_gap_consistency(
    resolution_gap: pd.DataFrame, report: S3Report
) -> None:
    """Δρ equals ``rho_domain - rho_residue`` on every gap row."""
    if resolution_gap.empty:
        report.add("resolution gap arithmetic consistent", "global", True, "empty")
        return
    for row in resolution_gap.itertuples(index=False):
        expected = float(row.rho_domain) - float(row.rho_residue)
        _assert_close(
            row.delta_rho_domain_residue,
            expected,
            f"resolution_gap delta for ({row.serotype}, {row.canon_label})",
        )
    report.add(
        "resolution gap arithmetic consistent (Δρ = domain − residue)",
        "global",
        True,
        f"{len(resolution_gap)} loci checked",
    )


def validate_monotonicity_audit_consistency(
    monotonicity_audit: pd.DataFrame, report: S3Report
) -> None:
    """``is_monotone`` iff ``n_violations == 0``; monotone → zero max_decrease."""
    if monotonicity_audit.empty:
        report.add("monotonicity audit self-consistent", "global", True, "empty")
        return
    for row in monotonicity_audit.itertuples(index=False):
        if bool(row.is_monotone) != (int(row.n_violations) == 0):
            raise ConsistencyError(
                f"monotonicity audit for ({row.serotype}, {row.canon_label}): "
                f"is_monotone={row.is_monotone} disagrees with "
                f"n_violations={row.n_violations}"
            )
        if bool(row.is_monotone):
            if float(row.max_decrease) != 0.0:
                raise ConsistencyError(
                    f"monotone locus ({row.serotype}, {row.canon_label}) has "
                    f"non-zero max_decrease={row.max_decrease}"
                )
            if int(row.first_violation_scale_index) != -1:
                raise ConsistencyError(
                    f"monotone locus ({row.serotype}, {row.canon_label}) has a "
                    f"violation index {row.first_violation_scale_index}"
                )
        else:
            if float(row.max_decrease) <= 0.0:
                raise ConsistencyError(
                    f"non-monotone locus ({row.serotype}, {row.canon_label}) "
                    f"has non-positive max_decrease={row.max_decrease}"
                )
    n_non_monotone = int((~monotonicity_audit["is_monotone"]).sum())
    report.add(
        "monotonicity audit self-consistent (flag matches violation count)",
        "global",
        True,
        f"{len(monotonicity_audit)} loci, {n_non_monotone} non-monotone",
    )


def validate_chain_contrast_totals(
    chain_contrast: pd.DataFrame, report: S3Report
) -> None:
    """Direction counts partition each chain's mechanisms."""
    if chain_contrast.empty:
        report.add("chain contrast totals partition", "global", True, "empty")
        return
    for row in chain_contrast.itertuples(index=False):
        if (
            int(row.n_increase) + int(row.n_decrease) + int(row.n_mixed)
            != int(row.n_mechanisms)
        ):
            raise ConsistencyError(
                f"chain ({row.serotype}, {row.chain}): increase + decrease + "
                f"mixed != n_mechanisms ({row.n_increase} + {row.n_decrease} + "
                f"{row.n_mixed} != {row.n_mechanisms})"
            )
        if int(row.n_increase) + int(row.n_decrease) != int(row.n_signed):
            raise ConsistencyError(
                f"chain ({row.serotype}, {row.chain}): increase + decrease != "
                f"n_signed ({row.n_increase} + {row.n_decrease} != "
                f"{row.n_signed})"
            )
    report.add(
        "chain contrast direction counts partition each chain's mechanisms",
        "global",
        True,
        f"{len(chain_contrast)} (serotype, chain) rows checked",
    )


def _rho_at(scale_curve: pd.DataFrame, scale_index: int) -> dict[tuple, float]:
    sub = scale_curve[scale_curve["scale_index"] == scale_index]
    return {
        (str(s), str(c)): float(r)
        for s, c, r in zip(
            sub["serotype"], sub["canon_label"], sub["rho"], strict=True
        )
    }


def _assert_close(got: float, expected: float, what: str) -> None:
    g = float(got)
    e = float(expected)
    if math.isnan(g) and math.isnan(e):
        return
    if math.isnan(g) or math.isnan(e) or abs(g - e) > _TOL:
        raise ConsistencyError(f"{what}: {got} != {expected}")


def _assert_unique(df: pd.DataFrame, key: list[str], what: str) -> None:
    if df.empty:
        return
    dup = df.duplicated(key, keep=False)
    if dup.any():
        examples = df.loc[dup, key].drop_duplicates().head(3).to_dict("records")
        raise ConsistencyError(
            f"{what} key {tuple(key)} not unique; examples: {examples}"
        )
