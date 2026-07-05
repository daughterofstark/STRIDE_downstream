"""Structural validation of the S6 replicate layer.

Structural / arithmetic checks only — no statistical assertions and no biological
claims (the gate is uncalibrated and residue-scale per-run outputs are
exploratory, §0.1, §5.3). Each check appends a
:class:`~stride_s6.models.ValidationCheck` to the report and raises
:class:`~stride_s6.models.errors.ConsistencyError` on failure, so later stages can
trust the layer.

Checks:

- every output table's key is unique;
- the **regime** table is arithmetically consistent — completeness counts are
  bounded (``0 ≤ n_positions_in_all_replicates ≤ n_positions``), ``frac_complete``
  matches the counts and lies in ``[0, 1]``, ``residue_claims_licensed`` matches
  ``n_replicates ≥ 5``, ``n_replicates_with_effects ≤ n_replicates``, and
  ``per_replicate_effects_available`` matches ``n_replicates_with_effects ≥ 2``;
- the **effect-spread** table is ordered and non-negative — ``n_obs ≥ 1``,
  ``theta_min ≤ theta_max``, ``theta_range`` matches ``theta_max − theta_min``, the
  range / sd / absolute-mean / pairwise-diff are ``≥ 0``, and the tier is
  ``exploratory``;
- the **concordance** table carries ``kendalls_w`` in ``[0, 1]`` (or ``nan``),
  ``mean_pairwise_spearman`` in ``[-1, 1]`` (or ``nan``), a ``concordance_class``
  in the vocabulary that matches the recomputed label, and the ``exploratory``
  tier; and
- the **blocked-analysis** ledger has a valid ``status``, an ``available`` flag
  that matches it, a non-empty reason / required-input / design reference on every
  row, and a leave-one-replicate-out row that is always blocked.
"""
from __future__ import annotations

import math

import pandas as pd

from ..build.replicate_concordance import _concordance_class
from ..models import S6Report
from ..models.errors import ConsistencyError
from ..models.schema import (
    ANALYSIS_LORO_STABILITY,
    CONCORDANCE_CLASSES,
    LEDGER_STATUSES,
    MIN_K_FOR_RESIDUE_LICENSE,
    MIN_REPLICATES_FOR_CONCORDANCE,
    ROUND_DECIMALS,
    STATUS_BLOCKED,
    TIER_EXPLORATORY,
)

_TOL = 10.0 ** (-(ROUND_DECIMALS - 1))


def _assert_unique(df: pd.DataFrame, key: list[str], what: str) -> None:
    if df.empty:
        return
    dupes = df.duplicated(subset=key, keep=False)
    if bool(dupes.any()):
        offending = df.loc[dupes, key].to_dict("records")
        raise ConsistencyError(f"{what} has non-unique key {key}: {offending}")


def validate_unique_keys(
    replicate_regime: pd.DataFrame,
    replicate_effect_spread: pd.DataFrame,
    replicate_concordance: pd.DataFrame,
    replicate_blocked_analyses: pd.DataFrame,
    report: S6Report,
) -> None:
    """Every output table's declared key is unique."""
    _assert_unique(replicate_regime, ["serotype"], "replicate_regime")
    _assert_unique(
        replicate_effect_spread,
        ["serotype", "canon_label"],
        "replicate_effect_spread",
    )
    _assert_unique(
        replicate_concordance, ["serotype"], "replicate_concordance"
    )
    _assert_unique(
        replicate_blocked_analyses,
        ["analysis_id"],
        "replicate_blocked_analyses",
    )
    report.add(
        "unique_keys",
        "all_tables",
        True,
        "all four output tables have unique keys",
    )


def validate_replicate_regime(
    replicate_regime: pd.DataFrame, report: S6Report
) -> None:
    """The regime table is arithmetically self-consistent."""
    for row in replicate_regime.to_dict("records"):
        serotype = row["serotype"]
        n_pos = int(row["n_positions"])
        n_all = int(row["n_positions_in_all_replicates"])
        k = int(row["n_replicates"])
        n_eff = int(row["n_replicates_with_effects"])
        if not (0 <= n_all <= n_pos):
            raise ConsistencyError(
                f"replicate_regime[{serotype}]: n_positions_in_all_replicates "
                f"{n_all} out of range [0, {n_pos}]"
            )
        expected_frac = (n_all / n_pos) if n_pos else 0.0
        if abs(float(row["frac_complete"]) - round(expected_frac, ROUND_DECIMALS)) > _TOL:
            raise ConsistencyError(
                f"replicate_regime[{serotype}]: frac_complete "
                f"{row['frac_complete']} != {expected_frac}"
            )
        if bool(row["residue_claims_licensed"]) != (
            k >= MIN_K_FOR_RESIDUE_LICENSE
        ):
            raise ConsistencyError(
                f"replicate_regime[{serotype}]: residue_claims_licensed "
                f"disagrees with K={k}"
            )
        if n_eff > k:
            raise ConsistencyError(
                f"replicate_regime[{serotype}]: n_replicates_with_effects "
                f"{n_eff} exceeds n_replicates {k}"
            )
        if bool(row["per_replicate_effects_available"]) != (
            n_eff >= MIN_REPLICATES_FOR_CONCORDANCE
        ):
            raise ConsistencyError(
                f"replicate_regime[{serotype}]: per_replicate_effects_available "
                f"disagrees with n_replicates_with_effects={n_eff}"
            )
    report.add(
        "regime_consistency",
        "replicate_regime",
        True,
        f"{len(replicate_regime)} regime row(s) arithmetically consistent",
    )


def validate_replicate_effect_spread(
    replicate_effect_spread: pd.DataFrame, report: S6Report
) -> None:
    """The effect-spread table is ordered, non-negative, and Tier B."""
    for row in replicate_effect_spread.to_dict("records"):
        tag = f"{row['serotype']}/{row['canon_label']}"
        if int(row["n_obs"]) < 1:
            raise ConsistencyError(
                f"replicate_effect_spread[{tag}]: n_obs < 1"
            )
        if float(row["theta_min"]) - float(row["theta_max"]) > _TOL:
            raise ConsistencyError(
                f"replicate_effect_spread[{tag}]: theta_min > theta_max"
            )
        expected_range = float(row["theta_max"]) - float(row["theta_min"])
        if abs(float(row["theta_range"]) - round(expected_range, ROUND_DECIMALS)) > _TOL:
            raise ConsistencyError(
                f"replicate_effect_spread[{tag}]: theta_range inconsistent"
            )
        for col in ("theta_range", "theta_sd", "abs_theta_mean", "max_pairwise_abs_diff"):
            if float(row[col]) < -_TOL:
                raise ConsistencyError(
                    f"replicate_effect_spread[{tag}]: {col} is negative"
                )
        if row["tier"] != TIER_EXPLORATORY:
            raise ConsistencyError(
                f"replicate_effect_spread[{tag}]: tier must be exploratory"
            )
    report.add(
        "effect_spread_consistency",
        "replicate_effect_spread",
        True,
        f"{len(replicate_effect_spread)} spread row(s) ordered and non-negative",
    )


def validate_replicate_concordance(
    replicate_concordance: pd.DataFrame, report: S6Report
) -> None:
    """The concordance table is bounded and its class matches the coefficient."""
    for row in replicate_concordance.to_dict("records"):
        serotype = row["serotype"]
        w = float(row["kendalls_w"])
        rho = float(row["mean_pairwise_spearman"])
        n_pos = int(row["n_positions_complete"])
        n_runs = int(row["n_replicates_with_effects"])
        if not math.isnan(w) and not (-_TOL <= w <= 1.0 + _TOL):
            raise ConsistencyError(
                f"replicate_concordance[{serotype}]: kendalls_w {w} outside [0, 1]"
            )
        if not math.isnan(rho) and not (-1.0 - _TOL <= rho <= 1.0 + _TOL):
            raise ConsistencyError(
                f"replicate_concordance[{serotype}]: mean_pairwise_spearman "
                f"{rho} outside [-1, 1]"
            )
        if row["concordance_class"] not in CONCORDANCE_CLASSES:
            raise ConsistencyError(
                f"replicate_concordance[{serotype}]: unknown concordance_class "
                f"{row['concordance_class']!r}"
            )
        expected = _concordance_class(w, n_pos, n_runs)
        if row["concordance_class"] != expected:
            raise ConsistencyError(
                f"replicate_concordance[{serotype}]: concordance_class "
                f"{row['concordance_class']!r} != recomputed {expected!r}"
            )
        if row["tier"] != TIER_EXPLORATORY:
            raise ConsistencyError(
                f"replicate_concordance[{serotype}]: tier must be exploratory"
            )
    report.add(
        "concordance_consistency",
        "replicate_concordance",
        True,
        f"{len(replicate_concordance)} concordance row(s) bounded and labelled",
    )


def validate_replicate_blocked_analyses(
    replicate_blocked_analyses: pd.DataFrame, report: S6Report
) -> None:
    """The ledger's statuses, flags, and mandatory fields are consistent."""
    loro_seen = False
    for row in replicate_blocked_analyses.to_dict("records"):
        aid = row["analysis_id"]
        if row["status"] not in LEDGER_STATUSES:
            raise ConsistencyError(
                f"replicate_blocked_analyses[{aid}]: unknown status "
                f"{row['status']!r}"
            )
        if bool(row["available"]) != (row["status"] != STATUS_BLOCKED):
            raise ConsistencyError(
                f"replicate_blocked_analyses[{aid}]: available flag disagrees "
                f"with status {row['status']!r}"
            )
        for col in ("description", "reason", "required_input", "design_ref"):
            if not str(row[col]).strip():
                raise ConsistencyError(
                    f"replicate_blocked_analyses[{aid}]: empty {col}"
                )
        if aid == ANALYSIS_LORO_STABILITY:
            loro_seen = True
            if row["status"] != STATUS_BLOCKED or bool(row["available"]):
                raise ConsistencyError(
                    "replicate_blocked_analyses: leave-one-replicate-out must "
                    "always be blocked"
                )
    if not loro_seen:
        raise ConsistencyError(
            "replicate_blocked_analyses: missing the leave-one-replicate-out row"
        )
    report.add(
        "blocked_ledger_consistency",
        "replicate_blocked_analyses",
        True,
        f"{len(replicate_blocked_analyses)} ledger row(s) consistent; "
        "LORO always blocked",
    )


def validate_all(
    replicate_regime: pd.DataFrame,
    replicate_effect_spread: pd.DataFrame,
    replicate_concordance: pd.DataFrame,
    replicate_blocked_analyses: pd.DataFrame,
    report: S6Report,
) -> None:
    """Run every structural check, populating ``report``."""
    validate_unique_keys(
        replicate_regime,
        replicate_effect_spread,
        replicate_concordance,
        replicate_blocked_analyses,
        report,
    )
    validate_replicate_regime(replicate_regime, report)
    validate_replicate_effect_spread(replicate_effect_spread, report)
    validate_replicate_concordance(replicate_concordance, report)
    validate_replicate_blocked_analyses(replicate_blocked_analyses, report)
