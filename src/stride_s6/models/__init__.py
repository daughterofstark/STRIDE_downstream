"""Report carriers for S6 (plain dataclasses serialised to replicate_summary.json)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ValidationCheck:
    name: str
    scope: str
    passed: bool
    detail: str = ""


@dataclass
class S6Report:
    """Aggregate S6 report: replicate-layer facts, provenance, and checks.

    ``provenance`` carries the design's required per-run header (§5.4): the input
    file digests, ``calibrated=false``, the K = 3 replication note, and the
    detected per-run-effect availability. ``blocked_analyses`` mirrors the ledger
    table so the summary JSON records exactly which replicate-level analyses are
    unavailable and why.
    """

    checks: list[ValidationCheck] = field(default_factory=list)
    n_replicate_regime: int = 0
    n_replicate_effect_spread: int = 0
    n_replicate_concordance: int = 0
    n_replicate_blocked_analyses: int = 0
    serotypes: list[str] = field(default_factory=list)
    per_replicate_effects_available: bool = False
    provenance: dict[str, Any] = field(default_factory=dict)
    facts: dict[str, Any] = field(default_factory=dict)
    blocked_analyses: list[dict[str, Any]] = field(default_factory=list)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, scope: str, passed: bool, detail: str = "") -> None:
        self.checks.append(ValidationCheck(name, scope, passed, detail))
