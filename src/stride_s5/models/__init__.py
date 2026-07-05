"""Report carriers for S5 (plain dataclasses serialised to conservation_summary.json)."""
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
class S5Report:
    """Aggregate S5 report: cross-serotype facts, provenance, and checks.

    ``provenance`` carries the design's required per-run header (§5.4): the input
    file digests, the provisional ρ*, ``calibrated=false``, and the n = 4
    replication note.
    """

    checks: list[ValidationCheck] = field(default_factory=list)
    n_position_conservation: int = 0
    n_direction_concordance: int = 0
    n_domain_serotype_matrix: int = 0
    n_cross_serotype_scorecard: int = 0
    serotypes: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    facts: dict[str, Any] = field(default_factory=dict)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, scope: str, passed: bool, detail: str = "") -> None:
        self.checks.append(ValidationCheck(name, scope, passed, detail))
