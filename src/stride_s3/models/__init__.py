"""Report carriers for S3 (plain dataclasses serialised to hierarchy_summary.json)."""
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
class S3Report:
    """Aggregate S3 report: hierarchy-reduction facts, provenance, and checks.

    ``provenance`` carries the design's required per-run header (§5.4): the input
    file digest, the provisional ρ*, ``calibrated=false``, and the K note.
    """

    checks: list[ValidationCheck] = field(default_factory=list)
    n_scale_curve: int = 0
    n_resolution_gap: int = 0
    n_monotonicity_audit: int = 0
    n_chain_contrast: int = 0
    serotypes: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    facts: dict[str, Any] = field(default_factory=dict)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, scope: str, passed: bool, detail: str = "") -> None:
        self.checks.append(ValidationCheck(name, scope, passed, detail))
