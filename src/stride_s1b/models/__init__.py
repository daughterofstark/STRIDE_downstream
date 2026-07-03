"""Report carriers for S1B (plain dataclasses serialised to annotation_summary.json)."""
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
class S1BReport:
    """Aggregate S1B report: annotation facts + every check outcome."""

    checks: list[ValidationCheck] = field(default_factory=list)
    n_residue_annotations: int = 0
    n_domain_annotations: int = 0
    n_hierarchy_annotations: int = 0
    n_serotype_annotations: int = 0
    serotypes: list[str] = field(default_factory=list)
    facts: dict[str, Any] = field(default_factory=dict)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, scope: str, passed: bool, detail: str = "") -> None:
        self.checks.append(ValidationCheck(name, scope, passed, detail))
