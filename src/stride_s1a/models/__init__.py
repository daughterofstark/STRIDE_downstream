"""Report carriers for S1A (plain dataclasses serialised to dataset_summary.json)."""
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
class S1AReport:
    """Aggregate S1A report: table facts + every check outcome."""

    checks: list[ValidationCheck] = field(default_factory=list)
    n_canonical_residues: int = 0
    n_domains: int = 0
    n_conserved_all: int = 0
    n_union: int = 0
    serotypes: list[str] = field(default_factory=list)
    facts: dict[str, Any] = field(default_factory=dict)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, scope: str, passed: bool, detail: str = "") -> None:
        self.checks.append(ValidationCheck(name, scope, passed, detail))
