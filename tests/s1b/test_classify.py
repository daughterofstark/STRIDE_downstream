"""Tests for the deterministic classification helpers."""
from __future__ import annotations

import pytest

from stride_s1b.build._classify import (
    availability_class,
    conservation_class,
    domain_status,
    is_resolved,
    secondary_structure_status,
)
from stride_s1b.models.schema import (
    AVAIL_ALL,
    AVAIL_NONE,
    AVAIL_SOME,
    CONSERVATION_PAN,
    CONSERVATION_PARTIAL,
    CONSERVATION_UNIQUE,
    DOMAIN_STATUS_ASSIGNED,
    DOMAIN_STATUS_UNASSIGNED,
    SS_STATUS_RESOLVED,
    SS_STATUS_UNRESOLVED,
)


@pytest.mark.parametrize(
    "value,expected",
    [
        ("Catalytic Triad", True),
        ("helix", True),
        ("unassigned", False),
        ("none", False),
        ("unknown", False),
        ("", False),
        ("  UNKNOWN  ", False),  # case/space-insensitive
        (None, False),
    ],
)
def test_is_resolved(value: object, expected: bool) -> None:
    assert is_resolved(value) is expected  # type: ignore[arg-type]


def test_domain_status() -> None:
    assert domain_status("Catalytic Triad") == DOMAIN_STATUS_ASSIGNED
    assert domain_status("unassigned") == DOMAIN_STATUS_UNASSIGNED


def test_secondary_structure_status() -> None:
    assert secondary_structure_status("helix") == SS_STATUS_RESOLVED
    assert secondary_structure_status("unknown") == SS_STATUS_UNRESOLVED


@pytest.mark.parametrize(
    "present,total,expected",
    [
        (3, 3, CONSERVATION_PAN),
        (2, 3, CONSERVATION_PARTIAL),
        (1, 3, CONSERVATION_UNIQUE),
        (1, 1, CONSERVATION_PAN),   # single-serotype dataset: trivially pan
        (0, 0, CONSERVATION_UNIQUE),  # degenerate guard
    ],
)
def test_conservation_class(present: int, total: int, expected: str) -> None:
    assert conservation_class(present, total) == expected


@pytest.mark.parametrize(
    "n,available,in_all,expected",
    [
        (3, True, True, AVAIL_ALL),
        (1, True, False, AVAIL_SOME),
        (0, False, False, AVAIL_NONE),
    ],
)
def test_availability_class(
    n: int, available: bool, in_all: bool, expected: str
) -> None:
    assert availability_class(n, available, in_all) == expected
