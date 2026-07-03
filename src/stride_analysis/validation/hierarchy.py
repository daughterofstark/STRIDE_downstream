"""Canonical hierarchy-path parsing (Level-2 grammar).

A STRIDE ``region_id`` encodes the hierarchy as a "/"-joined path::

    complex / protein / chain / domain / motif / secondary_structure / residue

A region at scale *s* is a prefix of the residue path truncated to that scale's
depth. This module validates depths and turns a residue-scale path into explicit
``h_<level>`` columns so downstream code never string-splits.
"""
from __future__ import annotations

from ..models.errors import HierarchyError
from ..models.schema import (
    HIERARCHY_COLUMNS,
    HIERARCHY_LEVELS,
    SCALE_LEVEL_PATH_DEPTH,
)


def split_region_id(region_id: str) -> list[str]:
    """Split a region_id, rejecting empty input or empty segments."""
    if region_id is None or region_id == "":
        raise HierarchyError("empty region_id")
    segments = region_id.split("/")
    for i, seg in enumerate(segments):
        if seg == "":
            raise HierarchyError(
                f"malformed hierarchy path {region_id!r}: empty segment at "
                f"position {i}"
            )
    return segments


def validate_path_depth(region_id: str, scale_level: str) -> None:
    """Assert ``region_id`` has the depth required for ``scale_level``."""
    expected = SCALE_LEVEL_PATH_DEPTH.get(scale_level)
    if expected is None:
        raise HierarchyError(
            f"unknown scale_level {scale_level!r} "
            f"(known: {sorted(SCALE_LEVEL_PATH_DEPTH)})"
        )
    depth = len(split_region_id(region_id))
    if depth != expected:
        raise HierarchyError(
            f"region_id {region_id!r} at scale {scale_level!r} has depth "
            f"{depth}, expected {expected}"
        )


def parse_residue_path(residue_region_id: str) -> dict[str, str]:
    """Parse a residue-scale region_id into ``{h_<level>: segment}``."""
    segments = split_region_id(residue_region_id)
    if len(segments) != len(HIERARCHY_LEVELS):
        raise HierarchyError(
            f"residue path {residue_region_id!r} has {len(segments)} segments, "
            f"expected {len(HIERARCHY_LEVELS)} ({'/'.join(HIERARCHY_LEVELS)})"
        )
    return dict(zip(HIERARCHY_COLUMNS, segments, strict=True))
