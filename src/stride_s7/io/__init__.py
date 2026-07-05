"""S7 IO subpackage: loaders for the S2–S6 inputs, writers for S7 artifacts."""
from __future__ import annotations

from .loaders import (
    INPUT_REQUIRED_COLUMNS,
    STAGE_INPUTS,
    file_digest,
    load_inputs,
    resolve_input_paths,
)
from .writers import (
    write_dataframe,
    write_manifest,
    write_markdown_table,
    write_summary,
    write_svg,
)

__all__ = [
    "STAGE_INPUTS",
    "INPUT_REQUIRED_COLUMNS",
    "resolve_input_paths",
    "load_inputs",
    "file_digest",
    "write_dataframe",
    "write_svg",
    "write_markdown_table",
    "write_manifest",
    "write_summary",
]
