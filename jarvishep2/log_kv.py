#!/usr/bin/env python3
"""Key/value log formatting helpers for Jarvis-HEP V2."""

from __future__ import annotations

import textwrap
from typing import Any, Iterable


def _wrap_cell(text: str, width: int) -> list[str]:
    wrapped = textwrap.wrap(
        text or "",
        width=width,
        break_long_words=True,
        break_on_hyphens=False,
        drop_whitespace=False,
    )
    return wrapped or [""]


def format_two_column_log(
    title: str,
    rows: Iterable[tuple[str, Any]],
    *,
    max_width: int = 80,
    min_value_width: int = 16,
) -> str:
    """Render a compact two-column key/value block for logger output."""
    max_width = max(40, int(max_width))
    normalized = [(str(k), str(v)) for k, v in rows]
    if not normalized:
        return "\n".join(_wrap_cell(str(title), max_width))

    key_header = "Field"
    value_header = "Value"
    separator = " | "
    separator_width = len(separator)

    max_key_seen = max(len(key_header), *(len(k) for k, _ in normalized))
    max_key_width = max(8, max_width - separator_width - min_value_width)
    key_width = min(max_key_seen, max_key_width)
    value_width = max(1, max_width - separator_width - key_width)

    lines = list(_wrap_cell(str(title), max_width))
    lines.append(f"{key_header:<{key_width}}{separator}{value_header:<{value_width}}")
    lines.append(f"{'-' * key_width}{separator}{'-' * value_width}")

    for key, value in normalized:
        key_chunks = _wrap_cell(key, key_width)
        value_chunks = _wrap_cell(value, value_width)
        nline = max(len(key_chunks), len(value_chunks))
        for idx in range(nline):
            key_line = key_chunks[idx] if idx < len(key_chunks) else ""
            value_line = value_chunks[idx] if idx < len(value_chunks) else ""
            lines.append(f"{key_line:<{key_width}}{separator}{value_line:<{value_width}}")

    return "\n".join(lines)