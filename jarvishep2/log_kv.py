#!/usr/bin/env python3
"""Key/value log formatting helpers for Jarvis-HEP V2."""

from __future__ import annotations

import textwrap
import time
from typing import Any, Iterable


def format_duration(seconds: float) -> str:
    """V1-compatible ``HH:MM:SS.msc`` duration for progress heartbeats."""
    total = max(0.0, float(seconds))
    hours = int(total // 3600)
    minutes = int((total % 3600) // 60)
    secs = total % 60
    millisec = secs % 1
    # Match V1: fractional part of seconds as three-digit milliseconds.
    frac = str(millisec)
    ms = frac[2:5] if "." in frac else "000"
    ms = (ms + "000")[:3]
    return f"{hours:02d}:{minutes:02d}:{int(secs):02d}.{ms}"


class PermilleProgress:
    """Emit V1-style ‰ progress lines; WARNING at exact 1% milestones."""

    def __init__(
        self,
        logger: Any,
        *,
        total: int,
        label: str = "samples finished",
        t0: float | None = None,
    ) -> None:
        self._logger = logger
        self.total = max(0, int(total))
        self.label = str(label)
        self.t0 = float(time.time() if t0 is None else t0)
        self._permille = -1

    def update(
        self,
        done: int,
        *,
        extra: str | None = None,
        force: bool = False,
    ) -> bool:
        """Log when the ‰ changes (or ``force``). Return True if a line was emitted."""
        if self.total <= 0:
            return False
        done_n = max(0, int(done))
        permille = min(1000, int(done_n / float(self.total) * 1000))
        if not force and permille == self._permille:
            return False
        self._permille = permille
        msg = "{}‰ of {}/{} {} in {}".format(
            permille,
            done_n,
            self.total,
            self.label,
            format_duration(time.time() - self.t0),
        )
        if extra:
            msg = f"{msg} ({extra})"
        if permille > 0 and permille % 10 == 0:
            self._logger.warning(msg)
        else:
            self._logger.info(msg)
        return True


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