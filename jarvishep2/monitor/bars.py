"""Locked progress bar: C rail, fill gradients navy → ice left to right."""

from __future__ import annotations

import re

from jarvishep2.monitor.styles import DARK, ICE, NAVY, paint

BLOCKS = " ▏▎▍▌▋▊▉█"
_MARKUP = re.compile(r"\[/?[^\]]*\]")
_NAVY = (0x13, 0x4A, 0x8D)
_ICE = (0x73, 0xB8, 0xF4)


def eighths(value: float, width: int, empty: str = "─") -> str:
    """Fill ``width`` cells with 1/8-block precision."""
    if width <= 0:
        return ""
    value = max(0.0, min(1.0, value))
    x = value * width
    full = min(width, int(x))
    if full >= width:
        return "█" * width
    frac = int(round((x - full) * 8))
    if frac >= 8:
        full += 1
        frac = 0
        if full >= width:
            return "█" * width
    partial = BLOCKS[frac] if frac > 0 else ""
    used = full + (1 if partial else 0)
    return "█" * full + partial + (empty * (width - used))


def _hex(rgb: tuple[int, int, int]) -> str:
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def _lerp(t: float) -> str:
    t = max(0.0, min(1.0, t))
    rgb = tuple(int(_NAVY[i] + (_ICE[i] - _NAVY[i]) * t) for i in range(3))
    return _hex(rgb)  # type: ignore[return-value]


def bar_rail(value: float, width: int) -> str:
    """C  ├████▌────────┤  fill colour: navy at 0% → ice at 100% of the track."""
    inner = max(1, width - 2)
    cells = eighths(value, inner, empty="─")
    parts = [paint("frame", "├")]
    last = ""
    run: list[str] = []

    def flush() -> None:
        nonlocal last, run
        if not run:
            return
        parts.append(f"[{last}]{''.join(run)}[/]")
        run = []

    for index, char in enumerate(cells):
        if char == "─":
            color = DARK
        else:
            color = _lerp(index / max(1, inner - 1))
        if color != last and run:
            flush()
        last = color
        run.append(char)
    flush()
    parts.append(paint("frame", "┤"))
    return "".join(parts)


def visible_width(markup: str) -> int:
    return len(_MARKUP.sub("", markup))


# Palette check: STYLES navy/ice stay the source of the gradient endpoints.
assert NAVY.lower() == "#134a8d"
assert ICE.lower() == "#73b8f4"


__all__ = ["bar_rail", "eighths", "visible_width"]
