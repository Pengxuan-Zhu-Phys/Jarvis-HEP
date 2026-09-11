"""Jarvis-HEP logo + ``Jarvis -v`` banner for the monitor splash.

Layout matches Jarvis-Agent's home hero: an 8×8 ⬤ monitor on the left and
the wordmark / tagline / authors on the right. Colours come from
``jarvishep2.versioning``.
"""

from __future__ import annotations

from dataclasses import dataclass

from jarvishep2.versioning import (
    ICON_TEMPLATE_RE,
    TEXT_GRADIENT_COLORS,
    default_logo_path,
    get_runtime_version,
)

FALLBACK_PATTERN = (
    "BBBW....",
    "BBWWY...",
    "BWWWY...",
    "BBWWYY..",
    "BBBWY...",
    "WWWWYYYY",
    "BBWWYY..",
    "BBBWY...",
)

LOGO_MONITOR_ROWS = 8
LOGO_MONITOR_COLS = 8
COLS_PER_SIDE = 4
REVEAL_FRAMES = 8
MONITOR_FRAMES = 18

LEFT_MONITOR_SERIES = (
    (1, 2, 1, 3, 2, 4, 1, 2),
    (2, 1, 3, 2, 4, 2, 3, 1),
    (1, 3, 4, 1, 2, 3, 2, 4),
    (3, 2, 1, 4, 3, 1, 4, 2),
    (2, 4, 2, 3, 1, 4, 2, 1),
    (4, 1, 2, 2, 3, 2, 1, 3),
)
RIGHT_MONITOR_SERIES = (
    (2, 1, 3, 1, 4, 2, 1, 3),
    (1, 3, 2, 4, 2, 1, 3, 2),
    (3, 2, 4, 1, 3, 2, 4, 1),
    (4, 1, 2, 3, 1, 4, 2, 3),
    (2, 4, 1, 2, 3, 1, 4, 2),
    (1, 2, 3, 4, 2, 3, 1, 4),
)

LEFT_BACKGROUND = "#2f7fd8"
LEFT_ACTIVE = "#ffffff"
RIGHT_BACKGROUND = "#134a8d"
RIGHT_ACTIVE = "#f6d33f"
INACTIVE = "#303745"
TAG_COLOR = "#f6d33f"


@dataclass(frozen=True)
class MonitorBranding:
    logo_pattern: tuple[str, ...]
    banner_lines: tuple[str, ...]
    version: str


def load_branding() -> MonitorBranding:
    version = get_runtime_version()
    pattern, banners = _parse_logo_file(default_logo_path(), version)
    return MonitorBranding(
        logo_pattern=pattern or FALLBACK_PATTERN,
        banner_lines=banners,
        version=version,
    )


def logo_widths(pattern: tuple[str, ...], side: str) -> tuple[int, ...]:
    if side == "left":
        return tuple(row[:4].count("W") for row in pattern)
    if side == "right":
        return tuple(row[4:8].count("Y") for row in pattern)
    raise ValueError(f"unknown logo side {side!r}")


def render_logo_monitor_frame(
    frame: int,
    pattern: tuple[str, ...],
    *,
    animate: bool = True,
) -> str:
    """Textual markup for the 8×8 ⬤ grid (Agent ``#logo-monitor``)."""
    left = logo_widths(pattern, "left")
    right = logo_widths(pattern, "right")
    if not animate:
        frame = REVEAL_FRAMES + MONITOR_FRAMES
    reveal_row = min(frame, LOGO_MONITOR_ROWS - 1)
    final_phase = frame >= REVEAL_FRAMES + MONITOR_FRAMES
    rows: list[str] = []
    for y in range(LOGO_MONITOR_ROWS):
        cells: list[str] = []
        for x in range(LOGO_MONITOR_COLS):
            color = _cell_color(
                x,
                y,
                frame,
                reveal_row,
                final_phase,
                pattern,
                left,
                right,
            )
            cells.append(f"[{color}]⬤[/]")
        rows.append(" ".join(cells))
    return "\n".join(rows)


def render_banner_markup(lines: tuple[str, ...]) -> str:
    """Gradient wordmark matching ``Jarvis -v`` / Agent home-panel."""
    rendered: list[str] = []
    for index, line in enumerate(lines):
        escaped = _escape_markup(line)
        if index < len(TEXT_GRADIENT_COLORS):
            bold = (
                " bold"
                if "Just a Robust" in line or "Author:" in line or "Version:" in line
                else ""
            )
            rendered.append(
                f"[{TEXT_GRADIENT_COLORS[index]}{bold}]{escaped}[/]"
            )
        else:
            rendered.append(escaped)
    return "\n".join(rendered)


def monitor_tag_markup() -> str:
    return (
        f"[{TAG_COLOR} bold]Jarvis Monitor[/]\n"
        "read-only runtime dashboard"
    )


def _parse_logo_file(
    path: str,
    version: str,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    try:
        text = open(path, encoding="utf-8").read()
    except OSError:
        return FALLBACK_PATTERN, (
            "JARVIS",
            "Just a Robust and Versatile Interface Suite for HEP",
            f"Author: Pengxuan Zhu, Erdong Guo.  Version:  {version}",
        )
    pattern: list[str] = []
    banners: list[str] = []
    for line in text.splitlines():
        matched = ICON_TEMPLATE_RE.match(line)
        if not matched:
            continue
        pattern.append(matched.group(1))
        rest = matched.group(3)
        if "Version:" in rest:
            rest = f"{rest.split('Version:', 1)[0]}Version:  {version}"
        if "Jarvis-HEP" in rest and "V2" not in rest:
            rest = rest.replace("Jarvis-HEP", "Jarvis-HEP V2", 1)
        banners.append(rest)
    return tuple(pattern[:8]), tuple(banners[:8])


def _cell_color(
    x: int,
    y: int,
    frame: int,
    reveal_row: int,
    final_phase: bool,
    pattern: tuple[str, ...],
    left_widths: tuple[int, ...],
    right_widths: tuple[int, ...],
) -> str:
    if y > reveal_row:
        return INACTIVE
    row = pattern[y] if y < len(pattern) else "." * LOGO_MONITOR_COLS
    cell = row[x] if x < len(row) else "."
    if final_phase:
        if x < COLS_PER_SIDE:
            return LEFT_ACTIVE if cell == "W" else LEFT_BACKGROUND
        return RIGHT_ACTIVE if cell == "Y" else RIGHT_BACKGROUND
    if x < COLS_PER_SIDE:
        widths = _current_widths("left", frame, left_widths, right_widths)
        active = x >= COLS_PER_SIDE - widths[y]
        return LEFT_ACTIVE if active else LEFT_BACKGROUND
    widths = _current_widths("right", frame, left_widths, right_widths)
    active = x - COLS_PER_SIDE < widths[y]
    return RIGHT_ACTIVE if active else RIGHT_BACKGROUND


def _current_widths(
    side: str,
    frame: int,
    left_widths: tuple[int, ...],
    right_widths: tuple[int, ...],
) -> tuple[int, ...]:
    series = LEFT_MONITOR_SERIES if side == "left" else RIGHT_MONITOR_SERIES
    index = max(0, frame - REVEAL_FRAMES)
    if index >= MONITOR_FRAMES - 6:
        return left_widths if side == "left" else right_widths
    return series[index % len(series)]


def _escape_markup(text: str) -> str:
    return text.replace("[", r"\[").replace("]", r"\]")


__all__ = [
    "MonitorBranding",
    "load_branding",
    "logo_widths",
    "monitor_tag_markup",
    "render_banner_markup",
    "render_logo_monitor_frame",
    "REVEAL_FRAMES",
    "MONITOR_FRAMES",
]
