"""Named TUI type ramp. Ask for a scheme name from docs/TUI/STYLES.txt."""

from __future__ import annotations

GOLD = "#f6d33f"
ICE = "#73b8f4"
BLUE = "#2f7fd8"
NAVY = "#134a8d"
DIM = "#8d93a1"
DARK = "#4a5160"
INK = "#101216"
PAPER = "#e6e8eb"

ROLES: dict[str, str] = {
    "hint-key": f"bold {GOLD}",
    "hint-meaning": DIM,
    "hint-sep": DARK,
    "topbar": f"bold {GOLD}",
    "hero-tag": f"bold {GOLD}",
    "banner": ICE,
    "dim": DIM,
    "body": PAPER,
    "panel-title": f"bold {ICE}",
    "tab-active": f"bold {GOLD}",
    "tab-idle": DIM,
    "frame": BLUE,
    "live": GOLD,
    "bar": BLUE,
    "cursor": f"bold {PAPER} on {BLUE}",
    "spark": PAPER,
}


def paint(role: str, text: str) -> str:
    spec = ROLES[role]
    escaped = text.replace("[", r"\[").replace("]", r"\]")
    return f"[{spec}]{escaped}[/]"


__all__ = [
    "BLUE",
    "DARK",
    "DIM",
    "GOLD",
    "ICE",
    "INK",
    "NAVY",
    "PAPER",
    "ROLES",
    "paint",
]
