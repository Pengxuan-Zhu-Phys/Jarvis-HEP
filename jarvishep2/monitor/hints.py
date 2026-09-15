"""Frozen bottom key-hint chrome.

Locked visual contract (do not restyle):

    Enter: attach  |  J/K: select  |  Q: quit

* Key name: bold, gold ``#f6d33f`` (``Enter``, ``J/K``, ``Q``, …).
* ``: meaning``: dim ``#8d93a1``.
* Groups joined by space-bar-space, bar in ``#4a5160``.

Every Monitor page footer goes through ``render_key_hint``. Bindings may
change; the format does not.
"""

from __future__ import annotations

from jarvishep2.monitor.styles import DARK, DIM, GOLD, paint

KEY_COLOR = GOLD
MEANING_COLOR = DIM
SEP_COLOR = DARK

SPLASH_KEYS: tuple[tuple[str, str], ...] = (
    ("Enter", "attach"),
    ("J/K", "select"),
    ("R", "refresh"),
    ("Q", "quit"),
)

SESSION_KEYS: tuple[tuple[str, str], ...] = (
    ("Esc", "chooser"),
    ("Q", "quit"),
)

OVERVIEW_KEYS: tuple[tuple[str, str], ...] = (
    ("Esc", "chooser"),
    ("Tab", "next"),
    ("←/→", "page"),
    ("1-7", "pages"),
    ("Q", "quit"),
)

OVERVIEW_ONLY_KEYS: tuple[tuple[str, str], ...] = (
    ("Esc", "chooser"),
    ("R", "refresh"),
    ("Q", "quit"),
)


def render_key_hint(pairs: tuple[tuple[str, str], ...] | list[tuple[str, str]]) -> str:
    """Bold gold key names, dim ``: meaning``, groups split by `` | ``."""
    groups: list[str] = []
    for key, meaning in pairs:
        groups.append(paint("hint-key", key) + paint("hint-meaning", f": {meaning}"))
    return f" {paint('hint-sep', '|')} ".join(groups)


def render_key_hint_plain(
    pairs: tuple[tuple[str, str], ...] | list[tuple[str, str]],
) -> str:
    """Colourless stand-in for the txt canvases."""
    return "  |  ".join(f"{key}: {meaning}" for key, meaning in pairs)


__all__ = [
    "OVERVIEW_KEYS",
    "OVERVIEW_ONLY_KEYS",
    "SESSION_KEYS",
    "SPLASH_KEYS",
    "render_key_hint",
    "render_key_hint_plain",
]
