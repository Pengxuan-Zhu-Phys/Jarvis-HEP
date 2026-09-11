"""Splash address bar: git branch, cwd path, clock."""

from __future__ import annotations

import subprocess
from datetime import datetime
from pathlib import Path


def git_branch(root: Path | None = None) -> str:
    cwd = Path.cwd() if root is None else root
    try:
        result = subprocess.run(
            ["git", "-C", str(cwd), "branch", "--show-current"],
            check=False,
            capture_output=True,
            text=True,
            timeout=0.8,
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    if result.returncode != 0:
        return ""
    return result.stdout.strip()


def home_relative_path(path: Path | None = None) -> str:
    resolved = (Path.cwd() if path is None else path).expanduser().resolve()
    home = Path.home().resolve()
    try:
        relative = resolved.relative_to(home)
    except ValueError:
        return str(resolved)
    return "~" if not relative.parts else f"~/{relative.as_posix()}"


def clock_label(now: datetime | None = None) -> str:
    dt = datetime.now() if now is None else now
    hour = dt.strftime("%I").lstrip("0") or "12"
    return f"{hour}:{dt.strftime('%M %p')}"


def compact_path(path: str, max_chars: int) -> str:
    if max_chars <= 0:
        return ""
    if len(path) <= max_chars:
        return path
    if max_chars <= 3:
        return path[:max_chars]
    return "..." + path[-(max_chars - 3) :]


def render_topbar_left(*, branch: str | None = None, path: str | None = None) -> str:
    """Path, plus ``⎇ branch ·`` only when a git branch exists.

    The clock owns the second ``·``. Do not put a trailing cdot here.
    """
    label = git_branch() if branch is None else branch
    shown = home_relative_path() if path is None else path
    if label:
        return f"⎇ {label} · {shown}"
    return shown


def render_topbar_clock(now: datetime | None = None) -> str:
    """Second cdot lives immediately before the clock."""
    return f"· {clock_label(now)}"


__all__ = [
    "clock_label",
    "compact_path",
    "git_branch",
    "home_relative_path",
    "render_topbar_clock",
    "render_topbar_left",
]
