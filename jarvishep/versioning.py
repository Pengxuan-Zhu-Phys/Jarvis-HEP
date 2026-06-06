#!/usr/bin/env python3
from __future__ import annotations

import os
import re


RESET = "\033[0m"
BOLD = "\033[1m"
ICON_PIXEL = "⬤ "
ICON_COLORS = {
    "B": "#2f7fd8",
    "W": "#ffffff",
    "Y": "#f6d33f",
    ".": "#134a8d",
}
TEXT_GRADIENT_START = ICON_COLORS["B"]
TEXT_GRADIENT_END = ICON_COLORS["Y"]
TEXT_GRADIENT_COLORS = [
    "#2f7fd8",
    "#73b8f4",
    "#aee4fc",
    "#b8ffe4",
    "#d5f884",
    "#f8f675",
    "#ffe95d",
    "#ffdc3f",
]
ICON_TEMPLATE_RE = re.compile(r"^([BWY.]{8})( {2})(.*)$")


def _hex_to_rgb(color: str) -> tuple[int, int, int]:
    color = color.removeprefix("#")
    return int(color[0:2], 16), int(color[2:4], 16), int(color[4:6], 16)


def _ansi_fg(color: str) -> str:
    r, g, b = _hex_to_rgb(color)
    return f"\033[38;2;{r};{g};{b}m"


def _render_icon_cells(cells: str) -> str:
    return "".join(f"{_ansi_fg(ICON_COLORS[cell])}{ICON_PIXEL}{RESET}" for cell in cells)


def _render_logo_template_line(line: str, gradient_index: int, gradient_total: int) -> str:
    matched = ICON_TEMPLATE_RE.match(line)
    if not matched:
        return line

    cells, separator, text = matched.groups()
    if gradient_index < len(TEXT_GRADIENT_COLORS):
        text_color = TEXT_GRADIENT_COLORS[gradient_index]
    elif gradient_total <= 1:
        text_color = TEXT_GRADIENT_START
    else:
        text_color = TEXT_GRADIENT_END
    text_style = BOLD if "Just a Robust" in text or "Author:" in text or "Version:" in text else ""
    return f"{_render_icon_cells(cells)}{separator}{_ansi_fg(text_color)}{text_style}{text}{RESET}"


def _render_logo_template_lines(lines: list[str]) -> list[str]:
    gradient_total = sum(1 for line in lines if ICON_TEMPLATE_RE.match(line))
    gradient_index = 0
    rendered: list[str] = []
    for line in lines:
        if ICON_TEMPLATE_RE.match(line):
            rendered.append(_render_logo_template_line(line, gradient_index, gradient_total))
            gradient_index += 1
        else:
            rendered.append(line)
    return rendered


def get_runtime_version() -> str:
    """Best-effort runtime package version for banner display."""
    try:
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("Jarvis-HEP")
        except PackageNotFoundError:
            return version("jarvis-hep")
    except Exception:
        pass

    pyproject = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "pyproject.toml"))
    if os.path.exists(pyproject):
        try:
            with open(pyproject, "r", encoding="utf-8") as f:
                text = f.read()
            matched = re.search(r'(?m)^version\s*=\s*"([^"]+)"$', text)
            if matched:
                return matched.group(1)
        except Exception:
            pass

    return "unknown"


def render_logo_with_version(logo_path: str) -> str:
    """Render logo text and inject runtime version line."""
    version = get_runtime_version()
    try:
        with open(logo_path, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
    except Exception:
        return (
            "=== Jarvis-HEP ===\n"
            "Author: Pengxuan Zhu, Erdong Guo.\n"
            f"Version: {version}"
        )

    # Keep banner compact: trim trailing blank lines from template.
    while lines and not lines[-1].strip():
        lines.pop()

    replaced = False
    for idx, line in enumerate(lines):
        if "Version:" in line:
            prefix = line.split("Version:", 1)[0]
            lines[idx] = f"{prefix}Version:  {version}"
            replaced = True
            break

    if not replaced:
        prefix = "          "
        for line in lines:
            if "Author:" in line:
                prefix = line.split("Author:", 1)[0]
                break
        lines.append(f"{prefix}Version:  {version}")

    lines = _render_logo_template_lines(lines)

    return "\n".join(lines)
