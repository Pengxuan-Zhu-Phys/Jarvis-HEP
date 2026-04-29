#!/usr/bin/env python3
from __future__ import annotations

import os
import re


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

    return "\n".join(lines)
