#!/usr/bin/env python3
"""Verify every canvas in docs/TUI/<name>-WxH/ is exactly W columns by H rows.

Run from repo root:  python3 docs/TUI/_check.py
Does not rewrite files. Canvases are the source of truth.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def expected_size(folder: str) -> tuple[int, int] | None:
    if "x" not in folder:
        return None
    tail = folder.rsplit("-", 1)[-1]
    if "x" not in tail:
        return None
    a, b = tail.split("x", 1)
    if not (a.isdigit() and b.isdigit()):
        return None
    return int(a), int(b)


def main() -> int:
    bad: list[str] = []
    checked = 0
    for path in sorted(ROOT.rglob("*.txt")):
        size = expected_size(path.parent.name)
        if size is None:
            continue
        width, height = size
        lines = path.read_text(encoding="utf-8").splitlines()
        rel = path.relative_to(ROOT)
        if len(lines) != height:
            bad.append(f"{rel}: {len(lines)} rows, expected {height}")
        for i, line in enumerate(lines, 1):
            if len(line) != width:
                bad.append(f"{rel}:{i}: {len(line)} cols, expected {width}")
                break
        checked += 1
    if bad:
        print("\n".join(bad))
        print(f"FAIL  {len(bad)} problem(s), {checked} file(s)")
        return 1
    print(f"OK  {checked} canvas(es) match their folder WxH")
    return 0


if __name__ == "__main__":
    sys.exit(main())
