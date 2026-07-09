#!/usr/bin/env python3
"""Thin bridge from Jarvis2 CLI to the Jarvis-PLOT engine.

Plot algorithms live in the JarvisPLOT package. This module only imports the
engine, injects argv, and maps failures to process-style exit codes.
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from typing import Any


class PlotBridgeError(RuntimeError):
    """Raised when the plot bridge cannot start or complete a plot run."""


def require_jarvisplot() -> Any:
    """Import ``jarvisplot`` or raise with an install hint."""
    try:
        import jarvisplot  # noqa: F401
        from jarvisplot.core import JarvisPLOT
    except ImportError as exc:
        raise ImportError(
            "Plotting requires JarvisPLOT. "
            "Install it with `pip install 'jarvishep2[plot]'` "
            "(or `pip install -e ../Jarvis-PLOT` for local development)."
        ) from exc
    return JarvisPLOT


def run_plot(
    plot_yaml: str,
    *,
    extra_args: Sequence[str] | None = None,
    plotter_cls: Any | None = None,
) -> int:
    """Run JarvisPLOT on ``plot_yaml``.

    Returns a process-style exit code (0 success, non-zero failure).
    ``plotter_cls`` is injectable for unit tests.
    """
    path = str(plot_yaml or "").strip()
    if not path:
        raise PlotBridgeError("plot YAML path must not be empty")

    if plotter_cls is None:
        plotter_cls = require_jarvisplot()

    argv = ["jplot", path, *[str(item) for item in (extra_args or ())]]
    previous_argv = sys.argv
    try:
        sys.argv = argv
        plotter = plotter_cls()
        plotter.init()
        return 0
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return int(code)
        return 1
    except Exception as exc:
        raise PlotBridgeError(f"JarvisPLOT failed on {path!r}: {exc}") from exc
    finally:
        sys.argv = previous_argv


__all__ = ["PlotBridgeError", "require_jarvisplot", "run_plot"]
