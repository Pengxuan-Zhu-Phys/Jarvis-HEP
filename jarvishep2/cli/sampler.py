#!/usr/bin/env python3
"""``Jarvis sampler`` forwards to the ``jsampler`` command.

The command surface lives in Jarvis-HEP-Sampler. This module is the same
entrypoint pattern as ``Jarvis plot`` and ``jplot``.
"""

from __future__ import annotations

import sys

from jarvishep2.run_outcome import EXIT_USAGE


def dispatch_sampler(argv: list[str] | None = None) -> int:
    """Handle ``Jarvis sampler``, which is ``jsampler`` with a Jarvis program name."""
    try:
        from jarvis_hep_sampler.cli import main as sampler_main
    except ImportError as exc:
        print(
            "Jarvis-HEP-Sampler is required for `Jarvis sampler`. "
            "Install it with `pip install -U Jarvis-HEP-Sampler`. "
            f"Detail: {exc}",
            file=sys.stderr,
        )
        return EXIT_USAGE
    return int(sampler_main(list(argv or []), prog="Jarvis sampler"))


__all__ = ["dispatch_sampler"]
