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
    from jarvishep2.sampler_sdk import SDK_DISTRIBUTION, SDK_SPECIFIER, SamplerSDKError, require_sampler_sdk

    try:
        require_sampler_sdk()
    except SamplerSDKError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_USAGE
    try:
        from jarvis_hep_sampler.cli import main as sampler_main
    except ImportError as exc:
        print(
            f"{SDK_DISTRIBUTION}{SDK_SPECIFIER} is required for `Jarvis sampler`. "
            f"Install it with `pip install -U \"{SDK_DISTRIBUTION}{SDK_SPECIFIER}\"`. "
            f"Detail: {exc}",
            file=sys.stderr,
        )
        return EXIT_USAGE
    return int(sampler_main(list(argv or []), prog="Jarvis sampler"))


__all__ = ["dispatch_sampler"]
