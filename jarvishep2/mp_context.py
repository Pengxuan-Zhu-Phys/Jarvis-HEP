#!/usr/bin/env python3
"""Spawn-only multiprocessing context (invariant #10)."""

from __future__ import annotations

import multiprocessing as mp

_SPAWN_CTX: mp.context.BaseContext | None = None


def get_spawn_context() -> mp.context.BaseContext:
    """Return the process-local spawn context."""
    global _SPAWN_CTX
    if _SPAWN_CTX is None:
        _SPAWN_CTX = mp.get_context("spawn")
    return _SPAWN_CTX


__all__ = ["get_spawn_context"]