#!/usr/bin/env python3
"""Redis calculator free-pool registration helpers (WP-D2.1)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.redis_queue import RedisQueue


def resolve_calculator_pools(worker_config: Mapping[str, Any] | None) -> dict[str, int]:
    """Derive per-calculator slot counts from explicit or module config."""
    if not isinstance(worker_config, Mapping):
        return {}

    explicit = worker_config.get("calculator_pools")
    if isinstance(explicit, Mapping) and explicit:
        return {
            str(name): max(1, int(count))
            for name, count in explicit.items()
            if str(name).strip()
        }

    global_parallel = worker_config.get("calculator_make_paraller")
    modules = worker_config.get("calculator_modules") or []
    if isinstance(modules, dict):
        modules = list(modules.values())

    pools: dict[str, int] = {}
    for item in modules:
        if not isinstance(item, Mapping):
            continue
        name = str(item.get("name", "")).strip()
        if not name:
            continue
        slots = item.get("make_paraller", global_parallel or 1)
        pools[name] = max(1, int(slots or 1))
    return pools


def register_calculator_pools(redis: RedisQueue, worker_config: Mapping[str, Any] | None) -> None:
    """Seed ``calc:free:<name>`` lists and ``hep:calculator:status`` counters."""
    for name, slots in resolve_calculator_pools(worker_config).items():
        redis.register_calc_pool(name, slots)


__all__ = ["register_calculator_pools", "resolve_calculator_pools"]