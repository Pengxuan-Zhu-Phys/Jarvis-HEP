#!/usr/bin/env python3
"""Redis calculator pool registration, including D20 shared mode packs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.calculator_modes import mode_info, shared_mode_groups
from jarvishep2.redis_queue import RedisQueue, format_calc_pack_id


def _modules(worker_config: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    raw = worker_config.get("calculator_modules") or []
    if isinstance(raw, Mapping):
        raw = list(raw.values())
    return [item for item in raw if isinstance(item, Mapping)] if isinstance(raw, Sequence) else []


def resolve_calculator_pools(worker_config: Mapping[str, Any] | None) -> dict[str, int]:
    """Return physical-pool sizes (one parent pool for all shared modes)."""
    if not isinstance(worker_config, Mapping):
        return {}
    explicit = worker_config.get("calculator_pools")
    explicit_pools = {
        str(name): max(1, int(count))
        for name, count in explicit.items()
        if str(name).strip()
    } if isinstance(explicit, Mapping) and explicit else {}
    global_parallel = worker_config.get("calculator_make_paraller")
    pools: dict[str, int] = {}
    for item in _modules(worker_config):
        info = mode_info(item)
        name = info[0] if info is not None else str(item.get("name", "")).strip()
        if not name:
            continue
        slots = explicit_pools.get(name, item.get("make_paraller", global_parallel or 1))
        pools[name] = max(1, int(slots or 1))
    # Preserve legacy explicit-only behavior for a worker blueprint with no
    # modules, while rejecting impossible mode keys at task-card validation.
    return pools or explicit_pools


def _shared_pack_modes(
    modules: list[Mapping[str, Any]],
    *,
    parent: str,
    slots: int,
    modes: list[str],
) -> dict[str, str]:
    """Rehydrate affinity from successful on-disk shared-pack stamps."""
    # Do not import from ``jarvishep2.Module`` at module import time: Module's
    # package initializer reaches the Worker runtime checkpoint, which imports
    # this pool module.  Registration runs after that import graph is complete.
    from jarvishep2.Module.runtime_preparer import read_install_stamp

    representative = next(
        (item for item in modules if (mode_info(item) or (None, None))[0] == parent),
        None,
    )
    if representative is None:
        return {}
    basepath = str(representative.get("path") or "")
    if "@PackID" not in basepath:
        return {}
    restored: dict[str, str] = {}
    for index in range(1, slots + 1):
        pack_id = format_calc_pack_id(index, slots=slots)
        runtime = basepath.replace("@PackID", pack_id)
        stamp = read_install_stamp(runtime)
        extra = stamp.get("extra") if isinstance(stamp, Mapping) else None
        mode = str(extra.get("mode") or "") if isinstance(extra, Mapping) else ""
        if mode in modes:
            restored[pack_id] = mode
    return restored


def register_calculator_pools(redis: RedisQueue, worker_config: Mapping[str, Any] | None) -> None:
    """Seed normal pools and shared-mode affinity pools before Workers start."""
    if not isinstance(worker_config, Mapping):
        return
    modules = _modules(worker_config)
    groups = shared_mode_groups(modules)
    sizes = resolve_calculator_pools(worker_config)
    for name, slots in sizes.items():
        if name not in groups:
            redis.register_calc_pool(name, slots)
            continue
        mode_names = groups[name]
        redis.register_shared_calc_pool(
            name,
            slots,
            modes=mode_names,
            pack_modes=_shared_pack_modes(
                modules, parent=name, slots=slots, modes=mode_names,
            ),
        )


__all__ = ["register_calculator_pools", "resolve_calculator_pools"]
