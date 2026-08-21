#!/usr/bin/env python3
"""Exclusive and shared calculator PackID pools (D25.4)."""

from __future__ import annotations

import time
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.redis_queue import (
    CALC_STATUS,
    OP_COUNT,
    SHARED_HELD_PREFIX,
    _redis_text,
    calc_busy_packs_key,
    calc_free_list_key,
    calc_shared_free_list_key,
    calc_shared_pack_mode_key,
    calc_shared_unassigned_list_key,
    calc_status_busy_field,
    calc_status_free_field,
    format_calc_pack_id,
    is_stable_calc_pack_id,
)

_ATOMIC_RELEASE_CALC_LUA = """
local removed = redis.call('HDEL', KEYS[1], ARGV[1])
if removed == 0 then
    return 0
end
redis.call('RPUSH', KEYS[2], ARGV[1])
redis.call('HINCRBY', KEYS[3], ARGV[2], 1)
redis.call('HINCRBY', KEYS[3], ARGV[3], -1)
redis.call('INCR', KEYS[4])
return 1
"""

_ATOMIC_RELEASE_SHARED_CALC_LUA = """
local removed = redis.call('HDEL', KEYS[1], ARGV[1])
if removed == 0 then
    return 0
end
redis.call('RPUSH', KEYS[2], ARGV[1])
if ARGV[4] == '' then
    redis.call('HDEL', KEYS[3], ARGV[1])
else
    redis.call('HSET', KEYS[3], ARGV[1], ARGV[4])
end
redis.call('HINCRBY', KEYS[4], ARGV[2], 1)
redis.call('HINCRBY', KEYS[4], ARGV[3], -1)
redis.call('INCR', KEYS[5])
return 1
"""


class _CalcPool:
    """Private RedisQueue mixin (D25.4)."""

    def register_calc_pool(self, name: str, n: int) -> None:
        """Register *n* exclusive PackID slots owned by Redis.

        Lifecycle (one sample / one worker per slot at a time)::

            free list  calc:free:<name>  = [001, 002, …, N]
            acquire    BLPOP free → busy[pack]=running   (exclusive)
            release    HDEL busy  → RPUSH free           (reuse same PackID)

        Always rebuilds the pool from scratch so leftover ``ready`` / UUID
        tokens from older builds cannot re-enter the free list.
        """
        self._require_client()
        if n <= 0:
            return
        slots = int(n)
        pool_key = calc_free_list_key(name)
        busy_key = calc_busy_packs_key(name)
        pack_ids = [format_calc_pack_id(i, slots=slots) for i in range(1, slots + 1)]
        pipe = self.r.pipeline(transaction=True)
        # Hard reset: wipe free + busy so no UUID/ready residue survives.
        pipe.delete(pool_key)
        pipe.delete(busy_key)
        for pack_id in pack_ids:
            pipe.rpush(pool_key, pack_id)
        pipe.hset(
            CALC_STATUS,
            mapping={
                calc_status_free_field(name): slots,
                calc_status_busy_field(name): 0,
            },
        )
        pipe.execute()

    def register_shared_calc_pool(
        self,
        name: str,
        n: int,
        *,
        modes: list[str],
        pack_modes: Mapping[str, str] | None = None,
    ) -> None:
        """Register one physical PackID pool, partitioned by its installed mode.

        A PackID lives in exactly one mode-affinity list while free.  A failed
        rebuild goes back to ``unassigned``; a successful rebuild is returned
        to the requested mode's list.  The active busy hash remains the normal
        ``calc:busy:<parent>`` key, so a physical directory can never be held
        by two modes at once.
        """
        self._require_client()
        slots = int(n)
        normalized_modes = list(dict.fromkeys(
            str(mode).strip() for mode in modes if str(mode).strip()
        ))
        if slots <= 0 or not normalized_modes:
            return
        prefix = f"calc:free:{name}:"
        stale = [_redis_text(key) for key in self.r.scan_iter(match=prefix + "*")]
        keys = [
            calc_free_list_key(name), calc_busy_packs_key(name),
            calc_shared_pack_mode_key(name), *stale,
        ]
        pack_modes = {
            str(pack): str(mode)
            for pack, mode in dict(pack_modes or {}).items()
            if str(mode) in normalized_modes
        }
        pack_ids = [format_calc_pack_id(index, slots=slots) for index in range(1, slots + 1)]
        pipe = self.r.pipeline(transaction=True)
        pipe.delete(*list(dict.fromkeys(keys)))
        for pack_id in pack_ids:
            mode = pack_modes.get(pack_id)
            if mode:
                pipe.rpush(calc_shared_free_list_key(name, mode), pack_id)
                pipe.hset(calc_shared_pack_mode_key(name), pack_id, mode)
            else:
                pipe.rpush(calc_shared_unassigned_list_key(name), pack_id)
        pipe.hset(
            CALC_STATUS,
            mapping={
                calc_status_free_field(name): slots,
                calc_status_busy_field(name): 0,
            },
        )
        pipe.execute()

    def shared_mode_free_counts(self, name: str, modes: Sequence[str]) -> dict[str, int]:
        """Return currently free affinity-slot counts for Worker greedy ordering."""
        self._require_client()
        normalized = list(dict.fromkeys(str(mode).strip() for mode in modes if str(mode).strip()))
        if not normalized:
            return {}
        pipe = self.r.pipeline(transaction=False)
        for mode in normalized:
            pipe.llen(calc_shared_free_list_key(name, mode))
        values = pipe.execute()
        return {mode: int(value or 0) for mode, value in zip(normalized, values)}

    def shared_mode_busy_counts(self, name: str, modes: Sequence[str]) -> dict[str, int]:
        """Count in-flight shared slots by target mode for cold-start balancing."""
        self._require_client()
        normalized = list(dict.fromkeys(str(mode).strip() for mode in modes if str(mode).strip()))
        counts = {mode: 0 for mode in normalized}
        for raw_status in (self.r.hvals(calc_busy_packs_key(name)) or []):
            status = _redis_text(raw_status)
            if not status.startswith("running:"):
                continue
            mode = status.removeprefix("running:")
            if mode in counts:
                counts[mode] += 1
        return counts

    def _claim_shared_pack(
        self,
        name: str,
        pack_id: str,
        current_mode: str | None,
        target_mode: str,
    ) -> tuple[str, str | None]:
        if not is_stable_calc_pack_id(pack_id):
            self._discard_stale_free_token(name, pack_id)
            raise ValueError(f"invalid shared calculator PackID {pack_id!r} for {name!r}")
        pipe = self.r.pipeline(transaction=True)
        # This label has no ownership semantics (the hash key and PackID do);
        # it lets concurrently starting Workers spread cold packs across modes.
        pipe.hset(calc_busy_packs_key(name), pack_id, f"running:{target_mode}")
        pipe.hincrby(CALC_STATUS, calc_status_free_field(name), -1)
        pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), 1)
        pipe.incr(OP_COUNT.format(kind="calculator"))
        pipe.execute()
        return pack_id, current_mode

    def acquire_shared_calc(
        self,
        name: str,
        mode: str,
        *,
        modes: Sequence[str],
        timeout: int = 30,
        affinity_wait_sec: float = 3.0,
    ) -> tuple[str, str | None] | None:
        """Acquire one parent PackID, preferring a pack already built for *mode*.

        This is the broker half of affinity scheduling.  Workers choose a
        preferred pending mode from :meth:`shared_mode_free_counts`; this method
        then takes an exact match first, an unassigned pack second, and borrows
        the most plentiful other mode only when necessary. If an exact warm
        pack is already running, wait briefly for it before rebuilding another
        mode's pack; this preserves affinity under saturated Worker contention.
        """
        self._require_client()
        parent = str(name or "").strip()
        target = str(mode or "").strip()
        all_modes = list(dict.fromkeys(str(item).strip() for item in modes if str(item).strip()))
        if not parent or not target or target not in all_modes:
            raise ValueError("shared calculator acquire requires a declared parent and mode")
        deadline = time.monotonic() + max(0.0, float(timeout))
        target_key = calc_shared_free_list_key(parent, target)
        unassigned_key = calc_shared_unassigned_list_key(parent)
        warm_wait_budget = max(0.0, float(affinity_wait_sec))
        while True:
            # Exact and never-built packs are always preferable and cheap.
            for key, current_mode in (
                (target_key, target),
                (unassigned_key, None),
            ):
                raw = self.r.lpop(key)
                if raw is None:
                    continue
                try:
                    return self._claim_shared_pack(
                        parent, _redis_text(raw).strip(), current_mode, target,
                    )
                except ValueError:
                    continue

            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return None

            # A peer Worker is already preparing/running this target mode. Let
            # its warm pack return instead of immediately destroying another
            # useful affinity. Do not wait when no target pack is in flight:
            # that would penalize pools smaller than the number of modes.
            target_busy = self.shared_mode_busy_counts(parent, [target]).get(target, 0)
            if target_busy > 0 and warm_wait_budget > 0:
                warm_wait = min(remaining, warm_wait_budget)
                raw = self._blpop_many([target_key], timeout=warm_wait)
                if raw is not None:
                    try:
                        return self._claim_shared_pack(
                            parent, _redis_text(raw[1]).strip(), target, target,
                        )
                    except ValueError:
                        continue
                if deadline - time.monotonic() <= 0:
                    return None

            other_modes = sorted(
                (item for item in all_modes if item != target),
                key=lambda item: (-int(self.r.llen(calc_shared_free_list_key(parent, item)) or 0), item),
            )
            candidates: list[tuple[str, str | None]] = [
                (target_key, target),
                (unassigned_key, None),
                *((calc_shared_free_list_key(parent, item), item) for item in other_modes),
            ]
            for key, current_mode in candidates:
                raw = self.r.lpop(key)
                if raw is None:
                    continue
                try:
                    return self._claim_shared_pack(
                        parent, _redis_text(raw).strip(), current_mode, target,
                    )
                except ValueError:
                    continue
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return None
            raw = self._blpop_many(
                [key for key, _ in candidates], timeout=remaining,
            )
            if raw is None:
                continue
            popped_key, pack = _redis_text(raw[0]), _redis_text(raw[1]).strip()
            current_mode = next((candidate_mode for key, candidate_mode in candidates if key == popped_key), None)
            try:
                return self._claim_shared_pack(parent, pack, current_mode, target)
            except ValueError:
                continue

    def release_shared_calc(self, name: str, pack_id: str, mode: str | None) -> bool:
        """Release a shared PackID to its successful mode list or unassigned."""
        self._require_client()
        parent = str(name or "").strip()
        pack = str(pack_id or "").strip()
        if not parent or not is_stable_calc_pack_id(pack):
            return False
        target_mode = str(mode or "").strip()
        free_key = (
            calc_shared_free_list_key(parent, target_mode)
            if target_mode else calc_shared_unassigned_list_key(parent)
        )
        return self._eval_atomic_release_shared_calc(
            parent, pack, free_key, target_mode,
        )

    def _eval_atomic_release_shared_calc(
        self,
        parent: str,
        pack_id: str,
        free_key: str,
        target_mode: str,
    ) -> bool:
        """Atomically return a shared PackID and record its trusted warm mode."""
        busy_key = calc_busy_packs_key(parent)
        mode_key = calc_shared_pack_mode_key(parent)
        try:
            result = self.r.eval(
                _ATOMIC_RELEASE_SHARED_CALC_LUA,
                5,
                busy_key,
                free_key,
                mode_key,
                CALC_STATUS,
                OP_COUNT.format(kind="calculator"),
                pack_id,
                calc_status_free_field(parent),
                calc_status_busy_field(parent),
                target_mode,
            )
        except Exception as exc:
            # Keep the same fakeredis compatibility contract as normal pools.
            if "unknown command 'eval'" not in str(exc).lower():
                raise
            if not self.r.hdel(busy_key, pack_id):
                return False
            pipe = self.r.pipeline(transaction=True)
            pipe.rpush(free_key, pack_id)
            if target_mode:
                pipe.hset(mode_key, pack_id, target_mode)
            else:
                pipe.hdel(mode_key, pack_id)
            pipe.hincrby(CALC_STATUS, calc_status_free_field(parent), 1)
            pipe.hincrby(CALC_STATUS, calc_status_busy_field(parent), -1)
            pipe.incr(OP_COUNT.format(kind="calculator"))
            pipe.execute()
            return True
        return bool(int(result or 0))

    def force_release_shared_calc(self, name: str, pack_id: str) -> bool:
        """Return an uncertain shared pack to unassigned after a Worker failure."""
        return self.release_shared_calc(name, pack_id, None)

    def _discard_stale_free_token(self, name: str, token: str) -> None:
        """Drop a non-PackID free-list entry (e.g. legacy ``ready``/UUID) permanently."""
        free_field = calc_status_free_field(name)
        try:
            current = int(self.r.hget(CALC_STATUS, free_field) or 0)
        except (TypeError, ValueError):
            current = 0
        if current > 0:
            self.r.hincrby(CALC_STATUS, free_field, -1)

    def acquire_calc(self, name: str, timeout: int = 30) -> str | None:
        """Claim one free PackID exclusively (blocks until a slot or timeout).

        Only stable numeric PackIDs (``001`` …) are returned. Junk free-list
        values (``ready``, UUIDs, empty) are discarded and the wait continues.
        Returns ``None`` only when *timeout* elapses with no valid free slot.
        """
        self._require_client()
        pool_key = calc_free_list_key(name)
        busy_key = calc_busy_packs_key(name)
        deadline = time.monotonic() + max(0.0, float(timeout))

        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return None
            # redis-py BLPOP timeout is whole seconds; keep waiting in chunks.
            wait_sec = max(1, int(remaining) if remaining >= 1 else 1)
            raw = self._blpop(pool_key, timeout=wait_sec)
            if raw is None:
                if time.monotonic() >= deadline:
                    return None
                continue

            pack_id = str(raw[1]).strip()
            if not is_stable_calc_pack_id(pack_id):
                # Legacy free-list pollution — drop and wait for a real slot.
                self._discard_stale_free_token(name, pack_id)
                continue

            # free → running (exclusive ownership for this worker/sample).
            pipe = self.r.pipeline(transaction=True)
            pipe.hset(busy_key, pack_id, "running")
            pipe.hincrby(CALC_STATUS, calc_status_free_field(name), -1)
            pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), 1)
            pipe.incr(OP_COUNT.format(kind="calculator"))
            pipe.execute()
            return pack_id

    def release_calc(self, name: str, pack_id: str) -> None:
        """Return a PackID to free after the sample finishes (running → free)."""
        if not pack_id or not str(pack_id).strip():
            raise ValueError("pack_id is required for release_calc")
        pack_id = str(pack_id).strip()
        self._require_client()

        # Never put non-stable ids back into the free list (closes the UUID loop).
        if not is_stable_calc_pack_id(pack_id):
            busy_key = calc_busy_packs_key(name)
            # Best-effort cleanup if an old worker still held a UUID slot.
            self.r.hdel(busy_key, pack_id)
            return

        busy_key = calc_busy_packs_key(name)
        result = self._eval_atomic_release_calc(name, pack_id)
        if not result:
            raise ValueError(f"unknown pack_id '{pack_id}' for calculator '{name}'")

    def _eval_atomic_release_calc(self, name: str, pack_id: str) -> bool:
        """Atomically transition one busy calculator slot back to free."""
        busy_key = calc_busy_packs_key(name)
        if not is_stable_calc_pack_id(pack_id):
            return bool(self.r.hdel(busy_key, pack_id))
        try:
            result = self.r.eval(
                _ATOMIC_RELEASE_CALC_LUA,
                4,
                busy_key,
                calc_free_list_key(name),
                CALC_STATUS,
                OP_COUNT.format(kind="calculator"),
                pack_id,
                calc_status_free_field(name),
                calc_status_busy_field(name),
            )
        except Exception as exc:
            # fakeredis versions used by the offline test suite may omit EVAL;
            # real Redis always takes the atomic Lua path above.
            if "unknown command 'eval'" not in str(exc).lower():
                raise
            if not self.r.hdel(busy_key, pack_id):
                return False
            pipe = self.r.pipeline(transaction=True)
            pipe.rpush(calc_free_list_key(name), pack_id)
            pipe.hincrby(CALC_STATUS, calc_status_free_field(name), 1)
            pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), -1)
            pipe.incr(OP_COUNT.format(kind="calculator"))
            pipe.execute()
            return True
        return bool(int(result or 0))

    def force_release_calc(self, name: str, pack_id: str) -> bool:
        """Best-effort PackID return for failure/cleanup paths (never raises).

        Returns True when this call transitioned the slot from busy → free.
        Safe under double-release races (second caller gets False).
        """
        if not pack_id or not str(pack_id).strip():
            return False
        pack_id = str(pack_id).strip()
        name = str(name or "").strip()
        if not name:
            return False
        self._require_client()
        # False means the slot was already free/not owned; Redis/EVAL errors
        # propagate so callers can distinguish them and retain local ownership.
        return self._eval_atomic_release_calc(name, pack_id)

    def sweep_held_calc_slots(self, held_packs: Mapping[str, Any]) -> int:
        """Release calculator slots recorded for a dead Worker (WP-D6.1)."""
        released = 0
        for name, pack_id in dict(held_packs or {}).items():
            if not pack_id or not str(pack_id).strip():
                continue
            try:
                label = str(name)
                if label.startswith(SHARED_HELD_PREFIX):
                    self.force_release_shared_calc(
                        label.removeprefix(SHARED_HELD_PREFIX), str(pack_id)
                    )
                else:
                    self.release_calc(label, str(pack_id))
            except ValueError:
                continue
            released += 1
        return released
