#!/usr/bin/env python3
"""Exclusive and shared calculator PackID pools (D25.4)."""

from __future__ import annotations

import time
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.queue.redis_queue import (
    CALC_STATUS,
    MONITOR_CALC_BUSY,
    MONITOR_SAMPLE_RUNNING,
    MONITOR_WANT,
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

# Ownership HSET is unconditional. CALC_STATUS / op_count / sidecar only when
# hep:monitor:want ch:calc is on. Returns {acquired, want_calc}.
_ATOMIC_CLAIM_CALC_LUA = """
redis.call('HSET', KEYS[1], ARGV[1], ARGV[5])
local want_calc = 0
if redis.call('EXISTS', KEYS[4]) == 1 then
    local ch_calc = redis.call('HGET', KEYS[4], 'ch:calc')
    if ch_calc and ch_calc ~= '0' and ch_calc ~= '' then
        want_calc = 1
        if ARGV[4] ~= '' then
            redis.call('HSET', KEYS[5], ARGV[1], ARGV[4])
        end
        redis.call('HINCRBY', KEYS[2], ARGV[2], -1)
        redis.call('HINCRBY', KEYS[2], ARGV[3], 1)
        redis.call('INCR', KEYS[3])
    end
    local ch_sample = redis.call('HGET', KEYS[4], 'ch:sample')
    if ch_sample and ch_sample ~= '0' and ch_sample ~= '' and ARGV[4] ~= '' and ARGV[6] ~= '' then
        redis.call('HSET', KEYS[6], ARGV[4], ARGV[6])
    end
end
return {1, want_calc}
"""

_ATOMIC_RELEASE_CALC_LUA = """
local removed = redis.call('HDEL', KEYS[1], ARGV[1])
if removed == 0 then
    return {0, 0}
end
redis.call('RPUSH', KEYS[2], ARGV[1])
local want_calc = 0
if redis.call('EXISTS', KEYS[5]) == 1 then
    local ch_calc = redis.call('HGET', KEYS[5], 'ch:calc')
    if ch_calc and ch_calc ~= '0' and ch_calc ~= '' then
        want_calc = 1
        redis.call('HDEL', KEYS[6], ARGV[1])
        redis.call('HINCRBY', KEYS[3], ARGV[2], 1)
        redis.call('HINCRBY', KEYS[3], ARGV[3], -1)
        redis.call('INCR', KEYS[4])
    end
end
return {1, want_calc}
"""

_ATOMIC_RELEASE_SHARED_CALC_LUA = """
local removed = redis.call('HDEL', KEYS[1], ARGV[1])
if removed == 0 then
    return {0, 0}
end
redis.call('RPUSH', KEYS[2], ARGV[1])
if ARGV[4] == '' then
    redis.call('HDEL', KEYS[3], ARGV[1])
else
    redis.call('HSET', KEYS[3], ARGV[1], ARGV[4])
end
local want_calc = 0
if redis.call('EXISTS', KEYS[6]) == 1 then
    local ch_calc = redis.call('HGET', KEYS[6], 'ch:calc')
    if ch_calc and ch_calc ~= '0' and ch_calc ~= '' then
        want_calc = 1
        redis.call('HDEL', KEYS[7], ARGV[1])
        redis.call('HINCRBY', KEYS[4], ARGV[2], 1)
        redis.call('HINCRBY', KEYS[4], ARGV[3], -1)
        redis.call('INCR', KEYS[5])
    end
end
return {1, want_calc}
"""


def _lua_eval_unavailable(exc: BaseException) -> bool:
    return "unknown command 'eval'" in str(exc).lower()


def _unpack_lua_pair(result: Any) -> tuple[bool, bool]:
    """Unpack Lua `{flag, want_calc}`. Never treat the tuple itself as a bool."""
    if isinstance(result, (list, tuple)) and len(result) >= 2:
        return bool(int(result[0] or 0)), bool(int(result[1] or 0))
    return bool(int(result or 0)), False


def _channel_on(raw: Any) -> bool:
    text = _redis_text(raw or "").strip()
    return text not in {"", "0"}


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

    def _monitor_channels(self) -> tuple[bool, bool]:
        """(ch:calc, ch:sample). Missing want → (False, False). Fail-closed."""
        try:
            ctrl = self._ctrl()
            if int(ctrl.exists(MONITOR_WANT) or 0) != 1:
                return False, False
            calc = ctrl.hget(MONITOR_WANT, "ch:calc")
            sample = ctrl.hget(MONITOR_WANT, "ch:sample")
            return _channel_on(calc), _channel_on(sample)
        except Exception:
            return False, False

    def _claim_exclusive_pack(
        self,
        name: str,
        pack_id: str,
        *,
        worker_id: str = "",
        overlay_json: str = "",
        busy_value: str = "running",
    ) -> bool:
        """HSET ownership; return want_calc. Sidecar/status only if ch:calc."""
        want_calc = self._eval_claim_calc(
            name,
            pack_id,
            worker_id=worker_id,
            overlay_json=overlay_json,
            busy_value=busy_value,
        )
        return want_calc

    def _eval_claim_calc(
        self,
        name: str,
        pack_id: str,
        *,
        worker_id: str = "",
        overlay_json: str = "",
        busy_value: str = "running",
    ) -> bool:
        busy_key = calc_busy_packs_key(name)
        sidecar = MONITOR_CALC_BUSY.format(name=name)
        if worker_id is None or worker_id == "":
            owner = ""
        else:
            try:
                owner = str(int(worker_id))
            except (TypeError, ValueError):
                owner = str(worker_id).strip()
        overlay = str(overlay_json or "")
        try:
            result = self.r.eval(
                _ATOMIC_CLAIM_CALC_LUA,
                6,
                busy_key,
                CALC_STATUS,
                OP_COUNT.format(kind="calculator"),
                MONITOR_WANT,
                sidecar,
                MONITOR_SAMPLE_RUNNING,
                pack_id,
                calc_status_free_field(name),
                calc_status_busy_field(name),
                owner,
                busy_value,
                overlay,
            )
        except Exception as exc:
            if not _lua_eval_unavailable(exc):
                raise
            return self._claim_calc_python(
                name,
                pack_id,
                worker_id=owner,
                overlay_json=overlay,
                busy_value=busy_value,
            )
        _acquired, want_calc = _unpack_lua_pair(result)
        return want_calc

    def _claim_calc_python(
        self,
        name: str,
        pack_id: str,
        *,
        worker_id: str,
        overlay_json: str,
        busy_value: str,
    ) -> bool:
        """fakeredis-without-EVAL fallback. Tests only."""
        want_calc, want_sample = self._monitor_channels()
        pipe = self.r.pipeline(transaction=True)
        pipe.hset(calc_busy_packs_key(name), pack_id, busy_value)
        if want_calc:
            if worker_id:
                pipe.hset(MONITOR_CALC_BUSY.format(name=name), pack_id, worker_id)
            pipe.hincrby(CALC_STATUS, calc_status_free_field(name), -1)
            pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), 1)
            pipe.incr(OP_COUNT.format(kind="calculator"))
        if want_sample and worker_id and overlay_json:
            pipe.hset(MONITOR_SAMPLE_RUNNING, worker_id, overlay_json)
        pipe.execute()
        return want_calc

    def _claim_shared_pack(
        self,
        name: str,
        pack_id: str,
        current_mode: str | None,
        target_mode: str,
        *,
        worker_id: str = "",
        overlay_json: str = "",
    ) -> tuple[str, str | None, bool]:
        if not is_stable_calc_pack_id(pack_id):
            self._discard_stale_free_token(name, pack_id)
            raise ValueError(f"invalid shared calculator PackID {pack_id!r} for {name!r}")
        want_calc = self._eval_claim_calc(
            name,
            pack_id,
            worker_id=worker_id,
            overlay_json=overlay_json,
            busy_value=f"running:{target_mode}",
        )
        return pack_id, current_mode, want_calc

    def acquire_shared_calc(
        self,
        name: str,
        mode: str,
        *,
        modes: Sequence[str],
        timeout: int = 30,
        affinity_wait_sec: float = 3.0,
        worker_id: str = "",
        overlay_json: str = "",
    ) -> tuple[str, str | None] | None:
        got = self._acquire_shared_calc(
            name,
            mode,
            modes=modes,
            timeout=timeout,
            affinity_wait_sec=affinity_wait_sec,
            worker_id=worker_id,
            overlay_json=overlay_json,
        )
        if got is None:
            return None
        pack_id, current_mode, _want = got
        return pack_id, current_mode

    def _acquire_shared_calc(
        self,
        name: str,
        mode: str,
        *,
        modes: Sequence[str],
        timeout: int = 30,
        affinity_wait_sec: float = 3.0,
        worker_id: str = "",
        overlay_json: str = "",
    ) -> tuple[str, str | None, bool] | None:
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
                        parent,
                        _redis_text(raw).strip(),
                        current_mode,
                        target,
                        worker_id=worker_id,
                        overlay_json=overlay_json,
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
                            parent,
                            _redis_text(raw[1]).strip(),
                            target,
                            target,
                            worker_id=worker_id,
                            overlay_json=overlay_json,
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
                        parent,
                        _redis_text(raw).strip(),
                        current_mode,
                        target,
                        worker_id=worker_id,
                        overlay_json=overlay_json,
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
                return self._claim_shared_pack(
                    parent,
                    pack,
                    current_mode,
                    target,
                    worker_id=worker_id,
                    overlay_json=overlay_json,
                )
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
        released, _want = self._eval_atomic_release_shared_calc(
            parent, pack, free_key, target_mode,
        )
        return released

    def _eval_atomic_release_shared_calc(
        self,
        parent: str,
        pack_id: str,
        free_key: str,
        target_mode: str,
    ) -> tuple[bool, bool]:
        """Atomically return a shared PackID and record its trusted warm mode.

        Returns ``(released, want_calc)``. Packmode + affinity RPUSH are
        unconditional; sidecar/status/op_count only when ch:calc is on.
        """
        busy_key = calc_busy_packs_key(parent)
        mode_key = calc_shared_pack_mode_key(parent)
        sidecar = MONITOR_CALC_BUSY.format(name=parent)
        try:
            result = self.r.eval(
                _ATOMIC_RELEASE_SHARED_CALC_LUA,
                7,
                busy_key,
                free_key,
                mode_key,
                CALC_STATUS,
                OP_COUNT.format(kind="calculator"),
                MONITOR_WANT,
                sidecar,
                pack_id,
                calc_status_free_field(parent),
                calc_status_busy_field(parent),
                target_mode,
            )
        except Exception as exc:
            if not _lua_eval_unavailable(exc):
                raise
            if not self.r.hdel(busy_key, pack_id):
                return False, False
            want_calc, _want_sample = self._monitor_channels()
            pipe = self.r.pipeline(transaction=True)
            pipe.rpush(free_key, pack_id)
            if target_mode:
                pipe.hset(mode_key, pack_id, target_mode)
            else:
                pipe.hdel(mode_key, pack_id)
            if want_calc:
                pipe.hdel(sidecar, pack_id)
                pipe.hincrby(CALC_STATUS, calc_status_free_field(parent), 1)
                pipe.hincrby(CALC_STATUS, calc_status_busy_field(parent), -1)
                pipe.incr(OP_COUNT.format(kind="calculator"))
            pipe.execute()
            return True, want_calc
        return _unpack_lua_pair(result)

    def force_release_shared_calc(self, name: str, pack_id: str) -> bool:
        """Return an uncertain shared pack to unassigned after a Worker failure."""
        return self.release_shared_calc(name, pack_id, None)

    def _discard_stale_free_token(self, name: str, token: str) -> None:
        """Drop a non-PackID free-list entry (e.g. legacy ``ready``/UUID) permanently."""
        want_calc, _want_sample = self._monitor_channels()
        if not want_calc:
            return
        free_field = calc_status_free_field(name)
        try:
            current = int(self.r.hget(CALC_STATUS, free_field) or 0)
        except (TypeError, ValueError):
            current = 0
        if current > 0:
            self.r.hincrby(CALC_STATUS, free_field, -1)

    def acquire_calc(
        self,
        name: str,
        timeout: int = 30,
        *,
        worker_id: str = "",
        overlay_json: str = "",
    ) -> str | None:
        """Claim one free PackID exclusively (blocks until a slot or timeout).

        Only stable numeric PackIDs (``001`` …) are returned. Junk free-list
        values (``ready``, UUIDs, empty) are discarded and the wait continues.
        Returns ``None`` only when *timeout* elapses with no valid free slot.
        """
        pack, _want = self._acquire_calc(
            name, timeout=timeout, worker_id=worker_id, overlay_json=overlay_json
        )
        return pack

    def _acquire_calc(
        self,
        name: str,
        timeout: int = 30,
        *,
        worker_id: str = "",
        overlay_json: str = "",
    ) -> tuple[str | None, bool]:
        """BLPOP + claim. Returns ``(pack, want_calc)``; pack is None on timeout."""
        self._require_client()
        pool_key = calc_free_list_key(name)
        deadline = time.monotonic() + max(0.0, float(timeout))

        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return None, False
            wait_sec = max(1, int(remaining) if remaining >= 1 else 1)
            raw = self._blpop(pool_key, timeout=wait_sec)
            if raw is None:
                if time.monotonic() >= deadline:
                    return None, False
                continue

            pack_id = str(raw[1]).strip()
            if not is_stable_calc_pack_id(pack_id):
                self._discard_stale_free_token(name, pack_id)
                continue

            want_calc = self._claim_exclusive_pack(
                name,
                pack_id,
                worker_id=worker_id,
                overlay_json=overlay_json,
            )
            return pack_id, want_calc

    def release_calc(self, name: str, pack_id: str) -> None:
        """Return a PackID to free after the sample finishes (running → free)."""
        if not pack_id or not str(pack_id).strip():
            raise ValueError("pack_id is required for release_calc")
        pack_id = str(pack_id).strip()
        self._require_client()

        # Never put non-stable ids back into the free list (closes the UUID loop).
        if not is_stable_calc_pack_id(pack_id):
            busy_key = calc_busy_packs_key(name)
            self.r.hdel(busy_key, pack_id)
            return

        released, _want = self._eval_atomic_release_calc(name, pack_id)
        if not released:
            raise ValueError(f"unknown pack_id '{pack_id}' for calculator '{name}'")

    def _eval_atomic_release_calc(self, name: str, pack_id: str) -> tuple[bool, bool]:
        """Atomically transition one busy calculator slot back to free.

        Returns ``(released, want_calc)``. Must not be used as ``if not result``.
        """
        busy_key = calc_busy_packs_key(name)
        sidecar = MONITOR_CALC_BUSY.format(name=name)
        if not is_stable_calc_pack_id(pack_id):
            return bool(self.r.hdel(busy_key, pack_id)), False
        try:
            result = self.r.eval(
                _ATOMIC_RELEASE_CALC_LUA,
                6,
                busy_key,
                calc_free_list_key(name),
                CALC_STATUS,
                OP_COUNT.format(kind="calculator"),
                MONITOR_WANT,
                sidecar,
                pack_id,
                calc_status_free_field(name),
                calc_status_busy_field(name),
            )
        except Exception as exc:
            if not _lua_eval_unavailable(exc):
                raise
            if not self.r.hdel(busy_key, pack_id):
                return False, False
            want_calc, _want_sample = self._monitor_channels()
            pipe = self.r.pipeline(transaction=True)
            pipe.rpush(calc_free_list_key(name), pack_id)
            if want_calc:
                pipe.hdel(sidecar, pack_id)
                pipe.hincrby(CALC_STATUS, calc_status_free_field(name), 1)
                pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), -1)
                pipe.incr(OP_COUNT.format(kind="calculator"))
            pipe.execute()
            return True, want_calc
        return _unpack_lua_pair(result)

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
        released, _want = self._eval_atomic_release_calc(name, pack_id)
        return released

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
