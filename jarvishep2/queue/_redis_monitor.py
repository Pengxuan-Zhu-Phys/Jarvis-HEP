#!/usr/bin/env python3
"""Opt-in inflight telemetry + occupancy reads (monitor inflight design).

Want is last-writer-wins HASH+EXPIRE on the control client. Never SET NX.
Never claims hep:control:lock. SCAN of hep:monitor:* is reset/resume/stop only.
"""

from __future__ import annotations

import json
import os
import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.queue.redis_queue import (
    MONITOR_CALC_BUSY,
    MONITOR_KEY_PATTERN,
    MONITOR_SAMPLE_RUNNING,
    MONITOR_WANT,
    MONITOR_WANT_TTL_SEC,
    _redis_text,
    calc_busy_packs_key,
    calc_free_list_key,
    is_stable_calc_pack_id,
)


_ATOMIC_SET_MONITOR_WANT_LUA = """
redis.call('HSET', KEYS[1],
    'ch:calc', ARGV[1],
    'ch:sample', ARGV[2],
    'ts', ARGV[3],
    'pid', ARGV[4])
redis.call('EXPIRE', KEYS[1], ARGV[5])
return 1
"""


def _lua_eval_unavailable(exc: BaseException) -> bool:
    return "unknown command 'eval'" in str(exc).lower()


class _MonitorTelemetry:
    """Private RedisQueue mixin: opt-in inflight telemetry + occupancy reads."""

    def monitor_calc_busy_key(self, name: str) -> str:
        return MONITOR_CALC_BUSY.format(name=str(name or "").strip())

    def set_monitor_want(
        self,
        *,
        ch_calc: bool,
        ch_sample: bool,
        ttl_sec: int = MONITOR_WANT_TTL_SEC,
        pid: int | None = None,
    ) -> None:
        """Atomic HSET+EXPIRE on r_ctrl. Last-writer-wins. Never SET NX."""
        self._require_client()
        ttl = max(1, int(ttl_sec))
        owner = str(os.getpid() if pid is None else int(pid))
        ts = str(time.time())
        calc = "1" if ch_calc else "0"
        sample = "1" if ch_sample else "0"
        ctrl = self._ctrl()
        try:
            ctrl.eval(
                _ATOMIC_SET_MONITOR_WANT_LUA,
                1,
                MONITOR_WANT,
                calc,
                sample,
                ts,
                owner,
                ttl,
            )
        except Exception as exc:
            if not _lua_eval_unavailable(exc):
                raise
            pipe = ctrl.pipeline(transaction=True)
            pipe.hset(
                MONITOR_WANT,
                mapping={
                    "ch:calc": calc,
                    "ch:sample": sample,
                    "ts": ts,
                    "pid": owner,
                },
            )
            pipe.expire(MONITOR_WANT, ttl)
            pipe.execute()

    def clear_monitor_want(self) -> None:
        """DEL hep:monitor:want on r_ctrl. Idempotent."""
        self._require_client()
        self._ctrl().delete(MONITOR_WANT)

    def drop_monitor_keys(self) -> int:
        """SCAN match hep:monitor:* then DEL. Reset / resume / Core stopping only."""
        self._require_client()
        ctrl = self._ctrl()
        keys = [_redis_text(key) for key in ctrl.scan_iter(match=MONITOR_KEY_PATTERN)]
        if not keys:
            return 0
        return int(ctrl.delete(*keys) or 0)

    def fetch_calc_occupancy(
        self,
        names: list[str],
        *,
        slots: Mapping[str, int] | None = None,
    ) -> dict[str, dict[str, int]]:
        """Return {name: {free, busy, slots}}. No SCAN.

        busy = HLEN(calc:busy:{name}) always.
        Exclusive pools: free = LLEN(calc:free:{name}); slots defaults to
        free+busy when *slots* omits the name.
        Shared pools: LLEN(calc:free:{name}) is 0 by construction. free =
        slots-busy when slots[name] is provided; otherwise free/slots are 0
        (busy stays honest).
        """
        self._require_client()
        slot_map = {
            str(name): int(value)
            for name, value in dict(slots or {}).items()
            if str(name).strip()
        }
        cleaned = [str(name).strip() for name in names if str(name).strip()]
        if not cleaned:
            return {}
        ctrl = self._ctrl()
        pipe = ctrl.pipeline(transaction=False)
        for name in cleaned:
            pipe.hlen(calc_busy_packs_key(name))
            pipe.llen(calc_free_list_key(name))
        rows = pipe.execute()
        out: dict[str, dict[str, int]] = {}
        for index, name in enumerate(cleaned):
            busy = int(rows[2 * index] or 0)
            listed_free = int(rows[2 * index + 1] or 0)
            if listed_free > 0:
                free = listed_free
                slot_n = slot_map.get(name, free + busy)
            elif name in slot_map:
                slot_n = slot_map[name]
                free = max(0, slot_n - busy)
            else:
                free = 0
                slot_n = 0
            out[name] = {"free": int(free), "busy": int(busy), "slots": int(slot_n)}
        return out

    def fetch_calc_busy_pack_ids(
        self,
        names: list[str],
    ) -> dict[str, tuple[int, ...]]:
        """Return occupied physical PackIDs for known pools in one pipeline.

        The caller supplies the bounded pool catalogue from runtime metadata;
        this never discovers keys and never reads calculator payloads.
        """
        self._require_client()
        cleaned = [str(name).strip() for name in names if str(name).strip()]
        if not cleaned:
            return {}
        pipe = self._ctrl().pipeline(transaction=False)
        for name in cleaned:
            pipe.hkeys(calc_busy_packs_key(name))
        rows = pipe.execute()
        return {
            name: tuple(
                sorted(
                    int(text)
                    for raw in (keys or ())
                    for text in [_redis_text(raw).strip()]
                    if is_stable_calc_pack_id(text)
                )
            )
            for name, keys in zip(cleaned, rows)
        }

    def publish_sample_overlay(
        self,
        worker_id: str,
        payload: Mapping[str, Any],
    ) -> None:
        """HSET hep:monitor:sample:running if ch:sample. No-op when want is off."""
        self._require_client()
        try:
            owner = str(int(worker_id))
        except (TypeError, ValueError):
            owner = str(worker_id or "").strip()
        if not owner:
            return
        ctrl = self._ctrl()
        if int(ctrl.exists(MONITOR_WANT) or 0) != 1:
            return
        if _redis_text(ctrl.hget(MONITOR_WANT, "ch:sample") or "").strip() in {"", "0"}:
            return
        ctrl.hset(
            MONITOR_SAMPLE_RUNNING,
            owner,
            json.dumps(dict(payload), separators=(",", ":")),
        )

    def clear_sample_overlay(self, worker_id: str) -> None:
        """HDEL overlay field if ch:sample. No-op when want is off."""
        self._require_client()
        try:
            owner = str(int(worker_id))
        except (TypeError, ValueError):
            owner = str(worker_id or "").strip()
        if not owner:
            return
        ctrl = self._ctrl()
        if int(ctrl.exists(MONITOR_WANT) or 0) != 1:
            return
        if _redis_text(ctrl.hget(MONITOR_WANT, "ch:sample") or "").strip() in {"", "0"}:
            return
        ctrl.hdel(MONITOR_SAMPLE_RUNNING, owner)

    def fetch_monitor_calc_busy(self, name: str) -> dict[str, str]:
        """HGETALL hep:monitor:calc:busy:{name}. Missing → {}."""
        self._require_client()
        text = str(name or "").strip()
        if not text:
            return {}
        raw = self._ctrl().hgetall(self.monitor_calc_busy_key(text)) or {}
        return {_redis_text(key): _redis_text(value) for key, value in dict(raw).items()}

    def fetch_monitor_sample_running(self) -> dict[str, dict[str, Any]]:
        """HGETALL hep:monitor:sample:running, JSON-decode values."""
        self._require_client()
        raw = self._ctrl().hgetall(MONITOR_SAMPLE_RUNNING) or {}
        decoded: dict[str, dict[str, Any]] = {}
        for key, value in dict(raw).items():
            worker = _redis_text(key)
            payload = _redis_text(value)
            try:
                parsed = json.loads(payload)
            except (TypeError, ValueError, json.JSONDecodeError):
                continue
            if isinstance(parsed, dict):
                decoded[worker] = parsed
        return decoded
