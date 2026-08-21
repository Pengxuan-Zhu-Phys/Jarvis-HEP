#!/usr/bin/env python3
"""Control lock + worker heartbeat store (D25.4)."""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

from jarvishep2.redis_queue import (
    CONTROL_LOCK,
    CONTROL_LOCK_TTL_SEC,
    OP_COUNT,
    WORKER_STATUS,
    _encode_heartbeat_value,
    _redis_text,
    decode_payload,
    encode_payload,
)

_ATOMIC_REFRESH_CONTROL_LOCK_LUA = """
local current = redis.call('GET', KEYS[1])
if not current then
    redis.call('SET', KEYS[1], ARGV[1], 'EX', ARGV[2])
    return 1
end
if current ~= ARGV[1] then
    return 0
end
redis.call('EXPIRE', KEYS[1], ARGV[2])
return 1
"""

_ATOMIC_RELEASE_CONTROL_LOCK_LUA = """
if redis.call('GET', KEYS[1]) ~= ARGV[1] then
    return 0
end
return redis.call('DEL', KEYS[1])
"""


class _ControlAndHeartbeat:
    """Private RedisQueue mixin (D25.4)."""

    def claim_control_lock(
        self,
        owner: str,
        *,
        ttl_sec: int = CONTROL_LOCK_TTL_SEC,
    ) -> bool:
        """Acquire exclusive control lease (SET NX EX). False if another owner holds it."""
        self._require_client()
        owner_text = str(owner or "").strip()
        if not owner_text:
            raise ValueError("control lock owner is required")
        ttl = max(5, int(ttl_sec))
        return bool(self.r.set(CONTROL_LOCK, owner_text, nx=True, ex=ttl))

    def refresh_control_lock(
        self,
        owner: str,
        *,
        ttl_sec: int = CONTROL_LOCK_TTL_SEC,
    ) -> bool:
        """Extend the lease only if *owner* still holds it."""
        self._require_client()
        owner_text = str(owner or "").strip()
        if not owner_text:
            raise ValueError("control lock owner is required")
        ttl = max(5, int(ttl_sec))
        try:
            result = self.r.eval(
                _ATOMIC_REFRESH_CONTROL_LOCK_LUA,
                1,
                CONTROL_LOCK,
                owner_text,
                ttl,
            )
        except Exception as exc:
            # fakeredis versions used by the offline suite may omit EVAL.
            # Real Redis always takes the atomic Lua path above.
            if "unknown command 'eval'" not in str(exc).lower():
                raise
            current = self.r.get(CONTROL_LOCK)
            if current is None:
                return bool(self.r.set(CONTROL_LOCK, owner_text, nx=True, ex=ttl))
            if _redis_text(current) != owner_text:
                return False
            return bool(self.r.expire(CONTROL_LOCK, ttl))
        return bool(int(result or 0))

    def release_control_lock(self, owner: str) -> bool:
        """Drop the control lease if we still own it."""
        self._require_client()
        owner_text = str(owner or "").strip()
        if not owner_text:
            return False
        try:
            result = self.r.eval(
                _ATOMIC_RELEASE_CONTROL_LOCK_LUA,
                1,
                CONTROL_LOCK,
                owner_text,
            )
        except Exception as exc:
            if "unknown command 'eval'" not in str(exc).lower():
                raise
            current = self.r.get(CONTROL_LOCK)
            if current is None or _redis_text(current) != owner_text:
                return False
            return bool(self.r.delete(CONTROL_LOCK))
        return bool(int(result or 0))

    def get_control_lock_owner(self) -> str | None:
        self._require_client()
        value = self.r.get(CONTROL_LOCK)
        return None if value is None else str(value)

    def heartbeat(self, worker_id: str, **fields: Any) -> None:
        self._require_client()
        if "last_heartbeat" not in fields and "ts" in fields:
            fields["last_heartbeat"] = fields["ts"]
        key = WORKER_STATUS.format(id=worker_id)
        mapping = {k: _encode_heartbeat_value(v) for k, v in fields.items()}
        pipe = self.r.pipeline(transaction=True)
        if mapping:
            pipe.hset(key, mapping=mapping)
        pipe.incr(OP_COUNT.format(kind="worker"))
        pipe.execute()

    def encode_task_for_heartbeat(self, task: Mapping[str, Any]) -> str:
        """Serialize an in-flight task for the Worker heartbeat hash."""
        return encode_payload(dict(task), codec=self._codec)

    def decode_heartbeat_task(self, heartbeat: Mapping[str, Any]) -> dict[str, Any] | None:
        """Decode a Worker heartbeat's serialized in-flight task payload."""
        raw = heartbeat.get("current_task")
        if raw is None or raw == "":
            return None
        if isinstance(raw, Mapping):
            return dict(raw)
        decoded = decode_payload(str(raw), codec=self._codec)
        return dict(decoded) if isinstance(decoded, dict) else None

    def decode_heartbeat_subprocess_pids(self, heartbeat: Mapping[str, Any]) -> list[int]:
        """Decode active calculator subprocess PIDs from a Worker heartbeat."""
        raw = heartbeat.get("active_subprocess_pids")
        if raw is None or raw == "":
            return []
        if isinstance(raw, (list, tuple)):
            decoded: Any = list(raw)
        else:
            try:
                decoded = json.loads(str(raw))
            except (TypeError, ValueError, json.JSONDecodeError):
                return []
        if not isinstance(decoded, list):
            return []
        pids: list[int] = []
        for item in decoded:
            try:
                pid = int(item)
            except (TypeError, ValueError):
                continue
            if pid > 0:
                pids.append(pid)
        return pids

    def decode_heartbeat_held_packs(self, heartbeat: Mapping[str, Any]) -> dict[str, str]:
        """Decode held calculator slots from a Worker heartbeat hash."""
        raw = heartbeat.get("held_calc_packs")
        if raw is None or raw == "":
            return {}
        if isinstance(raw, Mapping):
            return {str(key): str(value) for key, value in raw.items() if value}
        try:
            decoded = json.loads(str(raw))
        except (TypeError, ValueError, json.JSONDecodeError):
            return {}
        if not isinstance(decoded, dict):
            return {}
        return {str(key): str(value) for key, value in decoded.items() if value}
