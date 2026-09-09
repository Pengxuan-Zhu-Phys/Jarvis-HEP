#!/usr/bin/env python3
"""Broadcast process boards + inflight key helpers (D26.1)."""

from __future__ import annotations

import threading
import time
from typing import Any

from jarvishep2.queue.redis_queue import (
    ARCHIVER_BOARD_TTL_SEC,
    CodecError,
    INFLIGHT,
    PROC_ARCHIVER,
    PROC_BOARD_TTL_SEC,
    PROC_CHILDREN,
    PROC_CORE,
    PROC_REDIS,
    PROC_ROLES,
    PROC_WORKER,
    SAMPLE_STATS,
    TASK_QUEUE,
    _encode_heartbeat_value,
    _is_redis_timeout,
    _redis_text,
    _warn_control_timeout,
    decode_payload,
)


_BLMOVE_REQUIRED = (
    "Jarvis-HEP D26.1 requires BLMOVE (Redis/Valkey ≥ 6.2); "
    "Ubuntu 22.04 redis-server is 6.0. Install Redis ≥ 6.2 or valkey-server."
)
_BLMOVE_UNSUPPORTED = object()
_INFLIGHT_SCAN_MATCH = "hep:inflight:*"
_INFLIGHT_PREFIX = "hep:inflight:"

# Production occupancy after BLMOVE LEFT LEFT. Bounce extras until one owner.
_ATOMIC_OCCUPANCY_LUA = """
-- KEYS[1]=hep:task_queue  KEYS[2]=hep:inflight:{w}  KEYS[3]=hep:sample:stats
local n_start = redis.call('LLEN', KEYS[2])
if n_start == 0 then
    return {'empty'}
end
while redis.call('LLEN', KEYS[2]) > 1 do
    local extra = redis.call('LPOP', KEYS[2])
    redis.call('LPUSH', KEYS[1], extra)
end
if n_start == 1 then
    redis.call('HINCRBY', KEYS[3], 'running', 1)
end
return {'ok', redis.call('LINDEX', KEYS[2], 0)}
"""

# Test-only non-blocking steal. NEVER a production fallback.
_ATOMIC_STEAL_TO_INFLIGHT_LUA = """
if redis.call('LLEN', KEYS[2]) > 0 then
    return {'occupied'}
end
local payload = redis.call('LPOP', KEYS[1])
if not payload then
    return {'empty'}
end
redis.call('LPUSH', KEYS[2], payload)
redis.call('HINCRBY', KEYS[3], 'running', 1)
return {'ok', payload}
"""

_ATOMIC_ACK_INFLIGHT_LUA = """
-- KEYS[1]=hep:inflight:{w}  ARGV[1]=uuid
local payload = redis.call('LINDEX', KEYS[1], 0)
if not payload then return 0 end
local obj = cjson.decode(payload)
if tostring(obj['uuid']) ~= ARGV[1] then return 0 end
redis.call('LPOP', KEYS[1])
return 1
"""

_ATOMIC_RECLAIM_INFLIGHT_LUA = """
local items = redis.call('LRANGE', KEYS[1], 0, -1)
redis.call('DEL', KEYS[1])
return items
"""


def _board_ttl_sec(role: str, ttl_sec: int | None) -> int:
    ttl = PROC_BOARD_TTL_SEC if ttl_sec is None else max(1, int(ttl_sec))
    if str(role or "").strip().lower() == "archiver":
        ttl = max(ttl, ARCHIVER_BOARD_TTL_SEC)
    return int(ttl)


def _lua_unavailable(exc: BaseException) -> bool:
    text = str(exc).lower()
    return (
        "unknown command 'eval'" in text
        or 'unknown command "eval"' in text
        or "cjson" in text
        or "noscript" in text
    )


def _unknown_command(exc: BaseException) -> bool:
    text = str(exc).lower()
    return "unknown command" in text or "invalid command" in text


def _command_info_has_blmove(info: Any) -> bool:
    if info is None:
        return False
    if isinstance(info, (bytes, str)):
        return "blmove" in str(info).lower()
    if isinstance(info, dict):
        if any(str(key).lower() == "blmove" for key in info):
            return True
        return str(info.get("name") or "").lower() == "blmove"
    if isinstance(info, (list, tuple)):
        if not info:
            return False
        first = info[0]
        if first is None:
            return False
        if isinstance(first, (list, tuple)) and not first:
            return False
        return True
    return True


def _parse_redis_version(raw: Any) -> tuple[int, int] | None:
    text = str(raw or "").strip()
    if not text or text.startswith("<"):
        return None
    parts: list[int] = []
    for token in text.split("."):
        digits = ""
        for char in token:
            if char.isdigit():
                digits += char
            else:
                break
        if not digits:
            break
        parts.append(int(digits))
        if len(parts) >= 2:
            break
    if not parts:
        return None
    if len(parts) == 1:
        parts.append(0)
    return parts[0], parts[1]


class _ProcBoard:
    """Private RedisQueue mixin: broadcast boards + inflight ownership (D26.1)."""

    def _proc_board_key(self, role: str, *, owner_id: str | None = None) -> str:
        role_text = str(role or "").strip().lower()
        if role_text not in PROC_ROLES:
            raise ValueError(f"unknown proc board role: {role!r}")
        if role_text == "core":
            return PROC_CORE
        if role_text == "archiver":
            return PROC_ARCHIVER
        if role_text == "redis":
            return PROC_REDIS
        owner = "" if owner_id is None else str(owner_id).strip()
        if not owner:
            raise ValueError(f"{role_text} board requires owner_id")
        template = PROC_WORKER if role_text == "worker" else PROC_CHILDREN
        return template.format(id=owner)

    def publish_proc_board(
        self,
        role: str,
        /,
        *,
        owner_id: str | None = None,
        ttl_sec: int | None = None,
        **fields: Any,
    ) -> None:
        """HSET overlay + EXPIRE. Raises ValueError on unknown role.

        Worker/children require owner_id. Never called by Monitor.
        Must NOT INCR hep:worker:op_count.
        """
        self._require_client()
        key = self._proc_board_key(role, owner_id=owner_id)
        if "role" not in fields:
            fields["role"] = str(role or "").strip().lower()
        mapping = {
            str(name): _encode_heartbeat_value(value)
            for name, value in fields.items()
            if value is not None
        }
        ttl = _board_ttl_sec(role, ttl_sec)
        try:
            pipe = self._ctrl().pipeline(transaction=True)
            if mapping:
                pipe.hset(key, mapping=mapping)
            pipe.expire(key, ttl)
            pipe.execute()
        except Exception as exc:
            if not _is_redis_timeout(exc):
                raise
            _warn_control_timeout("publish_proc_board", exc)
            self.touch_proc_board(role, owner_id=owner_id, ttl_sec=ttl)

    def read_proc_board(
        self,
        role: str,
        *,
        owner_id: str | None = None,
    ) -> dict[str, Any]:
        """HGETALL; missing key → {} . Never raises on empty."""
        self._require_client()
        key = self._proc_board_key(role, owner_id=owner_id)
        try:
            return dict(self._ctrl().hgetall(key) or {})
        except Exception as exc:
            if not _is_redis_timeout(exc):
                raise
            _warn_control_timeout("read_proc_board", exc)
            return {}

    def read_proc_boards(
        self,
        role: str,
        *,
        owner_ids: list[str],
    ) -> dict[str, dict[str, Any]]:
        """Pipeline HGETALL for worker/children."""
        self._require_client()
        if not owner_ids:
            return {}
        ids = [str(owner_id) for owner_id in owner_ids]
        try:
            pipe = self._ctrl().pipeline(transaction=False)
            for owner_id in ids:
                pipe.hgetall(self._proc_board_key(role, owner_id=owner_id))
            rows = pipe.execute()
        except Exception as exc:
            if not _is_redis_timeout(exc):
                raise
            _warn_control_timeout("read_proc_boards", exc)
            return {owner_id: {} for owner_id in ids}
        return {
            owner_id: dict(row or {})
            for owner_id, row in zip(ids, rows)
        }

    def touch_proc_board(
        self,
        role: str,
        *,
        owner_id: str | None = None,
        ttl_sec: int | None = None,
    ) -> bool:
        """EXPIRE only; False if key missing."""
        self._require_client()
        key = self._proc_board_key(role, owner_id=owner_id)
        try:
            return bool(self._ctrl().expire(key, _board_ttl_sec(role, ttl_sec)))
        except Exception as exc:
            if not _is_redis_timeout(exc):
                raise
            _warn_control_timeout("touch_proc_board", exc)
            return False

    def drop_proc_board(self, role: str, *, owner_id: str | None = None) -> None:
        """Writer shutdown path."""
        self._require_client()
        try:
            self._ctrl().delete(self._proc_board_key(role, owner_id=owner_id))
        except Exception as exc:
            if not _is_redis_timeout(exc):
                raise
            _warn_control_timeout("drop_proc_board", exc)

    def publish_children_board(
        self,
        worker_id: str,
        *,
        file_operation_pid: int | None,
        file_operation_pgid: int | None,
        calc_pgids: list[int],
        reason: str = "heartbeat",
        ttl_sec: int | None = None,
    ) -> None:
        self.publish_proc_board(
            "children",
            owner_id=str(worker_id),
            ttl_sec=ttl_sec,
            role="children",
            file_operation_pid="" if file_operation_pid is None else int(file_operation_pid),
            file_operation_pgid="" if file_operation_pgid is None else int(file_operation_pgid),
            calc_pgids=[int(pid) for pid in calc_pgids],
            updated_reason=str(reason or "heartbeat"),
            ts=time.time(),
        )

    def read_children_board(self, worker_id: str) -> dict[str, Any]:
        return self.read_proc_board("children", owner_id=str(worker_id))

    def _inflight_key(self, worker_id: str) -> str:
        return INFLIGHT.format(worker=str(worker_id))

    def _eval_ctrl(self, script: str, numkeys: int, *keys_and_args: Any) -> Any:
        return self._ctrl().eval(script, numkeys, *keys_and_args)

    def _py_inflight_lock(self) -> threading.Lock:
        lock = getattr(self, "_inflight_py_lock", None)
        if lock is None:
            lock = threading.Lock()
            self._inflight_py_lock = lock
        return lock

    def _steal_allowed(self) -> bool:
        if bool(getattr(self, "_allow_steal_inflight", False)):
            return True
        return self.r is not None and self.r_ctrl is self.r

    def _decode_task_payload(self, raw: Any) -> dict[str, Any] | None:
        if raw is None:
            return None
        if isinstance(raw, dict):
            return dict(raw)
        decoded = decode_payload(raw, codec=self._codec)
        if not isinstance(decoded, dict):
            raise CodecError("task payload must decode to a dict")
        return decoded

    def _task_from_inflight_result(self, result: Any) -> dict[str, Any] | None:
        if result is None:
            return None
        if isinstance(result, dict):
            return dict(result)
        if not isinstance(result, (list, tuple)) or not result:
            return None
        status = _redis_text(result[0]).lower()
        if status in {"empty", "occupied"}:
            return None
        if status != "ok" or len(result) < 2:
            return None
        return self._decode_task_payload(result[1])

    def _occupy_inflight_python(self, inflight_key: str) -> list[Any]:
        with self._py_inflight_lock():
            ctrl = self._ctrl()
            n_start = int(ctrl.llen(inflight_key) or 0)
            if n_start == 0:
                return ["empty"]
            while int(ctrl.llen(inflight_key) or 0) > 1:
                extra = ctrl.lpop(inflight_key)
                if extra is None:
                    break
                ctrl.lpush(TASK_QUEUE, extra)
            if n_start == 1:
                ctrl.hincrby(SAMPLE_STATS, "running", 1)
            return ["ok", ctrl.lindex(inflight_key, 0)]

    def _steal_to_inflight_python(self, inflight_key: str) -> list[Any]:
        with self._py_inflight_lock():
            ctrl = self._ctrl()
            if int(ctrl.llen(inflight_key) or 0) > 0:
                return ["occupied"]
            payload = ctrl.lpop(TASK_QUEUE)
            if payload is None:
                return ["empty"]
            ctrl.lpush(inflight_key, payload)
            ctrl.hincrby(SAMPLE_STATS, "running", 1)
            return ["ok", payload]

    def _ack_inflight_python(self, inflight_key: str, uuid: str) -> int:
        with self._py_inflight_lock():
            ctrl = self._ctrl()
            raw = ctrl.lindex(inflight_key, 0)
            if raw is None:
                return 0
            try:
                payload = self._decode_task_payload(raw)
            except (CodecError, TypeError, ValueError):
                return 0
            if payload is None or str(payload.get("uuid")) != uuid:
                return 0
            popped = ctrl.lpop(inflight_key)
            if popped is None:
                return 0
            return 1

    def _reclaim_inflight_python(self, inflight_key: str) -> list[Any]:
        with self._py_inflight_lock():
            ctrl = self._ctrl()
            items = list(ctrl.lrange(inflight_key, 0, -1) or [])
            ctrl.delete(inflight_key)
            return items

    def _run_occupancy(self, inflight_key: str) -> dict[str, Any] | None:
        try:
            result = self._eval_ctrl(
                _ATOMIC_OCCUPANCY_LUA,
                3,
                TASK_QUEUE,
                inflight_key,
                SAMPLE_STATS,
            )
        except Exception as exc:
            if not _lua_unavailable(exc):
                raise
            result = self._occupy_inflight_python(inflight_key)
        return self._task_from_inflight_result(result)

    def _steal_to_inflight(self, inflight_key: str) -> dict[str, Any] | None:
        try:
            result = self._eval_ctrl(
                _ATOMIC_STEAL_TO_INFLIGHT_LUA,
                3,
                TASK_QUEUE,
                inflight_key,
                SAMPLE_STATS,
            )
        except Exception as exc:
            if not _lua_unavailable(exc):
                raise
            result = self._steal_to_inflight_python(inflight_key)
        return self._task_from_inflight_result(result)

    def _blmove_to_inflight(self, inflight_key: str, timeout: int) -> Any:
        blmove = getattr(self.r, "blmove", None)
        if not callable(blmove):
            return _BLMOVE_UNSUPPORTED
        try:
            return blmove(
                TASK_QUEUE,
                inflight_key,
                max(0, int(timeout)),
                src="LEFT",
                dest="LEFT",
            )
        except Exception as exc:
            if _is_redis_timeout(exc):
                return None
            if _unknown_command(exc):
                return _BLMOVE_UNSUPPORTED
            raise

    def require_blmove(self) -> None:
        """COMMAND INFO BLMOVE (or INFO redis_version ≥ 6.2). Fail fast."""
        self._require_client()
        client = self.r
        for name in ("BLMOVE", "blmove"):
            try:
                info = client.execute_command("COMMAND", "INFO", name)
            except Exception:
                info = None
            if _command_info_has_blmove(info):
                return
        version: tuple[int, int] | None = None
        for args in (("server",), ()):
            try:
                raw_info = client.info(*args)
            except Exception:
                continue
            if isinstance(raw_info, dict):
                version = _parse_redis_version(raw_info.get("redis_version"))
                if version is not None:
                    break
        if version is not None and version >= (6, 2):
            return
        raise RuntimeError(_BLMOVE_REQUIRED)

    def pull_task_to_inflight(
        self,
        worker_id: str,
        timeout: int = 5,
    ) -> dict[str, Any] | None:
        """Blocking BLMOVE LEFT LEFT then occupancy Lua.

        fakeredis without BLMOVE uses the test-only steal Lua (non-blocking).
        None only if the queue was empty. Occupancy returns the remaining head.
        """
        self._require_client()
        inflight_key = self._inflight_key(worker_id)
        moved = self._blmove_to_inflight(inflight_key, timeout)
        if moved is _BLMOVE_UNSUPPORTED:
            if not self._steal_allowed():
                raise RuntimeError(_BLMOVE_REQUIRED)
            return self._steal_to_inflight(inflight_key)
        if moved is None:
            return self.occupy_inflight_task(worker_id)
        return self._run_occupancy(inflight_key)

    def occupy_inflight_task(self, worker_id: str) -> dict[str, Any] | None:
        """Bounce extras until LLEN<=1, then return the remaining head.

        LLEN==1 does not run occupancy (would double-count ``running``).
        """
        self._require_client()
        inflight_key = self._inflight_key(worker_id)
        llen = int(self._ctrl().llen(inflight_key) or 0)
        if llen <= 0:
            return None
        if llen > 1:
            return self._run_occupancy(inflight_key)
        return self.get_inflight_task(worker_id)

    def get_inflight_task(self, worker_id: str) -> dict[str, Any] | None:
        """LINDEX 0; do not consume."""
        self._require_client()
        raw = self._ctrl().lindex(self._inflight_key(worker_id), 0)
        if raw is None:
            return None
        try:
            return self._decode_task_payload(raw)
        except (CodecError, TypeError, ValueError):
            return None

    def ack_inflight_task(self, worker_id: str, uuid: str) -> bool:
        """Uuid-matched LPOP of the list head. Never DEL."""
        self._require_client()
        inflight_key = self._inflight_key(worker_id)
        wanted = str(uuid)
        try:
            result = self._eval_ctrl(_ATOMIC_ACK_INFLIGHT_LUA, 1, inflight_key, wanted)
        except Exception as exc:
            if not _lua_unavailable(exc):
                raise
            result = self._ack_inflight_python(inflight_key, wanted)
        try:
            return bool(int(result or 0))
        except (TypeError, ValueError):
            return False

    def reclaim_inflight_task(self, worker_id: str) -> list[dict[str, Any]]:
        """Drain the whole inflight list. Caller retries and push_task each."""
        self._require_client()
        inflight_key = self._inflight_key(worker_id)
        try:
            items = self._eval_ctrl(_ATOMIC_RECLAIM_INFLIGHT_LUA, 1, inflight_key)
        except Exception as exc:
            if not _lua_unavailable(exc):
                raise
            items = self._reclaim_inflight_python(inflight_key)
        payloads: list[dict[str, Any]] = []
        for raw in items or []:
            try:
                decoded = self._decode_task_payload(raw)
            except (CodecError, TypeError, ValueError):
                continue
            if decoded is not None:
                payloads.append(decoded)
        return payloads

    def list_inflight_worker_ids(self) -> list[str]:
        """SCAN match hep:inflight:* (resume only)."""
        self._require_client()
        ids: list[str] = []
        seen: set[str] = set()
        for key in self._ctrl().scan_iter(match=_INFLIGHT_SCAN_MATCH):
            text = _redis_text(key)
            if not text.startswith(_INFLIGHT_PREFIX):
                continue
            worker_id = text[len(_INFLIGHT_PREFIX) :]
            if not worker_id or worker_id in seen:
                continue
            seen.add(worker_id)
            ids.append(worker_id)
        return ids
