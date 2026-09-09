#!/usr/bin/env python3
"""Control lock + worker heartbeat store (D25.4)."""

from __future__ import annotations

import json
import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.queue.redis_queue import (
    CONTROL_LOCK,
    CONTROL_LOCK_TTL_SEC,
    OP_COUNT,
    PROC_BOARD_TTL_SEC,
    PROC_CHILDREN,
    PROC_WORKER,
    WORKER_STATUS,
    _encode_heartbeat_value,
    _is_redis_timeout,
    _redis_text,
    _warn_control_timeout,
    decode_payload,
    encode_payload,
)


def classify_control_lock(current: Any, expected: str) -> str:
    """Classify a control-lock GET against the owner we expect.

    ``ok`` — we still hold it. ``missing`` — key is gone (TTL blip / eviction);
    Core's refresh is designed to reclaim an expired *own* lease, so this is
    transient. ``stolen`` — a different owner holds it (another Jarvis).
    """
    expected_text = str(expected or "").strip()
    if not expected_text:
        return "ok"
    if current is None:
        return "missing"
    current_text = _redis_text(current).strip()
    if not current_text:
        return "missing"
    if current_text != expected_text:
        return "stolen"
    return "ok"


def control_lock_missing_grace_sec(ttl_sec: int | float | None = None) -> float:
    """How long Archiver/Workers tolerate a missing lock before exiting."""
    ttl = CONTROL_LOCK_TTL_SEC if ttl_sec is None else max(5, int(ttl_sec))
    return float(max(60, 2 * ttl))


def next_control_lock_watch(
    status: str,
    *,
    missing_since: float | None,
    now: float,
    grace_sec: float,
) -> tuple[str, float | None]:
    """Advance lock-watch state.

    Returns ``(decision, missing_since)`` where decision is ``ok``,
    ``continue_missing``, ``shutdown_missing``, or ``shutdown_stolen``.
    """
    if status == "ok":
        return "ok", None
    if status == "stolen":
        return "shutdown_stolen", missing_since
    if missing_since is None:
        return "continue_missing", now
    if (now - missing_since) >= max(1.0, float(grace_sec)):
        return "shutdown_missing", missing_since
    return "continue_missing", missing_since

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


def _worker_board_fields(fields: Mapping[str, Any]) -> dict[str, Any]:
    """Display overlay for hep:proc:worker:{id}; never copies ownership blobs."""
    board: dict[str, Any] = {"role": "worker"}
    status = fields.get("status")
    if status is not None:
        board["status"] = status
    pid = fields.get("pid")
    if pid is not None:
        board["pid"] = pid
    ts = fields.get("ts", fields.get("last_heartbeat"))
    if ts is not None:
        board["ts"] = ts
    current = fields.get("current_sample")
    if current is None or current == "":
        current = fields.get("uuid")
    board["current_uuid"] = "" if current is None or current == "" else str(current)
    fo_pid = fields.get("file_operation_pid")
    if fo_pid is not None and fo_pid != "":
        board["file_operation_pid"] = fo_pid
    if "held_calc_n" in fields and fields.get("held_calc_n") is not None:
        try:
            board["held_calc_n"] = int(fields["held_calc_n"])
        except (TypeError, ValueError):
            pass
    else:
        raw = fields.get("held_calc_packs")
        if raw is None or raw == "":
            board["held_calc_n"] = 0
        elif isinstance(raw, Mapping):
            board["held_calc_n"] = len([value for value in raw.values() if value])
        else:
            try:
                decoded = json.loads(str(raw))
            except (TypeError, ValueError, json.JSONDecodeError):
                decoded = None
            if isinstance(decoded, dict):
                board["held_calc_n"] = len(
                    [value for value in decoded.values() if value]
                )
    interval = fields.get("heartbeat_interval_sec")
    if interval is not None:
        board["heartbeat_interval_sec"] = interval
    return board


# Ownership is hep:inflight:{id}; never copy the task payload onto status/board.
_HEARTBEAT_STATUS_OMIT = frozenset({"current_task", "u_coords", "execution_plan"})


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
        return bool(self._ctrl().set(CONTROL_LOCK, owner_text, nx=True, ex=ttl))

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
        ctrl = self._ctrl()
        try:
            result = ctrl.eval(
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
            current = ctrl.get(CONTROL_LOCK)
            if current is None:
                return bool(ctrl.set(CONTROL_LOCK, owner_text, nx=True, ex=ttl))
            if _redis_text(current) != owner_text:
                return False
            return bool(ctrl.expire(CONTROL_LOCK, ttl))
        return bool(int(result or 0))

    def release_control_lock(self, owner: str) -> bool:
        """Drop the control lease if we still own it."""
        self._require_client()
        owner_text = str(owner or "").strip()
        if not owner_text:
            return False
        ctrl = self._ctrl()
        try:
            result = ctrl.eval(
                _ATOMIC_RELEASE_CONTROL_LOCK_LUA,
                1,
                CONTROL_LOCK,
                owner_text,
            )
        except Exception as exc:
            if "unknown command 'eval'" not in str(exc).lower():
                raise
            current = ctrl.get(CONTROL_LOCK)
            if current is None or _redis_text(current) != owner_text:
                return False
            return bool(ctrl.delete(CONTROL_LOCK))
        return bool(int(result or 0))

    def get_control_lock_owner(self) -> str | None:
        self._require_client()
        try:
            value = self._ctrl().get(CONTROL_LOCK)
        except Exception as exc:
            if not _is_redis_timeout(exc):
                raise
            _warn_control_timeout("get_control_lock_owner", exc)
            return None
        if value is None:
            return None
        text = _redis_text(value).strip()
        return text or None

    def heartbeat(
        self,
        worker_id: str,
        *,
        board_ttl_sec: int | None = None,
        **fields: Any,
    ) -> None:
        self._require_client()
        if "last_heartbeat" not in fields and "ts" in fields:
            fields["last_heartbeat"] = fields["ts"]
        overlay_children = "file_operation_pgid" in fields or "calc_pgids" in fields
        children_pgid = fields.pop("file_operation_pgid", None)
        calc_pgids = fields.pop("calc_pgids", None)
        children_reason = fields.pop("children_updated_reason", None)
        wid = str(worker_id)
        status_key = WORKER_STATUS.format(id=wid)
        board_key = PROC_WORKER.format(id=wid)
        children_key = PROC_CHILDREN.format(id=wid)
        mapping = {
            k: _encode_heartbeat_value(v)
            for k, v in fields.items()
            if k not in _HEARTBEAT_STATUS_OMIT
        }
        board_mapping = {
            k: _encode_heartbeat_value(v)
            for k, v in _worker_board_fields(fields).items()
        }
        ttl = (
            max(1, int(board_ttl_sec))
            if board_ttl_sec is not None
            else PROC_BOARD_TTL_SEC
        )
        try:
            pipe = self._ctrl().pipeline(transaction=True)
            if mapping:
                pipe.hset(status_key, mapping=mapping)
            pipe.hdel(status_key, *_HEARTBEAT_STATUS_OMIT)
            if board_mapping:
                pipe.hset(board_key, mapping=board_mapping)
            pipe.expire(board_key, ttl)
            if overlay_children:
                fo_pid = fields.get("file_operation_pid")
                try:
                    fo_pgid_value: Any = (
                        "" if children_pgid is None or children_pgid == "" else int(children_pgid)
                    )
                    if fo_pgid_value != "" and int(fo_pgid_value) <= 0:
                        fo_pgid_value = ""
                except (TypeError, ValueError):
                    fo_pgid_value = ""
                if isinstance(calc_pgids, (list, tuple)):
                    calc_list = [int(pid) for pid in calc_pgids]
                else:
                    calc_list = []
                children_mapping = {
                    k: _encode_heartbeat_value(v)
                    for k, v in {
                        "role": "children",
                        "file_operation_pid": "" if fo_pid is None else fo_pid,
                        "file_operation_pgid": fo_pgid_value,
                        "calc_pgids": calc_list,
                        "updated_reason": str(children_reason or "heartbeat"),
                        "ts": fields.get("ts") or fields.get("last_heartbeat") or time.time(),
                    }.items()
                }
                pipe.hset(children_key, mapping=children_mapping)
            # Same TTL as the worker board so children cannot evaporate first.
            pipe.expire(children_key, ttl)
            pipe.incr(OP_COUNT.format(kind="worker"))
            pipe.execute()
        except Exception as exc:
            if not _is_redis_timeout(exc):
                raise
            _warn_control_timeout("heartbeat", exc)
            self.touch_proc_board("worker", owner_id=wid, ttl_sec=ttl)
            self.touch_proc_board("children", owner_id=wid, ttl_sec=ttl)

    def encode_task_for_heartbeat(self, task: Mapping[str, Any]) -> str:
        """Serialize an in-flight task for the Worker heartbeat hash."""
        return encode_payload(dict(task), codec=self._codec)

    def decode_heartbeat_task(self, heartbeat: Mapping[str, Any]) -> dict[str, Any] | None:
        """Decode a leftover heartbeat task blob. Heartbeat no longer stores it."""
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
