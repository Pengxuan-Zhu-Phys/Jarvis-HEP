#!/usr/bin/env python3
"""Broadcast process boards + inflight key helpers (D26.1)."""

from __future__ import annotations

import time
from typing import Any

from jarvishep2.queue.redis_queue import (
    ARCHIVER_BOARD_TTL_SEC,
    PROC_ARCHIVER,
    PROC_BOARD_TTL_SEC,
    PROC_CHILDREN,
    PROC_CORE,
    PROC_REDIS,
    PROC_ROLES,
    PROC_WORKER,
    _encode_heartbeat_value,
    _is_redis_timeout,
    _warn_control_timeout,
)


def _board_ttl_sec(role: str, ttl_sec: int | None) -> int:
    ttl = PROC_BOARD_TTL_SEC if ttl_sec is None else max(1, int(ttl_sec))
    if str(role or "").strip().lower() == "archiver":
        ttl = max(ttl, ARCHIVER_BOARD_TTL_SEC)
    return int(ttl)


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
