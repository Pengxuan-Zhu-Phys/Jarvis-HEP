#!/usr/bin/env python3
"""Read-only monitor snapshot reader for Jarvis-HEP V2 (WP-D5.2)."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import RedisQueue


@dataclass
class MonitorView:
    sampler: dict[str, Any] = field(default_factory=dict)
    factory: dict[str, Any] = field(default_factory=dict)
    workers: list[dict[str, Any]] = field(default_factory=list)
    calculators: dict[str, Any] = field(default_factory=dict)
    samples: dict[str, Any] = field(default_factory=dict)
    resources: dict[str, Any] = field(default_factory=dict)
    queues: dict[str, Any] = field(default_factory=dict)
    op_counts: dict[str, int] = field(default_factory=dict)
    proc_core: dict[str, Any] = field(default_factory=dict)
    proc_archiver: dict[str, Any] = field(default_factory=dict)
    proc_redis: dict[str, Any] = field(default_factory=dict)
    timestamp: float | None = None

    def has_active_scan(self) -> bool:
        if self.workers:
            return True
        if int(self.queues.get("task_queue_length", 0) or 0) > 0:
            return True
        if int(self.queues.get("archive_queue_length", 0) or 0) > 0:
            return True
        stats = self.samples or {}
        for key in ("running", "completed", "failed"):
            if int(stats.get(key, 0) or 0) > 0:
                return True
        return False

    @classmethod
    def from_snapshot(cls, snapshot: dict[str, Any]) -> MonitorView:
        sample_stats = dict(snapshot.get("sample_stats") or {})
        observed_at = _coerce_float(snapshot.get("timestamp"))
        return cls(
            sampler={
                "task_queue_length": int(snapshot.get("task_queue_length", 0) or 0),
            },
            factory={
                "workers_alive": int(snapshot.get("workers_alive", 0) or 0),
                "workers_total": int(snapshot.get("workers_total", 0) or 0),
            },
            workers=_project_workers(snapshot, now=observed_at),
            calculators=dict(snapshot.get("calculator_status") or {}),
            samples=sample_stats,
            resources={},
            queues={
                "task_queue_length": int(snapshot.get("task_queue_length", 0) or 0),
                "archive_queue_length": int(snapshot.get("archive_queue_length", 0) or 0),
            },
            op_counts={
                str(key): int(value or 0)
                for key, value in dict(snapshot.get("op_counts") or {}).items()
            },
            proc_core=dict(snapshot.get("proc_core") or {}),
            proc_archiver=dict(snapshot.get("proc_archiver") or {}),
            proc_redis=dict(snapshot.get("proc_redis") or {}),
            timestamp=observed_at,
        )


class SnapshotReader:
    """Build a :class:`MonitorView` without issuing Redis writes.

    Read-only is a client convention (no SET/HSET), not a Redis ACL.
    """

    def __init__(
        self,
        source: TaskFactory | RedisQueue,
        *,
        owner_ids: list[str] | None = None,
    ) -> None:
        self._source = source
        self._owner_ids = owner_ids

    def read(self) -> MonitorView:
        if isinstance(self._source, TaskFactory):
            return MonitorView.from_snapshot(self._source.get_monitor_snapshot())
        raw = self._source.snapshot_raw(owner_ids=self._owner_ids)
        raw.setdefault("workers", [])
        raw.setdefault("workers_alive", 0)
        raw.setdefault("workers_total", 0)
        raw["timestamp"] = None
        return MonitorView.from_snapshot(raw)


def attach_reader(
    *,
    factory: TaskFactory | None = None,
    redis: RedisQueue | None = None,
    owner_ids: list[str] | None = None,
) -> SnapshotReader:
    if factory is not None:
        return SnapshotReader(factory, owner_ids=owner_ids)
    if redis is not None:
        return SnapshotReader(redis, owner_ids=owner_ids)
    raise ValueError("attach_reader requires a TaskFactory or RedisQueue")


def format_monitor_view(view: MonitorView) -> str:
    lines = [
        "Jarvis Monitor Snapshot",
        f"timestamp: {view.timestamp}",
        f"workers: {view.factory.get('workers_alive', 0)}/{view.factory.get('workers_total', 0)} alive",
        f"task_queue_length: {view.queues.get('task_queue_length', 0)}",
        f"archive_queue_length: {view.queues.get('archive_queue_length', 0)}",
        f"sample_stats: {view.samples}",
        f"calculator_status: {view.calculators}",
        f"op_counts: {view.op_counts}",
        f"scan_mode: {view.proc_core.get('scan_mode')}",
    ]
    for worker in view.workers:
        lines.append(
            "worker "
            f"{worker.get('worker_id')}: pid={worker.get('pid')} "
            f"status={worker.get('status')} "
            f"current_uuid={worker.get('current_uuid')} "
            f"heartbeat_age_s={worker.get('heartbeat_age_s')} "
            f"alive={worker.get('alive')}"
        )
    return "\n".join(lines) + "\n"


def _project_workers(
    snapshot: Mapping[str, Any],
    *,
    now: float | None,
) -> list[dict[str, Any]]:
    heartbeats = dict(snapshot.get("worker_heartbeats") or {})
    proc_workers = dict(snapshot.get("proc_workers") or {})
    rows = list(snapshot.get("workers") or [])
    if not rows and proc_workers:
        rows = [{"worker_id": worker_id, **dict(row or {})} for worker_id, row in proc_workers.items()]
    projected: list[dict[str, Any]] = []
    for row in rows:
        worker_id = str(row.get("worker_id") if row.get("worker_id") is not None else "")
        merged = dict(row)
        merged.update(dict(heartbeats.get(worker_id) or {}))
        merged.update(dict(proc_workers.get(worker_id) or {}))
        projected.append(_worker_display(merged, now=now))
    return projected


def _worker_display(row: Mapping[str, Any], *, now: float | None) -> dict[str, Any]:
    ts = row.get("ts", row.get("last_heartbeat"))
    age = None
    if ts not in (None, "") and now is not None:
        try:
            age = max(0.0, float(now) - float(ts))
        except (TypeError, ValueError):
            age = None
    current = row.get("current_uuid")
    if current in (None, ""):
        current = row.get("current_sample")
    return {
        "worker_id": row.get("worker_id"),
        "pid": row.get("pid"),
        "status": row.get("status") or row.get("state"),
        "current_uuid": "" if current in (None, "") else str(current),
        "heartbeat_age_s": age,
        "alive": row.get("alive"),
    }


def _coerce_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


__all__ = [
    "MonitorView",
    "SnapshotReader",
    "attach_reader",
    "format_monitor_view",
]