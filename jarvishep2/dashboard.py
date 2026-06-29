#!/usr/bin/env python3
"""Read-only monitor snapshot reader for Jarvis-HEP V2 (WP-D5.2)."""

from __future__ import annotations

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
        return cls(
            sampler={
                "task_queue_length": int(snapshot.get("task_queue_length", 0) or 0),
            },
            factory={
                "workers_alive": int(snapshot.get("workers_alive", 0) or 0),
                "workers_total": int(snapshot.get("workers_total", 0) or 0),
            },
            workers=list(snapshot.get("workers") or []),
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
            timestamp=_coerce_float(snapshot.get("timestamp")),
        )


class SnapshotReader:
    """Build a :class:`MonitorView` without issuing Redis writes."""

    def __init__(self, source: TaskFactory | RedisQueue) -> None:
        self._source = source

    def read(self) -> MonitorView:
        if isinstance(self._source, TaskFactory):
            return MonitorView.from_snapshot(self._source.get_monitor_snapshot())
        raw = self._source.snapshot_raw()
        raw.setdefault("workers", [])
        raw.setdefault("workers_alive", 0)
        raw.setdefault("workers_total", 0)
        raw["timestamp"] = None
        return MonitorView.from_snapshot(raw)


def attach_reader(
    *,
    factory: TaskFactory | None = None,
    redis: RedisQueue | None = None,
) -> SnapshotReader:
    if factory is not None:
        return SnapshotReader(factory)
    if redis is not None:
        return SnapshotReader(redis)
    raise ValueError("attach_reader requires a TaskFactory or RedisQueue")


def format_monitor_view(view: MonitorView) -> str:
    lines = [
        "Jarvis2 Monitor Snapshot",
        f"timestamp: {view.timestamp}",
        f"workers: {view.factory.get('workers_alive', 0)}/{view.factory.get('workers_total', 0)} alive",
        f"task_queue_length: {view.queues.get('task_queue_length', 0)}",
        f"archive_queue_length: {view.queues.get('archive_queue_length', 0)}",
        f"sample_stats: {view.samples}",
        f"calculator_status: {view.calculators}",
        f"op_counts: {view.op_counts}",
    ]
    for worker in view.workers:
        lines.append(
            "worker "
            f"{worker.get('worker_id')}: pid={worker.get('pid')} alive={worker.get('alive')}"
        )
    return "\n".join(lines) + "\n"


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