"""Translate live Collector snapshots into the locked Overview presentation model."""

from __future__ import annotations

from datetime import UTC, datetime
import time
from typing import Any

from jarvishep2.monitor.collector import CollectorFrame
from jarvishep2.monitor.scans import ScanChoice
from jarvishep2.monitor.simu import OverviewFrame


def _int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError, OverflowError):
        return default


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError, OverflowError):
        return None


def _duration(seconds: float | None) -> str:
    if seconds is None:
        return "—"
    total = max(0, int(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def _age_fact(board: dict[str, Any], now: float) -> tuple[str, str]:
    timestamp = _float(board.get("ts"))
    if timestamp is None:
        return "heartbeat unknown", "hb —"
    age = max(0.0, now - timestamp)
    if age < 60:
        return f"hb {age:0.1f}s", f"{age:0.1f}s"
    minutes = age / 60.0
    return f"hb {minutes:0.1f}m", f"{minutes:0.1f}m"


class LiveOverviewProjector:
    """Stateful, read-only projection with window-sized local spark histories."""

    def __init__(
        self,
        choice: ScanChoice,
        *,
        redis_endpoint: str = "unavailable",
    ) -> None:
        self.choice = choice
        self.redis_endpoint = str(redis_endpoint or "unavailable")
        self._history_widths = (1, 1)
        self._samples: list[float | None] = [None]
        self._queues: list[float | None] = [None]
        self._last_health: tuple[tuple[str, str, str, str], ...] = ()

    def set_history_widths(self, samples: int, queues: int) -> None:
        widths = (max(1, int(samples)), max(1, int(queues)))
        if widths == self._history_widths:
            return
        self._history_widths = widths
        self._samples = self._resize(self._samples, widths[0])
        self._queues = self._resize(self._queues, widths[1])

    @staticmethod
    def _resize(values: list[float | None], width: int) -> list[float | None]:
        return [*values[:width], *([None] * max(0, width - len(values)))]

    @staticmethod
    def _push(values: list[float | None], value: float | None, width: int) -> None:
        values.insert(0, value)
        del values[width:]

    def project(self, frame: CollectorFrame) -> OverviewFrame:
        now = time.time()
        if not frame.stale:
            self._push(
                self._samples,
                frame.sample_rate_per_sec,
                self._history_widths[0],
            )
            self._push(
                self._queues,
                None
                if frame.queue_block.task_depth is None
                else float(frame.queue_block.task_depth),
                self._history_widths[1],
            )

        stats = frame.sample_stats
        done = max(0, _int(stats.get("completed")))
        running = max(0, _int(stats.get("running")))
        failed = max(0, _int(stats.get("failed")))
        sampler = frame.sampler_display
        progress = sampler.progress
        target = max(0, _int(progress.target))
        current = max(0, _int(progress.current))
        task_q = max(0, _int(frame.queue_block.task_depth))
        archive_q = max(0, _int(frame.queue_block.archive_depth))
        feedback_q = max(0, _int(frame.queue_block.feedback_depth))

        worker = frame.worker_block
        workers_total = max(0, _int(worker.expected))
        workers_alive = max(0, _int(worker.present))
        stale_workers = max(0, _int(worker.heartbeat_stale))
        busy = max(0, _int(worker.busy))

        host = frame.host
        cpu_reading = _float(host.get("cpu_percent"))
        memory_used = max(0, _int(host.get("memory_used")))
        memory_total = max(0, _int(host.get("memory_total")))
        resources_available = (
            bool(host.get("available"))
            and cpu_reading is not None
            and memory_total > 0
        )
        cpu = cpu_reading or 0.0
        processes = [row for row in host.get("processes", ()) if isinstance(row, dict)]
        fds = sum(max(0, _int(row.get("fds"))) for row in processes)
        fd_limits = [
            max(0, _int(row.get("fd_limit")))
            for row in processes
            if _int(row.get("fd_limit")) > 0
        ]
        fds_limit = sum(fd_limits) if fd_limits else 0

        core = frame.proc_core
        mode = str(core.get("scan_mode") or core.get("status") or "unknown").lower()
        started_at = _float(core.get("started_at"))
        elapsed = _duration(None if started_at is None else now - started_at)
        started = (
            "—"
            if started_at is None
            else datetime.fromtimestamp(started_at, tz=UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
        )
        health = self._health_items(frame, now=now)
        if not frame.stale:
            self._last_health = health
        elif self._last_health:
            health = self._last_health

        calculators = tuple(
            (
                str(name),
                max(0, _int(row.get("busy"))),
                max(0, _int(row.get("slots"))),
            )
            for name, row in frame.occupancy.items()
        )
        generation = current if progress.kind == "generation" else 0
        max_generations = target if progress.kind == "generation" else 0
        metrics = dict(sampler.metrics)
        return OverviewFrame(
            scan=self.choice.name,
            mode="stale" if frame.stale else mode,
            elapsed=elapsed,
            ref=self.choice.name,
            redis=self.redis_endpoint,
            hz="2.0 Hz",
            target=target,
            done=done,
            running=running,
            failed=failed,
            rate=frame.sample_rate,
            eta=str(progress.eta or "—"),
            avg="—",
            task_q=task_q,
            archive_q=archive_q,
            feedback_q=feedback_q,
            queue_block=frame.queue_block,
            worker_block=worker,
            workers_alive=workers_alive,
            workers_total=workers_total,
            stale=stale_workers,
            busy=busy,
            cpu=cpu,
            mem_g=memory_used / 1024**3,
            mem_total_g=memory_total / 1024**3,
            method=sampler.method,
            run_id=str(core.get("run_id") or "—"),
            host=str(core.get("host") or "—"),
            started=started,
            calculators=calculators,
            calculator_block=frame.calculator_block,
            spark_samples=tuple(self._samples),
            spark_queue=tuple(self._queues),
            fds=fds,
            fds_limit=fds_limit,
            sampler_generation=generation,
            sampler_max_generations=max_generations,
            sampler_radius=_float(metrics.get("radius")) or 0.0,
            sampler_core=max(0, _int(metrics.get("core"))),
            sampler_open=max(0, _int(metrics.get("open"))),
            sampler_state=sampler.state,
            sampler_saved=max(0, _int(sampler.saved)),
            sampler=sampler,
            resources_available=resources_available,
            health_items=health,
        )

    def _health_items(
        self,
        frame: CollectorFrame,
        *,
        now: float,
    ) -> tuple[tuple[str, str, str, str], ...]:
        core = frame.proc_core
        mode = str(core.get("scan_mode") or core.get("status") or "unknown").lower()
        outstanding = (
            bool(frame.queue_block.task_depth)
            or max(0, _int(frame.sample_stats.get("running"))) > 0
            or self._sampler_outstanding(frame)
        )

        if not core:
            core_item = ("CORE", "UNKNOWN", "board missing", "hb —")
        else:
            fact, compact = _age_fact(core, now)
            redis_status = str(frame.proc_redis.get("status") or "").lower()
            if redis_status in {"failed", "dead", "unreachable", "error"}:
                core_item = ("CORE", "CRITICAL", "redis unreachable", "redis")
            elif _float(core.get("ts")) is None:
                core_item = ("CORE", "UNKNOWN", "heartbeat unknown", "hb —")
            elif not frame.proc_redis:
                core_item = ("CORE", "UNKNOWN", "redis board missing", "redis —")
            elif redis_status not in {"running", "ready", "healthy", "ok"}:
                core_item = ("CORE", "DEGRADED", "redis state unknown", "redis ?")
            else:
                core_item = ("CORE", "HEALTHY", fact, compact)

        if not core:
            factory_item = ("FACTORY", "UNKNOWN", "board missing", "—")
        elif mode == "paused":
            factory_item = ("FACTORY", "CRITICAL", "watchdog paused", "paused")
        elif mode == "degraded":
            factory_item = ("FACTORY", "DEGRADED", "watchdog degraded", "degraded")
        elif mode not in {"running", "stopping", "draining"}:
            factory_item = ("FACTORY", "UNKNOWN", "mode unknown", "—")
        elif mode in {"stopping", "draining"}:
            factory_item = ("FACTORY", "DRAINING", "draining", "drain")
        elif not outstanding and frame.queue_block.archive_depth:
            factory_item = ("FACTORY", "DRAINING", "draining", "drain")
        elif not outstanding and self._sampler_complete(frame):
            factory_item = ("FACTORY", "COMPLETE", "complete", "done")
        else:
            factory_item = ("FACTORY", "HEALTHY", "dispatching", "run")

        worker = frame.worker_block
        expected, present = worker.expected, worker.present
        if mode == "paused":
            worker_item = ("WORKERS", "PAUSED", "watchdog", "paused")
        elif expected is None or present is None:
            worker_item = ("WORKERS", "UNKNOWN", "liveness unknown", "—")
        elif outstanding and expected > 0 and present <= 0:
            worker_item = ("WORKERS", "CRITICAL", "0 workers alive", f"0/{expected}")
        elif worker.heartbeat_stale or worker.process_unknown or present < expected:
            count = max(worker.heartbeat_stale, expected - present)
            worker_item = ("WORKERS", "DEGRADED", f"{count} stale", f"{present}/{expected}")
        else:
            worker_item = (
                "WORKERS",
                "HEALTHY",
                f"{present}/{expected} alive",
                f"{present}/{expected}",
            )

        archive = frame.proc_archiver
        archive_depth = frame.queue_block.archive_depth
        archive_status = str(archive.get("status") or "").lower()
        if not archive:
            state = "CRITICAL" if archive_depth else "UNKNOWN"
            fact = "archiver missing" if archive_depth else "board missing"
            archiver_item = ("ARCHIVER", state, fact, "arch —")
        elif archive_depth and archive_status in {"failed", "dead", "stopped", "error"}:
            archiver_item = ("ARCHIVER", "CRITICAL", "archiver dead", "dead")
        else:
            depth = max(0, _int(archive_depth))
            fact = "queue empty" if depth == 0 else f"queue {depth}"
            archiver_item = ("ARCHIVER", "HEALTHY", fact, f"q {depth}")

        return core_item, factory_item, worker_item, archiver_item

    @staticmethod
    def _sampler_outstanding(frame: CollectorFrame) -> bool:
        progress = frame.sampler_display.progress
        current = _float(progress.current)
        target = _float(progress.target)
        return current is not None and target is not None and current < target

    @staticmethod
    def _sampler_complete(frame: CollectorFrame) -> bool:
        progress = frame.sampler_display.progress
        current = _float(progress.current)
        target = _float(progress.target)
        return current is not None and target is not None and current >= target


__all__ = ["LiveOverviewProjector"]
