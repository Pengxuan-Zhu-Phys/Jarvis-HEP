"""Redis talker for the live TUI. The only process that writes hep:monitor:want."""

from __future__ import annotations

import os
import re
import time
from collections import deque
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field, replace
from typing import Any
from urllib.parse import urlsplit

from jarvishep2.dashboard import attach_reader
from jarvishep2.monitor.calculators import CalculatorPoolData, CalculatorsBlockData
from jarvishep2.monitor.scans import ScanChoice
from jarvishep2.monitor.metrics import format_sample_rate
from jarvishep2.monitor.queues import (
    QueueBlockData,
    QueueDirection,
    QueueDirectionWindow,
    feedback_queue_plan,
)
from jarvishep2.monitor.samples import SamplerDisplay, sampler_display_from_sources
from jarvishep2.monitor.workers import WorkerBlockData
from jarvishep2.redis_queue import MONITOR_WANT_TTL_SEC, RedisQueue
from jarvishep2.runtime_metadata import calculator_pools_from_metadata


def want_channels_for_page(slug: str) -> tuple[bool, bool] | None:
    """Return (ch:calc, ch:sample) or None when this page must not hold want."""
    if slug == "calculators":
        return True, False
    if slug == "samples":
        return False, True
    return None


@dataclass
class CollectorFrame:
    page: str
    occupancy: dict[str, dict[str, int]] = field(default_factory=dict)
    calc_busy: dict[str, str] = field(default_factory=dict)
    sample_running: dict[str, dict[str, Any]] = field(default_factory=dict)
    workers: list[dict[str, Any]] = field(default_factory=list)
    sample_stats: dict[str, Any] = field(default_factory=dict)
    queues: dict[str, Any] = field(default_factory=dict)
    queue_block: QueueBlockData = field(
        default_factory=lambda: QueueBlockData(None, None, feedback_mode="unknown")
    )
    worker_block: WorkerBlockData = field(default_factory=WorkerBlockData)
    factory: dict[str, Any] = field(default_factory=dict)
    sampler: dict[str, Any] = field(default_factory=dict)
    sampler_display: SamplerDisplay = field(default_factory=SamplerDisplay)
    sample_rate: str = "—"
    sample_rate_per_sec: float | None = None
    calculator_block: CalculatorsBlockData = field(default_factory=CalculatorsBlockData)
    op_counts: dict[str, int] = field(default_factory=dict)
    proc_core: dict[str, Any] = field(default_factory=dict)
    proc_archiver: dict[str, Any] = field(default_factory=dict)
    proc_redis: dict[str, Any] = field(default_factory=dict)
    host: dict[str, Any] = field(default_factory=dict)
    selected_calc: str | None = None
    want_on: bool = False
    stale: bool = False
    error: str | None = None


class Collector:
    """Want lifecycle + occupancy / sidecar reads. No Textual imports."""

    def __init__(
        self,
        redis: RedisQueue,
        *,
        owner_ids: list[str] | None = None,
        pid: int | None = None,
        pool_names: list[str] | None = None,
        slots: dict[str, int] | None = None,
        owns_client: bool = False,
        selected_calc: str | None = None,
        process_inventory: list[Mapping[str, Any]] | None = None,
        host_snapshotter: Callable[[], Mapping[str, Any]] | None = None,
        sampler_metadata: Mapping[str, Any] | None = None,
        monotonic: Callable[[], float] | None = None,
        pid_alive: Callable[[int], bool | None] | None = None,
    ) -> None:
        self._redis = redis
        self._owner_ids = [str(owner).strip() for owner in owner_ids or [] if str(owner).strip()]
        self._pid = int(os.getpid() if pid is None else pid)
        self._pool_names = [str(name).strip() for name in pool_names or [] if str(name).strip()]
        self._slots = dict(slots or {})
        self._owns_client = bool(owns_client)
        self._selected_calc = str(selected_calc).strip() if selected_calc else None
        if self._selected_calc not in self._pool_names:
            self._selected_calc = self._pool_names[0] if self._pool_names else None
        self._want_on = False
        self._process_inventory = [dict(row) for row in process_inventory or []]
        self._host_snapshotter = host_snapshotter or self._collect_host_snapshot
        self._sampler_metadata = dict(sampler_metadata or {})
        metadata_config = self._sampler_metadata.get("config")
        self._feedback_mode, self._feedback_chain_ids = feedback_queue_plan(
            str(self._sampler_metadata.get("method") or "Unknown"),
            metadata_config if isinstance(metadata_config, Mapping) else None,
        )
        self._monotonic = monotonic or time.monotonic
        self._pid_alive = pid_alive or self._pid_is_alive
        self._sample_rate_points: deque[tuple[float, int]] = deque()
        self._sampler_progress_key: tuple[str, str, float] | None = None
        self._sampler_progress_points: deque[tuple[float, float]] = deque()
        self._queue_directions = QueueDirectionWindow()
        self._host_cpu_sample_at: float | None = None
        self._host_cpu_percent: float = 0.0
        self._host_processes: dict[int, Any] = {}
        self.last_frame = CollectorFrame(page="overview")

    @property
    def redis_endpoint(self) -> str:
        """Return a credential-free endpoint label for STATUS."""
        config = dict(getattr(self._redis, "config", {}) or {})
        url = str(config.get("url") or "").strip()
        if url:
            parsed = urlsplit(url)
            host = parsed.hostname or "redis"
            port = parsed.port or 6379
            return f"{host}:{port}"
        host = str(config.get("host") or "localhost")
        try:
            port = int(config.get("port", 6379))
        except (TypeError, ValueError, OverflowError):
            port = 6379
        return f"{host}:{port}"

    def set_pools(
        self,
        names: list[str],
        slots: dict[str, int] | None = None,
    ) -> None:
        self._pool_names = [str(name).strip() for name in names if str(name).strip()]
        self._slots = dict(slots or {})
        if self._selected_calc not in self._pool_names:
            self._selected_calc = self._pool_names[0] if self._pool_names else None

    def select_calc(self, name: str | None) -> None:
        text = str(name or "").strip()
        self._selected_calc = text or None

    def cycle_calc(self, delta: int) -> str | None:
        """Move the Calculators-page selection without issuing a Redis command."""
        if not self._pool_names:
            self._selected_calc = None
            return None
        try:
            index = self._pool_names.index(self._selected_calc or "")
        except ValueError:
            index = 0
        self._selected_calc = self._pool_names[(index + int(delta)) % len(self._pool_names)]
        return self._selected_calc

    def tick(self, page: str) -> CollectorFrame:
        slug = str(page or "").strip() or "overview"
        try:
            self._sync_want(slug)
            occupancy: dict[str, dict[str, int]] = {}
            calc_busy: dict[str, str] = {}
            sample_running: dict[str, dict[str, Any]] = {}
            workers: list[dict[str, Any]] = []
            sample_stats: dict[str, Any] = {}
            queues: dict[str, Any] = {}
            factory: dict[str, Any] = {}
            sampler: dict[str, Any] = {}
            op_counts: dict[str, int] = {}
            proc_core: dict[str, Any] = {}
            proc_archiver: dict[str, Any] = {}
            proc_redis: dict[str, Any] = {}
            host: dict[str, Any] = {}
            queue_block = self.last_frame.queue_block
            worker_block = self.last_frame.worker_block
            sample_rate = self.last_frame.sample_rate
            sample_rate_per_sec = self.last_frame.sample_rate_per_sec
            sampler_display = self.last_frame.sampler_display
            calculator_block = self.last_frame.calculator_block
            if slug in {"overview", "workers", "factory", "sampler", "samples"}:
                view = attach_reader(
                    redis=self._redis,
                    owner_ids=self._owner_ids,
                    calculator_names=self._pool_names or None,
                    calculator_slots=self._slots or None,
                    feedback_chain_ids=list(self._feedback_chain_ids) or None,
                ).read()
                occupancy = view.occupancy
                workers = view.workers
                sample_stats = view.samples
                queues = view.queues
                factory = view.factory
                sampler = view.sampler
                op_counts = view.op_counts
                proc_core = view.proc_core
                proc_archiver = view.proc_archiver
                proc_redis = view.proc_redis
                observed_at = self._monotonic()
                sample_rate, sample_rate_per_sec = self._observe_sample_rate(
                    sample_stats,
                    observed_at=observed_at,
                )
                sampler_display = sampler_display_from_sources(
                    self._sampler_metadata,
                    proc_core,
                )
                sampler_display = self._observe_sampler_eta(
                    sampler_display,
                    observed_at=observed_at,
                )
                feedback_mode, feedback_chain_ids = feedback_queue_plan(
                    sampler_display.method,
                    sampler_display.metrics,
                )
                if feedback_mode != "unknown":
                    self._feedback_mode = feedback_mode
                    self._feedback_chain_ids = feedback_chain_ids
                queue_block = self._observe_queue_block(
                    queues,
                    observed_at=observed_at,
                    feedback_mode=self._feedback_mode,
                )
                worker_block = self._observe_worker_block(workers, proc_core)
                calculator_block = self._calculator_block(occupancy)
            if slug in {"overview", "host"}:
                host = dict(self._host_snapshotter())
            elif self._pool_names:
                occupancy = self._redis.fetch_calc_occupancy(
                    self._pool_names, slots=self._slots
                )
            if slug == "calculators" and self._selected_calc:
                calc_busy = self._redis.fetch_monitor_calc_busy(self._selected_calc)
            if slug == "samples":
                sample_running = self._redis.fetch_monitor_sample_running()
            frame = CollectorFrame(
                page=slug,
                occupancy=occupancy,
                calc_busy=calc_busy,
                sample_running=sample_running,
                workers=workers,
                sample_stats=sample_stats,
                queues=queues,
                queue_block=queue_block,
                worker_block=worker_block,
                factory=factory,
                sampler=sampler,
                sampler_display=sampler_display,
                sample_rate=sample_rate,
                sample_rate_per_sec=sample_rate_per_sec,
                calculator_block=calculator_block,
                op_counts=op_counts,
                proc_core=proc_core,
                proc_archiver=proc_archiver,
                proc_redis=proc_redis,
                host=host,
                selected_calc=self._selected_calc,
                want_on=self._want_on,
            )
        except Exception as exc:
            # Redis is observational here: retain the last good snapshot and
            # leave the next interval to refresh it.  A failed want write is
            # equally non-fatal; the missing/expired key is fail-closed.
            return replace(
                self.last_frame,
                page=slug,
                want_on=self._want_on,
                stale=True,
                queue_block=self._stale_queue_block(self.last_frame.queue_block),
                worker_block=replace(self.last_frame.worker_block, source_stale=True),
                error=str(exc) or exc.__class__.__name__,
            )
        self.last_frame = frame
        return frame

    def _observe_queue_block(
        self,
        queues: Mapping[str, Any],
        *,
        observed_at: float,
        feedback_mode: str,
    ) -> QueueBlockData:
        """Derive one bounded QueueBlockData from already-read queue lengths."""
        def depth(name: str) -> int | None:
            value = queues.get(name)
            if value is None:
                return None
            try:
                return max(0, int(value))
            except (TypeError, ValueError, OverflowError):
                return None

        task_depth = depth("task_queue_length")
        archive_depth = depth("archive_queue_length")
        mode = str(feedback_mode or "unknown").lower()
        raw_shards = queues.get("feedback_shards")
        shards = {
            str(key): max(0, int(value))
            for key, value in dict(raw_shards or {}).items()
            if value is not None
        }
        feedback_depth: int | None = None
        if mode == "shared":
            feedback_depth = depth("feedback_queue_length")
        elif mode == "sharded" and shards:
            feedback_depth = sum(shards.values())
        return QueueBlockData(
            task_depth=task_depth,
            archive_depth=archive_depth,
            feedback_depth=feedback_depth,
            feedback_mode=mode,
            feedback_shards=shards if mode == "sharded" else {},
            task_direction=self._queue_directions.observe(
                "task", task_depth, observed_at=observed_at
            ),
            archive_direction=self._queue_directions.observe(
                "archive", archive_depth, observed_at=observed_at
            ),
            feedback_direction=self._queue_directions.observe(
                "feedback",
                None if mode in {"inactive", "unknown"} else feedback_depth,
                observed_at=observed_at,
            ),
        )

    @staticmethod
    def _pid_is_alive(pid: int) -> bool | None:
        """Check one admitted PID; never search for children or process names."""
        try:
            import psutil
        except ImportError:
            return None
        try:
            return bool(psutil.Process(int(pid)).is_running())
        except psutil.NoSuchProcess:
            return False
        except (psutil.AccessDenied, TypeError, ValueError, OverflowError):
            return None

    def _observe_worker_block(
        self,
        workers: list[dict[str, Any]],
        proc_core: Mapping[str, Any],
    ) -> WorkerBlockData:
        """Build fleet evidence from approved Worker PIDs and their boards only."""
        def positive(value: Any) -> int | None:
            try:
                number = int(value)
            except (TypeError, ValueError, OverflowError):
                return None
            return number if number > 0 else None

        core_expected = positive(proc_core.get("workers_total"))
        inventory_pids = {
            pid
            for row in self._process_inventory
            if str(row.get("role") or "") == "worker"
            for pid in [positive(row.get("pid"))]
            if pid is not None
        }
        board_pids = {
            pid
            for row in workers
            for pid in [positive(row.get("pid"))]
            if pid is not None
        }
        worker_pids = inventory_pids | board_pids
        present = 0
        process_unknown = 0
        for pid in worker_pids:
            alive = self._pid_alive(pid)
            if alive is True:
                present += 1
            elif alive is None:
                process_unknown += 1

        busy = sum(str(row.get("status") or "").lower() == "busy" for row in workers)
        idle = sum(str(row.get("status") or "").lower() == "idle" for row in workers)
        assigned = sum(bool(str(row.get("current_uuid") or "").strip()) for row in workers)

        file_ops_expected = file_ops_alive = file_ops_missing = file_ops_unknown = file_ops_inline = 0
        heartbeat_recent = heartbeat_stale = heartbeat_unknown = 0
        for row in workers:
            mode = str(row.get("file_operation_mode") or "").strip().lower()
            fo_pid = positive(row.get("file_operation_pid"))
            if mode == "inline":
                file_ops_inline += 1
            elif mode == "process" or fo_pid is not None:
                file_ops_expected += 1
                if fo_pid is None:
                    file_ops_unknown += 1
                else:
                    alive = self._pid_alive(fo_pid)
                    if alive is True:
                        file_ops_alive += 1
                    elif alive is False:
                        file_ops_missing += 1
                    else:
                        file_ops_unknown += 1

            age = row.get("heartbeat_age_s")
            ttl = positive(row.get("board_ttl_sec"))
            try:
                age_number = max(0.0, float(age))
            except (TypeError, ValueError, OverflowError):
                age_number = None
            if age_number is None or ttl is None:
                heartbeat_unknown += 1
            elif age_number <= ttl:
                heartbeat_recent += 1
            else:
                heartbeat_stale += 1

        expected = core_expected if core_expected is not None else len(worker_pids)
        return WorkerBlockData(
            expected=expected,
            present=present if worker_pids else (0 if expected == 0 else None),
            process_unknown=process_unknown,
            busy=busy,
            idle=idle,
            assigned=assigned,
            file_ops_expected=file_ops_expected,
            file_ops_alive=file_ops_alive,
            file_ops_missing=file_ops_missing,
            file_ops_unknown=file_ops_unknown,
            file_ops_inline=file_ops_inline,
            heartbeat_recent=heartbeat_recent,
            heartbeat_stale=heartbeat_stale,
            heartbeat_unknown=heartbeat_unknown,
        )

    def _calculator_block(
        self,
        occupancy: Mapping[str, Mapping[str, Any]],
    ) -> CalculatorsBlockData:
        """Project real Redis PackID occupancy into the Overview constellation."""
        names = [str(name) for name in occupancy]
        busy_by_pool = self._redis.fetch_calc_busy_pack_ids(names)
        pools: list[CalculatorPoolData] = []
        for name in names:
            row = occupancy.get(name) or {}
            try:
                slots = max(0, int(row.get("slots", 0) or 0))
            except (TypeError, ValueError, OverflowError):
                slots = 0
            pools.append(
                CalculatorPoolData(
                    name=name,
                    slots=slots,
                    busy_packs=busy_by_pool.get(name, ()),
                )
            )
        return CalculatorsBlockData(tuple(pools))

    @staticmethod
    def _stale_queue_block(data: QueueBlockData) -> QueueBlockData:
        """Keep the last valid depths while making a failed observation explicit."""
        feedback = (
            data.feedback_direction
            if str(data.feedback_mode).lower() == "inactive"
            else QueueDirection("stale")
        )
        return replace(
            data,
            task_direction=QueueDirection("stale"),
            archive_direction=QueueDirection("stale"),
            feedback_direction=feedback,
        )

    def _observe_sample_rate(
        self,
        stats: Mapping[str, Any],
        *,
        observed_at: float,
    ) -> tuple[str, float | None]:
        """Track a local 60-second completed-sample window."""
        try:
            completed = max(0, int(stats.get("completed", 0) or 0))
        except (TypeError, ValueError):
            return "—", None
        now = float(observed_at)
        if self._sample_rate_points and completed < self._sample_rate_points[-1][1]:
            self._sample_rate_points.clear()
        self._sample_rate_points.append((now, completed))
        while (
            len(self._sample_rate_points) > 2
            and now - self._sample_rate_points[1][0] >= 60.0
        ):
            self._sample_rate_points.popleft()
        if len(self._sample_rate_points) < 2:
            return "—", None
        oldest_at, oldest_completed = self._sample_rate_points[0]
        elapsed = now - oldest_at
        if elapsed <= 0:
            return "—", None
        rate_per_sec = max(0.0, (completed - oldest_completed) / elapsed)
        return format_sample_rate(rate_per_sec), rate_per_sec

    def _observe_sampler_eta(
        self,
        display: SamplerDisplay,
        *,
        observed_at: float,
    ) -> SamplerDisplay:
        """Estimate ETA only from the sampler's own monotonic progress unit."""
        progress = display.progress
        if progress.kind not in {"finite", "generation", "iteration"}:
            self._sampler_progress_key = None
            self._sampler_progress_points.clear()
            return display
        try:
            current = max(0.0, float(progress.current))
            target = float(progress.target)
        except (TypeError, ValueError, OverflowError):
            self._sampler_progress_key = None
            self._sampler_progress_points.clear()
            return display
        if target <= 0:
            return display
        key = (display.method, progress.kind, target)
        if key != self._sampler_progress_key:
            self._sampler_progress_key = key
            self._sampler_progress_points.clear()
        if self._sampler_progress_points and current < self._sampler_progress_points[-1][1]:
            self._sampler_progress_points.clear()
        now = float(observed_at)
        self._sampler_progress_points.append((now, current))
        while (
            len(self._sampler_progress_points) > 2
            and now - self._sampler_progress_points[1][0] >= 60.0
        ):
            self._sampler_progress_points.popleft()
        if len(self._sampler_progress_points) < 2 or current >= target:
            return display
        oldest_at, oldest_current = self._sampler_progress_points[0]
        elapsed = now - oldest_at
        advanced = current - oldest_current
        if elapsed <= 0 or advanced <= 0:
            return display
        eta_sec = int(round((target - current) / (advanced / elapsed)))
        hours, remainder = divmod(max(0, eta_sec), 3600)
        minutes, seconds = divmod(remainder, 60)
        eta = f"{hours}:{minutes:02d}:{seconds:02d}"
        return replace(display, progress=replace(progress, eta=eta))

    def close(self) -> None:
        try:
            if self._want_on:
                self._redis.clear_monitor_want()
        except Exception:
            pass
        self._want_on = False
        if self._owns_client:
            try:
                self._redis.close()
            except Exception:
                pass

    def _sync_want(self, page: str) -> None:
        channels = want_channels_for_page(page)
        if channels is None:
            if self._want_on:
                self._redis.clear_monitor_want()
                self._want_on = False
            return
        ch_calc, ch_sample = channels
        self._redis.set_monitor_want(
            ch_calc=ch_calc,
            ch_sample=ch_sample,
            ttl_sec=MONITOR_WANT_TTL_SEC,
            pid=self._pid,
        )
        self._want_on = True

    def _collect_host_snapshot(self) -> dict[str, Any]:
        """Collect host totals and only the PIDs already admitted by the scan list."""
        try:
            import psutil
        except ImportError:
            return {"available": False, "reason": "psutil is not installed", "processes": []}

        memory = psutil.virtual_memory()
        swap = psutil.swap_memory()
        try:
            load = tuple(float(value) for value in os.getloadavg())
        except (AttributeError, OSError):
            load = ()
        processes: list[dict[str, Any]] = []
        for known in self._process_inventory:
            try:
                pid = int(known.get("pid") or 0)
            except (TypeError, ValueError):
                continue
            if pid <= 0:
                continue
            row = {"role": str(known.get("role") or "process"), "pid": pid}
            try:
                process = self._host_processes.get(pid)
                if process is None or not process.is_running():
                    process = psutil.Process(pid)
                    self._host_processes[pid] = process
                with process.oneshot():
                    rss = int(process.memory_info().rss)
                    row.update(
                        {
                            "alive": process.is_running(),
                            "cpu_percent": float(process.cpu_percent(interval=None)),
                            "rss": rss,
                            "ppid": process.ppid(),
                            "threads": process.num_threads(),
                            "cmdline": " ".join(process.cmdline()),
                            "created_at": float(process.create_time()),
                        }
                    )
                    try:
                        row["fds"] = int(process.num_fds())
                    except (AttributeError, OSError):
                        row["fds"] = None
                    try:
                        soft_limit, _hard_limit = process.rlimit(psutil.RLIMIT_NOFILE)
                        row["fd_limit"] = int(soft_limit)
                    except (AttributeError, OSError, ValueError):
                        row["fd_limit"] = None
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                self._host_processes.pop(pid, None)
                row["alive"] = False
            processes.append(row)
        # psutil measures CPU since its previous call. A manual refresh next
        # to a timer tick must not replace a valid reading with a near-zero
        # duration sample (psutil recommends at least 0.1 seconds).
        now = time.monotonic()
        if self._host_cpu_sample_at is None or now - self._host_cpu_sample_at >= 0.1:
            self._host_cpu_percent = float(psutil.cpu_percent(interval=None))
            self._host_cpu_sample_at = now
        return {
            "available": True,
            "timestamp": time.time(),
            "cpu_percent": self._host_cpu_percent,
            "memory_used": int(memory.used),
            "memory_total": int(memory.total),
            "swap_used": int(swap.used),
            "swap_total": int(swap.total),
            "load": load,
            "processes": processes,
        }


def open_collector(choice: ScanChoice) -> Collector | None:
    """Connect to the live scan's Redis. Simulated choices return None."""
    if choice.simulated:
        return None
    from jarvishep2.process_cleanup import (
        list_active_scans,
        resolve_scan_reference,
        runtime_metadata_for_scan,
    )

    try:
        scan = resolve_scan_reference(choice.reference, list_active_scans())
    except ValueError:
        return None
    metadata = runtime_metadata_for_scan(scan)
    if not metadata:
        return None
    redis_config = dict(metadata.get("redis") or {})
    if not redis_config:
        return None
    pool_names, slots = calculator_pools_from_metadata(metadata)
    redis = RedisQueue(redis_config)
    try:
        redis.connect()
    except Exception:
        return None
    return Collector(
        redis,
        owner_ids=_worker_owner_ids(scan),
        owns_client=True,
        pool_names=pool_names or [],
        slots=slots or {},
        process_inventory=_process_inventory(scan, control_pid=choice.control_pid),
        sampler_metadata=(
            metadata.get("sampler")
            if isinstance(metadata.get("sampler"), Mapping)
            else None
        ),
    )


def _worker_owner_ids(scan: Any) -> list[str]:
    """Project Worker ids from the already-approved OS process inventory."""
    ids: list[str] = []
    for process in getattr(scan, "processes", ()) or ():
        title = str(getattr(process, "command", "")).strip().split(None, 1)[0]
        matched = re.match(r"^Jarvis-Worker-(\d+)(?::|$)", title)
        if matched is None:
            continue
        owner = str(int(matched.group(1)))
        if owner not in ids:
            ids.append(owner)
    return ids


def _process_inventory(scan: Any, *, control_pid: int | None = None) -> list[dict[str, Any]]:
    """Classify only processes returned by the existing scan inventory."""
    rows: list[dict[str, Any]] = []
    for process in getattr(scan, "processes", ()) or ():
        try:
            pid = int(getattr(process, "pid"))
        except (TypeError, ValueError):
            continue
        command = str(getattr(process, "command", "")).strip()
        title = command.split(None, 1)[0] if command else ""
        if re.match(r"^Jarvis-Worker-\d+(?::|$)", title):
            role = "worker"
        elif "archiver" in title.lower():
            role = "archiver"
        elif title.startswith("redis-server"):
            role = "redis"
        elif control_pid is not None and pid == int(control_pid):
            role = "core"
        else:
            role = "process"
        rows.append({"role": role, "pid": pid, "command": command})
    return rows


__all__ = [
    "Collector",
    "CollectorFrame",
    "open_collector",
    "want_channels_for_page",
]
