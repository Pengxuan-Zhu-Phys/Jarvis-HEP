#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import statistics
import textwrap
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _utc_iso_from_epoch(epoch: float) -> str:
    return datetime.fromtimestamp(float(epoch), tz=timezone.utc).isoformat()


def _coerce_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except Exception:
        return None


def _coerce_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except Exception:
        return None


def _safe_div(numerator: float | int | None, denominator: float | int | None) -> float | None:
    numerator_f = _coerce_float(numerator)
    denominator_f = _coerce_float(denominator)
    if numerator_f is None or denominator_f is None or denominator_f <= 0:
        return None
    return numerator_f / denominator_f


def _mean(values: list[float]) -> float | None:
    if not values:
        return None
    return float(statistics.fmean(values))


def _median(values: list[float]) -> float | None:
    if not values:
        return None
    return float(statistics.median(values))


def _wrap_cell(text: str, width: int) -> list[str]:
    wrapped = textwrap.wrap(
        text or "",
        width=max(1, int(width)),
        break_long_words=True,
        break_on_hyphens=False,
        drop_whitespace=False,
    )
    return wrapped or [""]


class _NullResourceSampler:
    def start(self) -> None:
        return None

    def stop(self) -> None:
        return None

    def summary(self) -> dict[str, Any]:
        return {
            "cpu_percent_mean": None,
            "memory_rss_mb_peak": None,
            "open_file_peak": None,
        }


class _ProcessTreeResourceSampler:
    def __init__(self, poll_interval_sec: float = 1.0) -> None:
        self.poll_interval_sec = max(0.2, float(poll_interval_sec))
        self._lock = threading.Lock()
        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None
        self._running = False
        self._cpu_sum = 0.0
        self._cpu_samples = 0
        self._memory_rss_mb_peak = 0.0
        self._memory_seen = False
        self._open_file_peak = 0
        self._open_file_seen = False
        self._available = False

        try:
            import psutil  # type: ignore

            self._psutil = psutil
            self._root_proc = psutil.Process(os.getpid())
            self._available = True
        except Exception:
            self._psutil = None
            self._root_proc = None

    def start(self) -> None:
        if not self._available or self._running:
            return
        self._prime_cpu_counters()
        self._stop_event.clear()
        self._thread = threading.Thread(
            target=self._run,
            name="JarvisRunSummaryResourceSampler",
            daemon=True,
        )
        self._thread.start()
        self._running = True

    def stop(self) -> None:
        if not self._available:
            return
        self._sample_once()
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.poll_interval_sec + 1.0))
            self._thread = None
        self._running = False

    def summary(self) -> dict[str, Any]:
        with self._lock:
            cpu_mean = None
            if self._cpu_samples > 0:
                cpu_mean = self._cpu_sum / float(self._cpu_samples)
            return {
                "cpu_percent_mean": cpu_mean,
                "memory_rss_mb_peak": self._memory_rss_mb_peak if self._memory_seen else None,
                "open_file_peak": self._open_file_peak if self._open_file_seen else None,
            }

    def _run(self) -> None:
        while not self._stop_event.wait(self.poll_interval_sec):
            self._sample_once()

    def _iter_processes(self):
        if not self._available or self._root_proc is None:
            return []
        try:
            procs = [self._root_proc]
            procs.extend(self._root_proc.children(recursive=True))
            return procs
        except Exception:
            return []

    def _prime_cpu_counters(self) -> None:
        for proc in self._iter_processes():
            try:
                proc.cpu_percent(None)
            except Exception:
                continue

    def _sample_once(self) -> None:
        if not self._available:
            return

        total_cpu = 0.0
        cpu_seen = False
        total_rss = 0.0
        rss_seen = False
        total_open_files = 0
        open_files_seen = False

        for proc in self._iter_processes():
            try:
                total_cpu += float(proc.cpu_percent(None))
                cpu_seen = True
            except Exception:
                pass

            try:
                total_rss += float(proc.memory_info().rss)
                rss_seen = True
            except Exception:
                pass

            try:
                total_open_files += len(proc.open_files())
                open_files_seen = True
            except Exception:
                pass

        with self._lock:
            if cpu_seen:
                self._cpu_sum += total_cpu
                self._cpu_samples += 1
            if rss_seen:
                self._memory_seen = True
                self._memory_rss_mb_peak = max(
                    self._memory_rss_mb_peak,
                    total_rss / (1024.0 * 1024.0),
                )
            if open_files_seen:
                self._open_file_seen = True
                self._open_file_peak = max(self._open_file_peak, int(total_open_files))


class RunSummaryCollector:
    def __init__(
        self,
        *,
        output_dir: str,
        project_name: str | None,
        sampler_name: str | None,
        configured_workers: int | None = None,
        run_label: str | None = None,
        run_id: str | None = None,
        resource_sampler: Any | None = None,
        resource_poll_interval_sec: float = 1.0,
    ) -> None:
        self.output_dir = os.path.abspath(str(output_dir))
        self.project_name = project_name
        self.sampler_name = sampler_name
        self.configured_workers = _coerce_int(configured_workers)
        self.run_label = str(run_label or project_name or sampler_name or "run")
        self.run_id = str(run_id or self._build_run_id(self.run_label))

        if resource_sampler is None:
            resource_sampler = _ProcessTreeResourceSampler(
                poll_interval_sec=resource_poll_interval_sec
            )
        self.resource_sampler = resource_sampler or _NullResourceSampler()

        self._lock = threading.Lock()
        self._factory = None
        self._scheduler = None
        self._started = False
        self._finished = False
        self._start_epoch: float | None = None
        self._end_epoch: float | None = None
        self._start_time_iso: str | None = None
        self._end_time_iso: str | None = None
        self._resource_summary: dict[str, Any] | None = None
        self._external_tool_time_sec = 0.0
        self._crashed_subtasks = 0
        self._skipped_subtasks: int | None = None

    @staticmethod
    def _build_run_id(run_label: str) -> str:
        ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        pid = os.getpid()
        return f"{run_label}-{ts}-pid{pid}"

    def attach_factory(self, factory) -> None:
        self._factory = factory

    def attach_scheduler(self, scheduler) -> None:
        self._scheduler = scheduler

    def start(self) -> None:
        if self._started:
            return
        now = time.time()
        self._start_epoch = now
        self._start_time_iso = _utc_iso_from_epoch(now)
        self._started = True
        try:
            self.resource_sampler.start()
        except Exception:
            self.resource_sampler = _NullResourceSampler()

    def finish(self) -> dict[str, Any]:
        if not self._started:
            self.start()
        if not self._finished:
            now = time.time()
            self._end_epoch = now
            self._end_time_iso = _utc_iso_from_epoch(now)
            try:
                self.resource_sampler.stop()
            finally:
                try:
                    self._resource_summary = dict(self.resource_sampler.summary() or {})
                except Exception:
                    self._resource_summary = _NullResourceSampler().summary()
            self._finished = True
        return self.build_summary()

    def record_external_command(
        self,
        *,
        duration_sec: float | None,
        ok: bool,
        timed_out: bool = False,
    ) -> None:
        duration = _coerce_float(duration_sec)
        with self._lock:
            if duration is not None and duration >= 0.0:
                self._external_tool_time_sec += duration
            if (not ok) or bool(timed_out):
                self._crashed_subtasks += 1

    def record_skipped_subtasks(self, count: int) -> None:
        count_i = _coerce_int(count)
        if count_i is None:
            return
        with self._lock:
            if self._skipped_subtasks is None:
                self._skipped_subtasks = 0
            self._skipped_subtasks += max(0, count_i)

    def build_summary(self) -> dict[str, Any]:
        if self._start_epoch is None:
            self.start()
        end_epoch = self._end_epoch if self._end_epoch is not None else time.time()
        start_epoch = self._start_epoch if self._start_epoch is not None else end_epoch
        wall_time_sec = max(0.0, float(end_epoch - start_epoch))

        factory_metrics = {}
        if self._factory is not None and hasattr(self._factory, "get_run_metrics"):
            try:
                factory_metrics = dict(self._factory.get_run_metrics() or {})
            except Exception:
                factory_metrics = {}

        resource_summary = dict(self._resource_summary or {})
        if not resource_summary:
            try:
                resource_summary = dict(self.resource_sampler.summary() or {})
            except Exception:
                resource_summary = _NullResourceSampler().summary()

        with self._lock:
            external_tool_time_sec = float(self._external_tool_time_sec)
            crashed_subtasks = int(self._crashed_subtasks)
            skipped_subtasks = self._skipped_subtasks

        total_points_submitted = _coerce_int(factory_metrics.get("submitted")) or 0
        total_points_finished = _coerce_int(factory_metrics.get("ok")) or 0
        total_points_failed = _coerce_int(factory_metrics.get("failed")) or 0
        completed_durations_sec = list(factory_metrics.get("completed_durations_sec") or [])
        total_point_eval_sec = _coerce_float(factory_metrics.get("total_point_eval_sec")) or 0.0

        retry_count = factory_metrics.get("retry_count")
        if retry_count is not None:
            retry_count = _coerce_int(retry_count)

        time_in_framework_sec = max(0.0, float(total_point_eval_sec - external_tool_time_sec))
        framework_overhead_fraction = _safe_div(
            time_in_framework_sec,
            time_in_framework_sec + external_tool_time_sec,
        )

        summary = {
            "run_id": self.run_id,
            "project_name": self.project_name,
            "sampler_name": self.sampler_name,
            "start_time": self._start_time_iso or _utc_iso_from_epoch(start_epoch),
            "end_time": self._end_time_iso or _utc_iso_from_epoch(end_epoch),
            "wall_time_sec": wall_time_sec,
            "total_points_submitted": total_points_submitted,
            "total_points_finished": total_points_finished,
            "total_points_failed": total_points_failed,
            "success_rate": _safe_div(total_points_finished, total_points_submitted),
            "throughput_points_per_min": _safe_div(total_points_finished * 60.0, wall_time_sec),
            "configured_workers": self.configured_workers
            if self.configured_workers is not None
            else _coerce_int(factory_metrics.get("configured_workers")),
            "peak_active_workers": _coerce_int(factory_metrics.get("peak_active_workers")),
            "mean_active_workers": _coerce_float(factory_metrics.get("mean_active_workers")),
            "avg_point_eval_sec": _mean(completed_durations_sec),
            "median_point_eval_sec": _median(completed_durations_sec),
            "time_in_external_tools_sec": external_tool_time_sec,
            "time_in_framework_sec": time_in_framework_sec,
            "framework_overhead_fraction": framework_overhead_fraction,
            "retry_count": retry_count,
            "crashed_subtasks": crashed_subtasks,
            "skipped_subtasks": skipped_subtasks,
            "cpu_percent_mean": _coerce_float(resource_summary.get("cpu_percent_mean")),
            "memory_rss_mb_peak": _coerce_float(resource_summary.get("memory_rss_mb_peak")),
            "open_file_peak": _coerce_int(resource_summary.get("open_file_peak")),
        }
        return summary


class RunSummaryRenderer:
    _EXPLANATORY_LINE = (
        "Tracked cumulative timing across completed points; not equal to wall-clock time."
    )

    _SUMMARY_FIELD_ORDER = [
        "run_id",
        "project_name",
        "sampler_name",
        "start_time",
        "end_time",
        "wall_time_sec",
        "total_points_submitted",
        "total_points_finished",
        "total_points_failed",
        "success_rate",
        "throughput_points_per_min",
        "configured_workers",
        "peak_active_workers",
        "mean_active_workers",
        "avg_point_eval_sec",
        "median_point_eval_sec",
        "time_in_external_tools_sec",
        "time_in_framework_sec",
        "framework_overhead_fraction",
        "retry_count",
        "crashed_subtasks",
        "skipped_subtasks",
        "cpu_percent_mean",
        "memory_rss_mb_peak",
        "open_file_peak",
    ]

    _SECTIONS = (
        (
            "Run Overview",
            None,
            [
                ("run_id", "Run ID"),
                ("project_name", "Project"),
                ("sampler_name", "Sampler"),
                ("start_time", "Start Time"),
                ("end_time", "End Time"),
                ("wall_time_sec", "Wall Time (s)"),
                ("configured_workers", "Workers Configured"),
                ("peak_active_workers", "Peak Active Workers"),
                ("mean_active_workers", "Mean Active Workers"),
                ("total_points_submitted", "Points Submitted"),
                ("total_points_finished", "Points Completed"),
                ("total_points_failed", "Points Failed"),
                ("success_rate", "Success Rate"),
                ("throughput_points_per_min", "Throughput (pts/min)"),
            ],
        ),
        (
            "Execution Breakdown",
            _EXPLANATORY_LINE,
            [
                ("avg_point_eval_sec", "Mean Point Time (s)"),
                ("median_point_eval_sec", "Median Point Time (s)"),
                ("time_in_external_tools_sec", "External Tool Time, Total (s)"),
                ("time_in_framework_sec", "Internal Workflow Time, Total (s)"),
                ("framework_overhead_fraction", "Internal Time Fraction"),
                ("retry_count", "Retry Count"),
                ("crashed_subtasks", "Crashed Tasks"),
                ("skipped_subtasks", "Skipped Tasks"),
                ("cpu_percent_mean", "Mean CPU Usage (%)"),
                ("memory_rss_mb_peak", "Peak RSS (MB)"),
                ("open_file_peak", "Peak Open Files"),
            ],
        ),
    )

    def render(self, summary: dict[str, Any]) -> str:
        blocks = []
        for title, note, rowspec in self._SECTIONS:
            rows = [
                (label, *self._format_value(key, summary.get(key)))
                for key, label in rowspec
            ]
            blocks.append(self._render_table(title, rows, note=note))
        return "\n\n".join(blocks) + "\n"

    def write_outputs(
        self,
        summary: dict[str, Any],
        output_dir: str,
        *,
        rendered_text: str | None = None,
    ) -> dict[str, str]:
        output_root = Path(output_dir).expanduser().resolve()
        output_root.mkdir(parents=True, exist_ok=True)

        json_path = output_root / "run_summary.json"
        csv_path = output_root / "run_summary.csv"
        txt_path = output_root / "run_summary.txt"

        rendered_text = rendered_text if rendered_text is not None else self.render(summary)

        ordered_summary = {
            key: summary.get(key)
            for key in self._SUMMARY_FIELD_ORDER
        }

        with open(json_path, "w", encoding="utf-8") as f1:
            json.dump(ordered_summary, f1, indent=2, ensure_ascii=False)
            f1.write("\n")

        with open(csv_path, "w", encoding="utf-8", newline="") as f1:
            writer = csv.DictWriter(f1, fieldnames=self._SUMMARY_FIELD_ORDER)
            writer.writeheader()
            writer.writerow(ordered_summary)

        with open(txt_path, "w", encoding="utf-8") as f1:
            f1.write(rendered_text)

        return {
            "json": str(json_path),
            "csv": str(csv_path),
            "txt": str(txt_path),
        }

    def _render_table(
        self,
        title: str,
        rows: list[tuple[str, str, bool]],
        *,
        note: str | None = None,
        max_width: int = 78,
        min_value_width: int = 11,
    ) -> str:
        metric_header = "Metric"
        value_header = "Value"
        metric_width = max(len(metric_header), *(len(metric) for metric, _, _ in rows))
        metric_width = min(metric_width, max_width - 7 - max(1, int(min_value_width)))
        value_width = max(int(min_value_width), max_width - metric_width - 7)
        longest_value = max(len(value_header), *(len(value) for _, value, _ in rows))
        value_width = max(value_width, longest_value)

        border = f"+{'-' * (metric_width + 2)}+{'-' * (value_width + 2)}+"
        lines = [f"[{title}]"]
        if note:
            lines.append(str(note))
        lines.append(border)
        lines.append(
            f"| {metric_header:<{metric_width}} | {value_header:<{value_width}} |"
        )
        lines.append(border)

        for metric, value, value_is_numeric in rows:
            metric_chunks = _wrap_cell(metric, metric_width)
            value_chunks = _wrap_cell(value, value_width)
            height = max(len(metric_chunks), len(value_chunks))
            for idx in range(height):
                metric_chunk = metric_chunks[idx] if idx < len(metric_chunks) else ""
                value_chunk = value_chunks[idx] if idx < len(value_chunks) else ""
                lines.append(
                    f"| {metric_chunk:<{metric_width}} | {value_chunk:>{value_width}} |"
                )

        lines.append(border)
        return "\n".join(lines)

    def _format_value(self, key: str, value: Any) -> tuple[str, bool]:
        if value is None:
            return "N/A", False

        if key in {"success_rate", "framework_overhead_fraction"}:
            return f"{float(value) * 100.0:.2f}%", True

        if key == "cpu_percent_mean":
            return f"{float(value):.1f}%", True

        if key in {
            "wall_time_sec",
            "throughput_points_per_min",
            "mean_active_workers",
            "avg_point_eval_sec",
            "median_point_eval_sec",
            "time_in_external_tools_sec",
            "time_in_framework_sec",
            "memory_rss_mb_peak",
        }:
            return f"{float(value):.3f}", True

        if key in {
            "configured_workers",
            "peak_active_workers",
            "total_points_submitted",
            "total_points_finished",
            "total_points_failed",
            "retry_count",
            "crashed_subtasks",
            "skipped_subtasks",
            "open_file_peak",
        }:
            return f"{int(value)}", True

        if isinstance(value, float):
            return f"{value:.3f}", True

        return str(value), False
