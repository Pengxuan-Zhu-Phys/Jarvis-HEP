#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import threading
import time
from typing import Any


DEFAULT_BENCHMARK_SECONDS = 30.0
TIMING_ENABLED = False
_STAGE_LOCK = threading.Lock()
_STAGE_TOTALS: dict[str, float] = {}
_STAGE_COUNTS: dict[str, int] = {}


def monotonic_seconds() -> float:
    return time.perf_counter()


def reset_stage_timers() -> None:
    with _STAGE_LOCK:
        _STAGE_TOTALS.clear()
        _STAGE_COUNTS.clear()


def enable_stage_timers(*, reset: bool = True) -> None:
    global TIMING_ENABLED
    if reset:
        reset_stage_timers()
    TIMING_ENABLED = True


def disable_stage_timers() -> None:
    global TIMING_ENABLED
    TIMING_ENABLED = False


def record_stage(stage: str, elapsed_seconds: float) -> None:
    if elapsed_seconds <= 0.0:
        return
    with _STAGE_LOCK:
        _STAGE_TOTALS[stage] = _STAGE_TOTALS.get(stage, 0.0) + float(elapsed_seconds)
        _STAGE_COUNTS[stage] = _STAGE_COUNTS.get(stage, 0) + 1


def snapshot_stage_timers(
    *,
    completed_samples: int = 0,
    wall_time_sec: float = 0.0,
) -> dict[str, Any]:
    with _STAGE_LOCK:
        totals = dict(_STAGE_TOTALS)
        counts = dict(_STAGE_COUNTS)

    completed = max(0, int(completed_samples or 0))
    wall_time = max(0.0, float(wall_time_sec or 0.0))
    stages = {}
    cumulative_total = 0.0
    cumulative_per_sample = 0.0
    for name in sorted(totals.keys()):
        total = float(totals[name])
        count = int(counts.get(name, 0))
        per_sample = (total / float(completed)) if completed > 0 else None
        cumulative_total += total
        if per_sample is not None:
            cumulative_per_sample += per_sample
        stages[name] = {
            "total_sec": total,
            "count": count,
            "avg_call_sec": (total / float(count)) if count > 0 else None,
            "per_sample_sec": per_sample,
            "wall_fraction": (total / wall_time) if wall_time > 0.0 else None,
        }

    return {
        "enabled": bool(stages),
        "completed_samples": completed,
        "wall_time_sec": wall_time,
        "stages": stages,
        "cumulative_total_sec": cumulative_total,
        "cumulative_per_sample_sec": cumulative_per_sample if completed > 0 else None,
        "cumulative_wall_fraction": (cumulative_total / wall_time) if wall_time > 0.0 else None,
    }


def normalize_benchmark_seconds(value: Any) -> float:
    if value is None:
        return DEFAULT_BENCHMARK_SECONDS
    try:
        seconds = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid benchmark duration: {value!r}") from exc
    if seconds <= 0:
        raise ValueError("--benchmark duration must be greater than zero seconds")
    return seconds


class BenchmarkDeadline:
    """Stop sampler submission after a fixed wall-clock duration."""

    def __init__(self, sampler: Any, seconds: float) -> None:
        self.sampler = sampler
        self.seconds = normalize_benchmark_seconds(seconds)
        self.started_at_monotonic: float | None = None
        self.ended_at_monotonic: float | None = None
        self.deadline_monotonic: float | None = None
        self._original_next_sample = None

    def __enter__(self):
        self.started_at_monotonic = time.monotonic()
        self.deadline_monotonic = self.started_at_monotonic + self.seconds
        self._original_next_sample = getattr(self.sampler, "next_sample", None)
        if not callable(self._original_next_sample):
            raise RuntimeError(
                f"Benchmark mode requires sampler {self.sampler!r} to expose next_sample()."
            )
        enable_stage_timers(reset=True)

        def _timed_next_sample(*args, **kwargs):
            if time.monotonic() >= float(self.deadline_monotonic):
                raise StopIteration
            return self._original_next_sample(*args, **kwargs)

        self.sampler.next_sample = _timed_next_sample
        return self

    def __exit__(self, exc_type, exc, tb):
        self.ended_at_monotonic = time.monotonic()
        if self._original_next_sample is not None:
            self.sampler.next_sample = self._original_next_sample
        disable_stage_timers()
        return False

    @property
    def elapsed_seconds(self) -> float:
        if self.started_at_monotonic is None:
            return 0.0
        end = self.ended_at_monotonic if self.ended_at_monotonic is not None else time.monotonic()
        return max(0.0, float(end - self.started_at_monotonic))


def collect_benchmark_report(
    *,
    factory: Any,
    requested_seconds: float,
    wall_time_sec: float,
    run_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    metrics = {}
    if factory is not None and hasattr(factory, "get_run_metrics"):
        try:
            metrics = dict(factory.get_run_metrics() or {})
        except Exception:
            metrics = {}

    submitted = int(metrics.get("submitted", 0) or 0)
    completed = int(metrics.get("ok", 0) or 0)
    failed = int(metrics.get("failed", 0) or 0)
    wall_time = max(0.0, float(wall_time_sec))
    samples_per_sec = (float(completed) / wall_time) if wall_time > 0.0 else 0.0

    report = {
        "requested_seconds": float(requested_seconds),
        "wall_time_sec": wall_time,
        "samples_per_sec": samples_per_sec,
        "total_submitted": submitted,
        "total_completed": completed,
        "total_failed": failed,
        "stage_timers": snapshot_stage_timers(
            completed_samples=completed,
            wall_time_sec=wall_time,
        ),
    }
    if isinstance(run_summary, dict):
        report["run_summary"] = dict(run_summary)
        for key, value in run_summary.items():
            report.setdefault(key, value)
    return report


def write_benchmark_report(report: dict[str, Any], output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, "benchmark.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return path
