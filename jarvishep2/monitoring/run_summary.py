#!/usr/bin/env python3
"""Run summary projection for Jarvis-HEP V2 (WP-D5.2, frozen contract)."""

from __future__ import annotations

import csv
import json
import statistics
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


RUN_SUMMARY_FIELD_ORDER: tuple[str, ...] = (
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
)


def _utc_iso_from_epoch(epoch: float) -> str:
    return datetime.fromtimestamp(float(epoch), tz=timezone.utc).isoformat()


def _coerce_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _coerce_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
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


def validate_run_summary(summary: Mapping[str, Any]) -> None:
    """Raise ValueError when *summary* is missing a frozen-contract field."""
    missing = [key for key in RUN_SUMMARY_FIELD_ORDER if key not in summary]
    if missing:
        raise ValueError(f"run_summary missing fields: {', '.join(missing)}")


def build_run_summary(
    *,
    factory_metrics: Mapping[str, Any] | None = None,
    project_name: str | None = None,
    sampler_name: str | None = None,
    run_id: str | None = None,
    start_epoch: float | None = None,
    end_epoch: float | None = None,
    configured_workers: int | None = None,
    resource_summary: Mapping[str, Any] | None = None,
    external_tool_time_sec: float = 0.0,
    crashed_subtasks: int = 0,
    skipped_subtasks: int | None = None,
) -> dict[str, Any]:
    """Project Redis/factory counters into the frozen run-summary schema."""
    metrics = dict(factory_metrics or {})
    end = float(end_epoch if end_epoch is not None else datetime.now(timezone.utc).timestamp())
    start = float(start_epoch if start_epoch is not None else end)
    wall_time_sec = max(0.0, end - start)

    total_points_submitted = _coerce_int(metrics.get("submitted")) or 0
    total_points_finished = _coerce_int(metrics.get("ok")) or 0
    total_points_failed = _coerce_int(metrics.get("failed")) or 0
    completed_durations_sec = [
        float(item)
        for item in (metrics.get("completed_durations_sec") or [])
        if _coerce_float(item) is not None
    ]
    total_point_eval_sec = _coerce_float(metrics.get("total_point_eval_sec")) or 0.0
    external_tool_time_sec = max(0.0, float(external_tool_time_sec))
    time_in_framework_sec = max(0.0, float(total_point_eval_sec - external_tool_time_sec))
    framework_overhead_fraction = _safe_div(
        time_in_framework_sec,
        time_in_framework_sec + external_tool_time_sec,
    )

    retry_count = metrics.get("retry_count")
    if retry_count is not None:
        retry_count = _coerce_int(retry_count)

    resource = dict(resource_summary or {})
    summary = {
        "run_id": str(run_id or metrics.get("run_id") or "jarvis2-run"),
        "project_name": project_name,
        "sampler_name": sampler_name,
        "start_time": _utc_iso_from_epoch(start),
        "end_time": _utc_iso_from_epoch(end),
        "wall_time_sec": wall_time_sec,
        "total_points_submitted": total_points_submitted,
        "total_points_finished": total_points_finished,
        "total_points_failed": total_points_failed,
        "success_rate": _safe_div(total_points_finished, total_points_submitted),
        "throughput_points_per_min": _safe_div(total_points_finished * 60.0, wall_time_sec),
        "configured_workers": (
            configured_workers
            if configured_workers is not None
            else _coerce_int(metrics.get("configured_workers"))
        ),
        "peak_active_workers": _coerce_int(metrics.get("peak_active_workers")),
        "mean_active_workers": _coerce_float(metrics.get("mean_active_workers")),
        "avg_point_eval_sec": _mean(completed_durations_sec),
        "median_point_eval_sec": _median(completed_durations_sec),
        "time_in_external_tools_sec": external_tool_time_sec,
        "time_in_framework_sec": time_in_framework_sec,
        "framework_overhead_fraction": framework_overhead_fraction,
        "retry_count": retry_count,
        "crashed_subtasks": int(crashed_subtasks),
        "skipped_subtasks": skipped_subtasks,
        "cpu_percent_mean": _coerce_float(resource.get("cpu_percent_mean")),
        "memory_rss_mb_peak": _coerce_float(resource.get("memory_rss_mb_peak")),
        "open_file_peak": _coerce_int(resource.get("open_file_peak")),
    }
    validate_run_summary(summary)
    return summary


class RunSummaryRenderer:
    """Write run_summary.{json,csv,txt} using the frozen field order."""

    def render(self, summary: Mapping[str, Any]) -> str:
        lines = ["[Run Summary]"]
        for key in RUN_SUMMARY_FIELD_ORDER:
            value = summary.get(key)
            lines.append(f"{key}: {value}")
        return "\n".join(lines) + "\n"

    def write_outputs(
        self,
        summary: Mapping[str, Any],
        output_dir: str,
        *,
        rendered_text: str | None = None,
    ) -> dict[str, str]:
        output_root = Path(output_dir).expanduser().resolve()
        output_root.mkdir(parents=True, exist_ok=True)
        ordered = {key: summary.get(key) for key in RUN_SUMMARY_FIELD_ORDER}
        rendered = rendered_text if rendered_text is not None else self.render(summary)

        json_path = output_root / "run_summary.json"
        csv_path = output_root / "run_summary.csv"
        txt_path = output_root / "run_summary.txt"

        with open(json_path, "w", encoding="utf-8") as handle:
            json.dump(ordered, handle, indent=2, ensure_ascii=False)
            handle.write("\n")

        with open(csv_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(RUN_SUMMARY_FIELD_ORDER))
            writer.writeheader()
            writer.writerow(ordered)

        with open(txt_path, "w", encoding="utf-8") as handle:
            handle.write(rendered)

        return {
            "json": str(json_path),
            "csv": str(csv_path),
            "txt": str(txt_path),
        }


__all__ = [
    "RUN_SUMMARY_FIELD_ORDER",
    "RunSummaryRenderer",
    "build_run_summary",
    "validate_run_summary",
]