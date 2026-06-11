#!/usr/bin/env python3
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import threading
import time 

from jarvishep import benchmark
from jarvishep.log_kv import format_two_column_log
from jarvishep.sample import materialize_failure_artifacts


class WorkerFactory:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(WorkerFactory, cls).__new__(cls)
                cls._instance.module_manager = args[0] if args else None
                cls._instance.executor = None
                cls._instance._max_workers = 4
        return cls._instance

    def configure(self, module_manager=None, max_workers=4):
        if not hasattr(self, "initialized"):
            self.module_manager = module_manager if module_manager else self.module_manager
            self._max_workers = max_workers
            self.executor = ThreadPoolExecutor(
                max_workers=max_workers,
                thread_name_prefix="jarvis-hep-factory",
            )
            self.task_count = 0
            self.task_time = 0
            self.completed_ok = 0
            self.completed_failed = 0
            self._status_start_time = time.time()
            self.last_time_mark = self._status_start_time
            self._last_status_count = 0
            self._info_status_interval = 100
            self._last_info_status_count = 0
            self._last_info_time_mark = self._status_start_time
            self._status_lock = threading.Lock()
            # Regular progress logs use two fixed cadences in normal-throughput mode:
            #   INFO    -> every 100 submitted tasks
            #   WARNING -> every 1000 submitted tasks
            self._status_interval = 1000

            # After long-running scans, switch to minute-level summary logs.
            self._factory_start_ts = time.monotonic()
            self._low_throughput_after_seconds = 5 * 60.0
            self._low_throughput_summary_seconds = 60.0
            self._low_throughput_enabled = False
            self._last_low_throughput_log_ts = self._factory_start_ts
            self.subprocess_scheduler = None
            self.run_summary_collector = None
            self._active_workers = 0
            self._peak_active_workers = 0
            self._active_worker_integral = 0.0
            self._active_worker_last_ts = self._factory_start_ts
            self._total_point_eval_sec = 0.0
            self._completed_durations_sec = []
            self._retry_metric_supported = False
            self.retry_count = 0
            self.initialized = True

    def get_executor(self):
        if self.executor is None:
            self.executor = ThreadPoolExecutor(
                max_workers=self._max_workers,
                thread_name_prefix="jarvis-hep-factory",
            )
        return self.executor

    def set_config(self, config):
        self.config = config

    def set_logger(self, logger):
        self.logger = logger
        self.logger.warning("Building the factory for workers ...")
        self.info   = {}

    def set_subprocess_scheduler(self, scheduler):
        self.subprocess_scheduler = scheduler

    def get_subprocess_scheduler(self):
        return getattr(self, "subprocess_scheduler", None)

    def set_run_summary_collector(self, collector):
        self.run_summary_collector = collector

    def _update_active_worker_integral_locked(self, now_monotonic):
        now_monotonic = float(now_monotonic)
        elapsed = max(0.0, now_monotonic - self._active_worker_last_ts)
        self._active_worker_integral += elapsed * float(self._active_workers)
        self._active_worker_last_ts = now_monotonic

    @staticmethod
    def _extract_retry_count(sample_info):
        if not isinstance(sample_info, dict):
            return None

        candidates = []
        if "NAttempt" in sample_info:
            candidates.append(sample_info.get("NAttempt"))

        nuisance = sample_info.get("nuisance")
        if isinstance(nuisance, dict) and "NAttempt" in nuisance:
            candidates.append(nuisance.get("NAttempt"))

        attempts = None
        for candidate in candidates:
            try:
                value = int(candidate)
            except Exception:
                continue
            attempts = value if attempts is None else max(attempts, value)

        if attempts is None:
            return None
        return max(0, attempts - 1)

    def _module_failure_policy(self) -> str:
        manager = getattr(self, "module_manager", None)
        if manager is not None and hasattr(manager, "_module_failure_policy"):
            try:
                policy = str(manager._module_failure_policy()).strip().lower()
            except Exception:
                return "fail-fast"
            if policy in {"continue", "continue-on-error"}:
                return "continue"
        return "fail-fast"

    def _execute_workflow_tracked(self, sample_info):
        start_monotonic = time.monotonic()
        retry_count = self._extract_retry_count(sample_info)

        with self._status_lock:
            self._update_active_worker_integral_locked(start_monotonic)
            self._active_workers += 1
            if self._active_workers > self._peak_active_workers:
                self._peak_active_workers = self._active_workers

        success = False
        try:
            result = self.module_manager.execute_workflow(sample_info)
            if isinstance(sample_info, dict) and str(sample_info.get("status", "")).strip().lower() == "failed":
                materialize_failure_artifacts(sample_info, error="workflow module failure")
            success = True
            return result
        except Exception as exc:
            if isinstance(sample_info, dict):
                materialize_failure_artifacts(sample_info, error=exc)
            if self._module_failure_policy() == "continue":
                sample_uuid = (sample_info or {}).get("uuid", "UNKNOWN")
                if isinstance(sample_info, dict):
                    observables = sample_info.setdefault("observables", {})
                    if isinstance(observables, dict):
                        observables.setdefault("uuid", sample_uuid)
                    sample_info["status"] = "Failed"
                message = format_two_column_log(
                    "Workflow exception downgraded by ModuleFailurePolicy=continue",
                    [("uuid", sample_uuid), ("error", exc)],
                )
                self.logger.error(message)
                return 1.0
            raise
        finally:
            end_monotonic = time.monotonic()
            duration_sec = max(0.0, end_monotonic - start_monotonic)
            with self._status_lock:
                self._update_active_worker_integral_locked(end_monotonic)
                self._active_workers = max(0, self._active_workers - 1)
                self._total_point_eval_sec += duration_sec
                if success:
                    self.completed_ok += 1
                    self._completed_durations_sec.append(duration_sec)
                else:
                    self.completed_failed += 1
                if retry_count is not None:
                    self._retry_metric_supported = True
                    self.retry_count += retry_count

    @staticmethod
    def _compute_task_count_interval(task_count):
        return 1_000

    def _update_status_interval(self, task_count):
        target = self._compute_task_count_interval(task_count)
        with self._status_lock:
            if target > self._status_interval:
                self._status_interval = target

    def _update_low_throughput_mode(self, now_monotonic=None):
        now_monotonic = time.monotonic() if now_monotonic is None else float(now_monotonic)
        with self._status_lock:
            if (
                not self._low_throughput_enabled
                and (now_monotonic - self._factory_start_ts) >= self._low_throughput_after_seconds
            ):
                self._low_throughput_enabled = True
                # Allow first summary output immediately after mode switch.
                self._last_low_throughput_log_ts = now_monotonic - self._low_throughput_summary_seconds
            enabled = self._low_throughput_enabled

        return enabled

    def submit_task(self, sample_info):
        timing_enabled = benchmark.TIMING_ENABLED
        timing_start = benchmark.monotonic_seconds() if timing_enabled else None
        try:
            if self.executor is None:
                self.get_executor()

            sample_uuid = (sample_info or {}).get("uuid", "UNKNOWN")
            if isinstance(sample_info, dict) and self.run_summary_collector is not None:
                sample_info.setdefault("run_summary_collector", self.run_summary_collector)
            future = self.executor.submit(self._execute_workflow_tracked, sample_info)

            def _on_done(done_future):
                try:
                    exc = done_future.exception()
                except Exception as done_exc:
                    exc = done_exc
                if exc is not None:
                    message = format_two_column_log(
                        "future exception consumed",
                        [("uuid", sample_uuid), ("error", exc)],
                    )
                    self.logger.error(message)

            future.add_done_callback(_on_done)

            with self._status_lock:
                self.task_count += 1
                current_count = self.task_count
            self._update_status_interval(current_count)
            self.print_status(current_count)
            return future
        finally:
            if timing_enabled and timing_start is not None:
                benchmark.record_stage(
                    "submit_future_locks",
                    benchmark.monotonic_seconds() - timing_start,
                )

    def _log_status(
        self,
        task_count,
        elapsed_seconds,
        total_seconds,
        window_tasks,
        level="WARNING",
        low_throughput_mode=False,
        completed_ok=0,
        completed_failed=0,
        unfinished=0,
    ):
        from jarvishep.utils import format_duration

        if low_throughput_mode:
            window_rate = (float(window_tasks) / float(elapsed_seconds)) if elapsed_seconds > 0 else 0.0
            total_rate = (float(task_count) / float(total_seconds)) if total_seconds > 0 else 0.0
            self.logger.warning(
                "Factory low-throughput summary ->\n"
                "\twindow_time     -> {}\n"
                "\twindow_tasks    -> {}\n"
                "\tsubmitted_total -> {}\n"
                "\tok_total        -> {}\n"
                "\tfailed_total    -> {}\n"
                "\tunfinished      -> {}\n"
                "\twindow_rate     -> {:.2f} tasks/s\n"
                "\ttotal_rate      -> {:.2f} tasks/s".format(
                    format_duration(elapsed_seconds),
                    int(window_tasks),
                    int(task_count),
                    int(completed_ok),
                    int(completed_failed),
                    int(unfinished),
                    window_rate,
                    total_rate,
                )
            )
            return

        log_method = self.logger.warning
        if str(level).upper() == "INFO":
            log_method = self.logger.info

        log_method(
            "Submitted {} tasks, time for last {} tasks: {}, total time: {}".format(
                task_count,
                window_tasks,
                format_duration(elapsed_seconds),
                format_duration(total_seconds),
            )
        )

    def _build_shutdown_summary(self):
        now = time.time()
        with self._status_lock:
            task_count = int(self.task_count)
            completed_ok = int(self.completed_ok)
            completed_failed = int(self.completed_failed)
            last_status_count = int(self._last_status_count)
            total_seconds = float(max(0.0, now - self._status_start_time))

        tail_tasks = max(0, task_count - last_status_count)
        unfinished = max(0, task_count - completed_ok - completed_failed)
        avg_rate = (float(task_count) / total_seconds) if total_seconds > 0 else 0.0
        return {
            "submitted": task_count,
            "ok": completed_ok,
            "failed": completed_failed,
            "unfinished": unfinished,
            "tail_tasks": tail_tasks,
            "total_seconds": total_seconds,
            "avg_rate": avg_rate,
            "subprocess": (
                self.get_subprocess_scheduler().snapshot()
                if self.get_subprocess_scheduler() is not None
                else None
            ),
        }

    def get_run_metrics(self):
        now_monotonic = time.monotonic()
        with self._status_lock:
            self._update_active_worker_integral_locked(now_monotonic)
            elapsed = max(0.0, now_monotonic - self._factory_start_ts)
            mean_active_workers = (
                self._active_worker_integral / elapsed
                if elapsed > 0.0
                else 0.0
            )
            retry_count = self.retry_count if self._retry_metric_supported else None
            completed_durations_sec = list(self._completed_durations_sec)

            return {
                "submitted": int(self.task_count),
                "ok": int(self.completed_ok),
                "failed": int(self.completed_failed),
                "configured_workers": int(self._max_workers),
                "peak_active_workers": int(self._peak_active_workers),
                "mean_active_workers": float(mean_active_workers),
                "total_point_eval_sec": float(self._total_point_eval_sec),
                "completed_durations_sec": completed_durations_sec,
                "retry_count": retry_count,
            }

    def _log_shutdown_summary(self, summary):
        from jarvishep.utils import format_duration

        self.logger.warning(
            format_two_column_log(
                "Factory shutdown summary",
                [
                    ("submitted", summary["submitted"]),
                    ("ok", summary["ok"]),
                    ("failed", summary["failed"]),
                    ("unfinished", summary["unfinished"]),
                    ("tail_tasks", summary["tail_tasks"]),
                    ("total_time", format_duration(summary["total_seconds"])),
                    ("avg_rate", f"{summary['avg_rate']:.2f} tasks/s"),
                ],
            )
        )
        sp = summary.get("subprocess")
        if isinstance(sp, dict):
            self.logger.warning(
                format_two_column_log(
                    "Subprocess scheduler summary",
                    [
                        ("submitted", sp.get("submitted", 0)),
                        ("completed", sp.get("completed", 0)),
                        ("failed", sp.get("failed", 0)),
                        ("timed_out", sp.get("timed_out", 0)),
                        ("peak_running", sp.get("peak_running", 0)),
                        ("pending", sp.get("pending", 0)),
                        ("fd_count", sp.get("fd_count", "N/A")),
                        (
                            "rss_mb",
                            "N/A"
                            if sp.get("rss_mb") is None
                            else f"{sp.get('rss_mb'):.2f}",
                        ),
                    ],
                )
            )

    def print_status(self, task_count=None):
        task_count = self.task_count if task_count is None else task_count
        now_monotonic = time.monotonic()
        low_throughput_mode = self._update_low_throughput_mode(now_monotonic=now_monotonic)

        with self._status_lock:
            interval = self._status_interval
            previous_count = self._last_status_count
            info_previous_count = self._last_info_status_count

            if low_throughput_mode:
                if (now_monotonic - self._last_low_throughput_log_ts) < self._low_throughput_summary_seconds:
                    return
            else:
                info_due = (
                    (task_count - info_previous_count) >= self._info_status_interval
                    and task_count > info_previous_count
                )
                warning_due = (
                    (task_count - previous_count) >= interval
                    and task_count > previous_count
                )
                if not info_due and not warning_due:
                    return

            end_time = time.time()
            completed_ok = self.completed_ok
            completed_failed = self.completed_failed
            unfinished = max(0, task_count - completed_ok - completed_failed)
            if low_throughput_mode:
                elapsed_seconds = end_time - self.last_time_mark
                total_seconds = max(0.0, end_time - self._status_start_time)
                window_tasks = max(0, task_count - previous_count)
                self.last_time_mark = end_time
                if task_count > self._last_status_count:
                    self._last_status_count = task_count
                self._last_low_throughput_log_ts = now_monotonic
                info_payload = None
                warning_payload = (
                    task_count,
                    elapsed_seconds,
                    total_seconds,
                    window_tasks,
                    "WARNING",
                    True,
                    completed_ok,
                    completed_failed,
                    unfinished,
                )
            else:
                total_seconds = max(0.0, end_time - self._status_start_time)
                info_payload = None
                warning_payload = None

                if info_due and task_count > self._last_info_status_count:
                    info_elapsed_seconds = end_time - self._last_info_time_mark
                    info_window_tasks = max(0, task_count - self._last_info_status_count)
                    self._last_info_time_mark = end_time
                    self._last_info_status_count = task_count
                    info_payload = (
                        task_count,
                        info_elapsed_seconds,
                        total_seconds,
                        info_window_tasks,
                        "INFO",
                        False,
                        completed_ok,
                        completed_failed,
                        unfinished,
                    )

                if warning_due and task_count > self._last_status_count:
                    warning_elapsed_seconds = end_time - self.last_time_mark
                    warning_window_tasks = max(0, task_count - self._last_status_count)
                    self.last_time_mark = end_time
                    self._last_status_count = task_count
                    warning_payload = (
                        task_count,
                        warning_elapsed_seconds,
                        total_seconds,
                        warning_window_tasks,
                        "WARNING",
                        False,
                        completed_ok,
                        completed_failed,
                        unfinished,
                    )

        payloads = []
        if info_payload is not None:
            payloads.append(info_payload)
        if warning_payload is not None:
            payloads.append(warning_payload)

        if not payloads:
            return

        for payload in payloads:
            self._log_status(*payload)

    def shutdown(self, wait=True, cancel_futures=True):
        if getattr(self, "executor", None) is not None:
            self.executor.shutdown(wait=wait, cancel_futures=cancel_futures)
            self.executor = None

        summary = None
        if hasattr(self, "logger"):
            summary = self._build_shutdown_summary()

        if self.get_subprocess_scheduler() is not None:
            try:
                self.get_subprocess_scheduler().shutdown(wait=wait, timeout=60.0)
            except Exception as exc:
                if hasattr(self, "logger"):
                    self.logger.error(f"Subprocess scheduler shutdown failed -> {exc}")
            finally:
                self.subprocess_scheduler = None

        if summary is not None:
            self._log_shutdown_summary(summary)
