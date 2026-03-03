#!/usr/bin/env python3
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import threading
import time 


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
            self.executor = ThreadPoolExecutor(max_workers=max_workers)
            self.log_executor = ThreadPoolExecutor(max_workers=1)
            self.task_count = 0
            self.task_time = 0
            self.completed_ok = 0
            self.completed_failed = 0
            self.last_time_mark = time.time()
            self._last_status_count = 0
            self._status_lock = threading.Lock()
            # Status interval is monotonic: only degrades to larger windows.
            # Base policy:
            #   < 10k tasks  -> 100
            #   >=10k        -> 1000
            #   >=50k        -> 5000
            #   >=100k       -> 10000
            # then continue similarly with x5/x2 growth.
            self._status_interval = 100

            # After long-running scans, switch to minute-level summary logs.
            self._factory_start_ts = time.monotonic()
            self._low_throughput_after_seconds = 5 * 60.0
            self._low_throughput_summary_seconds = 60.0
            self._low_throughput_enabled = False
            self._last_low_throughput_log_ts = self._factory_start_ts
            self.initialized = True

    def get_executor(self):
        if self.executor is None:
            self.executor = ThreadPoolExecutor(max_workers=self._max_workers)
        return self.executor

    def set_config(self, config):
        self.config = config

    def set_logger(self, logger):
        self.logger = logger
        if getattr(self, "log_executor", None) is not None:
            self.log_executor.submit(self.logger.warning, "Building the factory for workers ...")
        else:
            self.logger.warning("Building the factory for workers ...")
        self.info   = {}

    @staticmethod
    def _compute_task_count_interval(task_count):
        task_count = int(task_count)
        if task_count < 10_000:
            return 100

        threshold = 10_000
        interval = 1_000
        multipliers = (5, 2)
        level = 0

        while True:
            mult = multipliers[level % len(multipliers)]
            next_threshold = threshold * mult
            if task_count < next_threshold:
                return int(interval)
            threshold = next_threshold
            interval *= mult
            level += 1

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
        if self.executor is None:
            self.get_executor()

        sample_uuid = (sample_info or {}).get("uuid", "UNKNOWN")
        future = self.executor.submit(self.module_manager.execute_workflow, sample_info)

        def _on_done(done_future):
            exc = done_future.exception()
            with self._status_lock:
                if exc is not None:
                    self.completed_failed += 1
                else:
                    self.completed_ok += 1
            if exc is not None:
                self.log_executor.submit(
                    self.logger.error,
                    f"future exception consumed: uuid={sample_uuid} error={exc}",
                )

        future.add_done_callback(_on_done)

        with self._status_lock:
            self.task_count += 1
            current_count = self.task_count
        self._update_status_interval(current_count)
        self.print_status(current_count)
        return future

    def _log_status(
        self,
        task_count,
        elapsed_seconds,
        total_seconds,
        window_tasks,
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

        self.logger.warning(
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
            total_seconds = float(self.task_time + max(0.0, now - self.last_time_mark))

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
        }

    def _log_shutdown_summary(self, summary):
        from jarvishep.utils import format_duration

        self.logger.warning(
            "Factory shutdown summary ->\n"
            "\tsubmitted   -> {}\n"
            "\tok          -> {}\n"
            "\tfailed      -> {}\n"
            "\tunfinished  -> {}\n"
            "\ttail_tasks  -> {}\n"
            "\ttotal_time  -> {}\n"
            "\tavg_rate    -> {:.2f} tasks/s".format(
                summary["submitted"],
                summary["ok"],
                summary["failed"],
                summary["unfinished"],
                summary["tail_tasks"],
                format_duration(summary["total_seconds"]),
                summary["avg_rate"],
            )
        )

    def print_status(self, task_count=None):
        task_count = self.task_count if task_count is None else task_count
        now_monotonic = time.monotonic()
        low_throughput_mode = self._update_low_throughput_mode(now_monotonic=now_monotonic)

        with self._status_lock:
            interval = self._status_interval
            previous_count = self._last_status_count

            if low_throughput_mode:
                if (now_monotonic - self._last_low_throughput_log_ts) < self._low_throughput_summary_seconds:
                    return
            elif (task_count - previous_count) < interval:
                return

            if not low_throughput_mode and task_count <= self._last_status_count:
                return

            end_time = time.time()
            elapsed_seconds = end_time - self.last_time_mark
            self.task_time += elapsed_seconds
            total_seconds = self.task_time
            window_tasks = max(0, task_count - previous_count)
            self.last_time_mark = end_time
            if task_count > self._last_status_count:
                self._last_status_count = task_count
            completed_ok = self.completed_ok
            completed_failed = self.completed_failed
            unfinished = max(0, task_count - completed_ok - completed_failed)
            if low_throughput_mode:
                self._last_low_throughput_log_ts = now_monotonic

        if getattr(self, "log_executor", None) is not None:
            self.log_executor.submit(
                self._log_status,
                task_count,
                elapsed_seconds,
                total_seconds,
                window_tasks,
                low_throughput_mode,
                completed_ok,
                completed_failed,
                unfinished,
            )
        else:
            self._log_status(
                task_count,
                elapsed_seconds,
                total_seconds,
                window_tasks,
                low_throughput_mode,
                completed_ok,
                completed_failed,
                unfinished,
            )

    def shutdown(self, wait=True, cancel_futures=True):
        if getattr(self, "executor", None) is not None:
            self.executor.shutdown(wait=wait, cancel_futures=cancel_futures)
            self.executor = None

        summary = None
        if hasattr(self, "logger"):
            summary = self._build_shutdown_summary()

        if getattr(self, "log_executor", None) is not None:
            if summary is not None:
                # Must be emitted before log executor shutdown even when
                # `cancel_futures=True` is requested by caller.
                self._log_shutdown_summary(summary)
            self.log_executor.shutdown(wait=wait, cancel_futures=cancel_futures)
            self.log_executor = None
        elif summary is not None:
            self._log_shutdown_summary(summary)
