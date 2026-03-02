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
            self.last_time_mark = time.time()
            self._last_status_count = 0
            self._status_lock = threading.Lock()
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

    def submit_task(self, sample_info):
        if self.executor is None:
            self.get_executor()

        sample_uuid = (sample_info or {}).get("uuid", "UNKNOWN")
        self.log_executor.submit(
            self.logger.info,
            f"SUBMIT start -> uuid={sample_uuid}",
        )
        future = self.executor.submit(self.module_manager.execute_workflow, sample_info)

        def _on_done(done_future):
            exc = done_future.exception()
            if exc is not None:
                self.log_executor.submit(
                    self.logger.error,
                    f"future exception consumed: uuid={sample_uuid} error={exc}",
                )
                self.log_executor.submit(
                    self.logger.info,
                    f"SUBMIT end -> uuid={sample_uuid} status=failed",
                )
            else:
                self.log_executor.submit(
                    self.logger.info,
                    f"SUBMIT end -> uuid={sample_uuid} status=ok",
                )

        future.add_done_callback(_on_done)

        with self._status_lock:
            self.task_count += 1
            current_count = self.task_count
        self.print_status(current_count)
        return future

    def _log_status(self, task_count, elapsed_seconds, total_seconds):
        from jarvishep.utils import format_duration

        self.logger.warning(
            "Submitted {} tasks, time for last 100 tasks: {}, total time: {}".format(
                task_count,
                format_duration(elapsed_seconds),
                format_duration(total_seconds),
            )
        )

    def print_status(self, task_count=None):
        task_count = self.task_count if task_count is None else task_count
        if task_count % 100 != 0:
            return

        with self._status_lock:
            if task_count <= self._last_status_count:
                return
            end_time = time.time()
            elapsed_seconds = end_time - self.last_time_mark
            self.task_time += elapsed_seconds
            total_seconds = self.task_time
            self.last_time_mark = end_time
            self._last_status_count = task_count

        if getattr(self, "log_executor", None) is not None:
            self.log_executor.submit(self._log_status, task_count, elapsed_seconds, total_seconds)
        else:
            self._log_status(task_count, elapsed_seconds, total_seconds)

    def shutdown(self, wait=True, cancel_futures=True):
        if getattr(self, "executor", None) is not None:
            self.executor.shutdown(wait=wait, cancel_futures=cancel_futures)
            self.executor = None

        if getattr(self, "log_executor", None) is not None:
            self.log_executor.shutdown(wait=wait, cancel_futures=cancel_futures)
            self.log_executor = None
