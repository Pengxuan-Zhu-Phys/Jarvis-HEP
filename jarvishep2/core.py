#!/usr/bin/env python3
"""Jarvis2 control-process orchestrator for the distributed runtime."""

from __future__ import annotations

import os
import time
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.archiver import ArchiverProcess, SimpleArchiver
from jarvishep2.command_parser import CommandParser
from jarvishep2.factory import TaskFactory
from jarvishep2.command_parser import prepare_calculator_modules
from jarvishep2.worker_config import build_command_parser, build_worker_config
from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging
from jarvishep2.redis_queue import RedisQueue, make_fakeredis_queue
from jarvishep2.runtime_config import (
    get_archiver_config,
    get_delete_method,
    get_runtime_block,
    get_watchdog_config,
)
from jarvishep2.dashboard import SnapshotReader, format_monitor_view
from jarvishep2.monitoring.run_summary import RunSummaryRenderer, build_run_summary
from jarvishep2.sample import Sample


class Jarvis2Core:
    """Minimal distributed run orchestrator for the D1.1 opera-only MVP."""

    def __init__(self, config: Mapping[str, Any] | None = None) -> None:
        self.config = dict(config or {})
        self.runtime = get_runtime_block(self.config)
        self.info: dict[str, Any] = {}
        self.redis: RedisQueue | None = None
        self.factory: TaskFactory | None = None
        self.archiver: SimpleArchiver | ArchiverProcess | None = None
        self.sampler: Any = None
        self.command_parser: CommandParser | None = None
        self._logger = get_jarvis_logger("core")

    def init_logger(self) -> None:
        setup_jarvis_logging(role="core")

    def is_redis_runtime(self) -> bool:
        """Return True when the distributed Redis path should be used."""
        return str(self.runtime.get("mode", "auto")).strip().lower() == "redis"

    def init_redis(self, *, client: Any = None) -> RedisQueue:
        redis_config = dict(self.runtime.get("redis") or {})
        if client is not None:
            self.redis = RedisQueue(redis_config, client=client)
        elif not redis_config:
            self.redis = make_fakeredis_queue()
        else:
            self.redis = RedisQueue(redis_config)
        self.redis.connect()
        return self.redis

    def init_command_parser(self) -> CommandParser:
        """Run Phase-1 static command resolution for the loaded task config."""
        self.command_parser = build_command_parser(self.config)
        return self.command_parser

    @staticmethod
    def _command_parser_payload(parser: CommandParser) -> dict[str, Any]:
        return {
            "project_root": parser.project_root,
            "scan_name": parser.scan_name,
            "libdeps_paths": dict(parser.libdeps_paths),
            "registered": {
                name: {
                    "name": item.name,
                    "path": item.path,
                    "resolution": item.resolution,
                }
                for name, item in parser.registered.items()
            },
            "registered_symlink_root": parser.registered_symlink_root,
        }

    def _apply_command_parser_to_worker_config(self, worker_config: dict[str, Any]) -> dict[str, Any]:
        if self.command_parser is None:
            self.init_command_parser()
        merged = dict(worker_config)
        calculator_modules = merged.get("calculator_modules")
        if calculator_modules:
            merged["calculator_modules"] = prepare_calculator_modules(
                calculator_modules,
                self.command_parser,
            )
        merged["command_parser"] = self._command_parser_payload(self.command_parser)
        return merged

    def build_worker_config(self, **overrides: Any) -> dict[str, Any]:
        """Build a picklable Worker blueprint with Phase-1 command resolution applied."""
        task_result_dir = str(
            overrides.pop("task_result_dir", None)
            or self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        return build_worker_config(
            self.config,
            task_result_dir=task_result_dir,
            parser=self.command_parser,
            calculator_modules=overrides.pop("calculator_modules", None),
            likelihood_expressions=overrides.pop("likelihood_expressions", None),
            opera_modules=overrides.pop("opera_modules", None),
            sample_dirs=overrides.pop("sample_dirs", None),
            extra=overrides or None,
        )

    def init_factory(self, worker_config: Mapping[str, Any] | None = None) -> TaskFactory | None:
        if not self.is_redis_runtime():
            self._logger.info("Runtime.mode != redis; skipping TaskFactory bring-up")
            return None

        redis_config = dict(self.runtime.get("redis") or {})
        self.factory = TaskFactory.get_instance(redis_config)
        if self.redis is not None:
            self.factory.redis = self.redis
        else:
            self.factory.init_redis()

        workers = int(self.runtime.get("workers", 1) or 1)
        if workers <= 0:
            workers = 1

        if worker_config is None:
            if self.command_parser is None:
                self.init_command_parser()
            merged_config = self.build_worker_config()
        else:
            merged_config = dict(worker_config)
            if "command_parser" not in merged_config:
                merged_config = self._apply_command_parser_to_worker_config(merged_config)
        self.factory.start_workers(workers, **merged_config)
        self.factory.start_monitor(update_hz=120.0)
        watchdog = get_watchdog_config(self.config)
        self.factory.start_watchdog(**watchdog)
        self._logger.info("TaskFactory started with %d worker(s)", workers)
        return self.factory

    def init_archiver(self, db_path: str | None = None) -> SimpleArchiver | ArchiverProcess:
        if self.redis is None:
            raise RuntimeError("init_redis() must run before init_archiver()")
        task_result_dir = str(
            self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        database_dir = os.path.join(task_result_dir, "DATABASE")
        sample_root = os.path.join(task_result_dir, "SAMPLE")
        os.makedirs(database_dir, exist_ok=True)
        os.makedirs(sample_root, exist_ok=True)
        resolved_db_path = db_path or os.path.join(database_dir, "samples.hdf5")
        archiver_config = get_archiver_config(self.config)
        delete_method = get_delete_method(self.config)
        redis_config = dict(self.runtime.get("redis") or self.redis.connection_config())

        if str(archiver_config.get("mode", "thread")).strip().lower() == "process":
            self.archiver = ArchiverProcess(
                redis_config,
                db_path=resolved_db_path,
                sample_root=sample_root,
                delete_method=delete_method,
                archiver_config=archiver_config,
            )
            self.archiver.start()
        else:
            self.archiver = SimpleArchiver(
                self.redis,
                resolved_db_path,
                sample_root=sample_root,
                delete_method=delete_method,
                archiver_config=archiver_config,
            )
            self.archiver.start()
        return self.archiver

    def _archiver_records_written(self) -> int:
        archiver = self.archiver
        if archiver is None:
            return 0
        counter = getattr(archiver, "records_written", 0)
        if hasattr(counter, "value"):
            return int(counter.value)
        return int(counter)

    def set_sampler(self, sampler: Any) -> None:
        self.sampler = sampler
        if self.redis is not None and hasattr(sampler, "set_redis"):
            sampler.set_redis(self.redis)

    def submit_samples(self, samples: Sequence[Sample]) -> None:
        if self.sampler is None:
            raise RuntimeError("sampler is not configured")
        if hasattr(self.sampler, "_submit_group"):
            self.sampler._submit_group(list(samples))
        else:
            for sample in samples:
                self.sampler._submit(sample)

    def wait_for_results(
        self,
        expected: int,
        *,
        timeout: float = 30.0,
        poll_interval: float = 0.1,
    ) -> None:
        if self.archiver is None:
            raise RuntimeError("archiver is not configured")
        deadline = time.monotonic() + max(0.1, float(timeout))
        while time.monotonic() < deadline:
            if self._archiver_records_written() >= expected:
                return
            time.sleep(poll_interval)
        raise TimeoutError(
            f"timed out waiting for {expected} archived results; "
            f"got {self._archiver_records_written()}"
        )

    def get_monitor_snapshot(self) -> dict[str, Any]:
        if self.factory is None:
            return {}
        return self.factory.get_monitor_snapshot()

    def monitor_once(self) -> str:
        if self.factory is None:
            return ""
        view = SnapshotReader(self.factory).read()
        return format_monitor_view(view)

    def write_run_summary(self, output_dir: str | None = None) -> dict[str, str]:
        if self.factory is None:
            raise RuntimeError("factory is not configured")
        task_result_dir = str(
            output_dir
            or self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        metrics = dict(self.factory.get_run_metrics())
        started = metrics.pop("run_started_at", None)
        scan = self.config.get("Scan") or {}
        scan_name = scan.get("name") if isinstance(scan, Mapping) else None
        sampler_name = type(self.sampler).__name__ if self.sampler is not None else None
        summary = build_run_summary(
            factory_metrics=metrics,
            project_name=str(self.config.get("project_name") or scan_name or ""),
            sampler_name=sampler_name,
            run_id=str(self.info.get("run_id") or metrics.get("run_id") or "jarvis2-run"),
            start_epoch=float(started) if started is not None else None,
            configured_workers=int(self.runtime.get("workers", 0) or len(self.factory.workers)),
        )
        return RunSummaryRenderer().write_outputs(summary, task_result_dir)

    def shutdown(self, *, wait: bool = True, write_run_summary: bool = False) -> None:
        if self.archiver is not None:
            if isinstance(self.archiver, ArchiverProcess):
                self.archiver.stop(wait=wait)
            else:
                self.archiver.stop(wait=wait, drain=True)
            self.archiver = None
        if write_run_summary and self.factory is not None:
            try:
                self.write_run_summary()
            except Exception as exc:
                self._logger.warning("run_summary write failed -> %s", exc)
        if self.factory is not None:
            self.factory.shutdown(wait=wait)
            self.factory = None
        elif self.redis is not None:
            self.redis.close()
        self.redis = None


__all__ = ["Jarvis2Core"]