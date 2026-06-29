#!/usr/bin/env python3
"""Redis-backed Worker process for Jarvis-HEP V2."""

from __future__ import annotations

import concurrent.futures
import json
import os
import signal
import threading
import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.async_subprocess import AsyncSubprocessScheduler, SubprocessRuntimeConfig
from jarvishep2.command_parser import CommandParser
from jarvishep2.env_setup import EnvCapture, resolve_env_setup_sources
from jarvishep2.archive_handoff import list_product_names, resolve_staging_dir, stage_sample_dir
from jarvishep2.file_ops import DEFAULT_DELETE_METHOD, delete_paths, normalize_delete_method
from jarvishep2.likelihood import LogLikelihoodEvaluator
from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging
from jarvishep2.mapper import build_mapper
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.operas import preload_operas
from jarvishep2.redis_queue import RedisQueue
from jarvishep2.sample import Sample, materialize_failure_artifacts
from jarvishep2.sample import ExecutionStep
from jarvishep2.workflow import group_by_layer, max_layer_width, resolve_module_layers

_SPAWN_CTX = get_spawn_context()
Process = _SPAWN_CTX.Process


class Worker(Process):
    """Long-lived process that pulls one Sample at a time from Redis.

    ``__init__`` accepts either a live :class:`RedisQueue` (from the Factory) or
    a picklable connection mapping. Only the connection **settings** are stored
    for spawn; the child process builds its own Redis client in ``run()``.
    """

    def __init__(
        self,
        worker_id: int,
        redis: RedisQueue | Mapping[str, Any],
        worker_config: Mapping[str, Any],
    ) -> None:
        super().__init__(name=f"HEP2-Worker-{worker_id}", daemon=False)
        self.worker_id = int(worker_id)
        self.redis_config = RedisQueue.extract_connection_config(redis)
        self.worker_config = dict(worker_config)
        self._redis: RedisQueue | None = None
        self._mapper = None
        self._operas: dict[str, Any] = {}
        self._calculators: dict[str, CalculatorModule] = {}
        self._likelihood: LogLikelihoodEvaluator | None = None
        self._scheduler: AsyncSubprocessScheduler | None = None
        self._command_parser: CommandParser | None = None
        self._delete_method = DEFAULT_DELETE_METHOD
        self._staging_dir = ""
        self._handoff_to_staging = True
        self._observables_lock: threading.Lock | None = None
        self._is_running = True
        self._current_sample_uuid: str | None = None

    def _handle_signal(self, signum: int, _frame: Any) -> None:
        self._is_running = False

    def _init_redis(self) -> None:
        """Connect to Redis in the child process (spawn-safe)."""
        self._redis = RedisQueue(self.redis_config)
        self._redis.connect()
        self._redis.heartbeat(
            str(self.worker_id),
            status="idle",
            pid=self.pid,
        )

    def _init_runtime(self) -> None:
        self._observables_lock = threading.Lock()
        self._delete_method = normalize_delete_method(
            self.worker_config.get("delete_method", DEFAULT_DELETE_METHOD)
        )
        self._staging_dir = str(self.worker_config.get("staging_dir") or "").strip()
        if not self._staging_dir:
            sample_config = self.worker_config.get("sample_config") or {}
            task_result_dir = str(sample_config.get("task_result_dir") or "").strip()
            if task_result_dir:
                self._staging_dir = resolve_staging_dir(task_result_dir)
        if "handoff_to_staging" in self.worker_config:
            self._handoff_to_staging = bool(self.worker_config.get("handoff_to_staging"))
        else:
            cleanup = self.worker_config.get("cleanup_config") or {}
            strategy = str(cleanup.get("strategy", "mv_to_staging")).strip().lower()
            self._handoff_to_staging = strategy == "mv_to_staging"
        self._mapper = build_mapper(self.worker_config.get("mapper"))
        opera_configs = self.worker_config.get("opera_modules") or {}
        if isinstance(opera_configs, list):
            opera_configs = {
                str(item.get("name", f"Operas{index}")): item
                for index, item in enumerate(opera_configs)
            }
        self._operas = preload_operas(opera_configs)
        calculator_configs = self.worker_config.get("calculator_modules") or []
        if isinstance(calculator_configs, dict):
            calculator_configs = list(calculator_configs.values())
        self._calculators = CalculatorModule.from_config_list(calculator_configs)
        layer_width = max(
            1,
            int(self.worker_config.get("subprocess_max_concurrency", 0) or 0),
            max_layer_width(resolve_module_layers(calculator_configs)),
        )
        self._scheduler = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(
                max_concurrency=layer_width,
                log_policy="quiet",
                progress_interval_sec=3600.0,
            ),
            logger=get_jarvis_logger("worker").bind(worker_id=self.worker_id),
        )
        self._scheduler.start()
        self._command_parser = CommandParser.from_picklable(self.worker_config.get("command_parser"))
        for module in self._calculators.values():
            module.attach_scheduler(self._scheduler)
            module.attach_command_parser(self._command_parser)
            sources = resolve_env_setup_sources(
                module.env_setup,
                command_parser=self._command_parser,
            )
            if sources:
                module.bind_env(EnvCapture.merged_env(sources))
        likelihood_exprs = self.worker_config.get("likelihood_expressions") or []
        self._likelihood = LogLikelihoodEvaluator(likelihood_exprs)

    def _heartbeat(self, status: str) -> None:
        if self._redis is None:
            return
        self._redis.heartbeat(
            str(self.worker_id),
            status=status,
            pid=self.pid,
            current_sample=self._current_sample_uuid,
            ts=time.time(),
        )

    def _run_calculator_step(self, step_name: str, sample: Sample) -> None:
        module = self._calculators.get(step_name)
        if module is None:
            raise KeyError(f"unknown calculator module '{step_name}'")
        if self._redis is None:
            raise RuntimeError("redis is not initialized in worker process")

        timeout = int(self.worker_config.get("calc_acquire_timeout", 30))
        pack_id = self._redis.acquire_calc(step_name, timeout=timeout)
        if pack_id is None:
            raise TimeoutError(f"timed out acquiring calculator slot for '{step_name}'")
        try:
            module.acquire_pack_id(pack_id)
            module.prepare_runtime(sample.info)
            updated = module.execute(sample.info, runtime_prepared=True)
            lock = self._observables_lock
            if lock is not None:
                with lock:
                    self._merge_calculator_observables(sample, step_name, pack_id, updated)
            else:
                self._merge_calculator_observables(sample, step_name, pack_id, updated)
        finally:
            self._redis.release_calc(step_name, pack_id)

    def _merge_calculator_observables(
        self,
        sample: Sample,
        step_name: str,
        pack_id: str,
        updated: Mapping[str, Any],
    ) -> None:
        sample.observables.update(updated)
        if isinstance(sample.info, dict):
            pack_ids = dict(sample.info.get("pack_ids") or {})
            pack_ids[step_name] = pack_id
            sample.info["pack_ids"] = pack_ids
            sample.info["pack_id"] = pack_id
            sample.info["observables"] = dict(sample.observables)

    def _force_serial_layers(self) -> bool:
        """Rollback switch: run same-layer calculators one after another."""
        return bool(self.worker_config.get("force_serial_layers", False))

    def _run_calculator_steps(self, calc_steps: list[ExecutionStep], sample: Sample) -> None:
        if len(calc_steps) > 1 and not self._force_serial_layers():
            with concurrent.futures.ThreadPoolExecutor(max_workers=len(calc_steps)) as pool:
                futures = [
                    pool.submit(self._run_calculator_step, step.name, sample)
                    for step in calc_steps
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()
            return
        for step in calc_steps:
            self._run_calculator_step(step.name, sample)

    def _run_layer(self, layer: list[ExecutionStep], sample: Sample) -> None:
        """Run one execution-plan layer; fan out same-layer calculators concurrently."""
        calc_steps = [step for step in layer if step.type == "calculator"]
        inline_steps = [step for step in layer if step.type != "calculator"]
        self._run_calculator_steps(calc_steps, sample)

        for step in inline_steps:
            if step.type == "opera":
                self._run_opera_step(step.name, sample)
            elif step.type == "likelihood":
                self._run_likelihood(sample)

    def _run_opera_step(self, step_name: str, sample: Sample) -> None:
        module = self._operas.get(step_name)
        if module is None:
            raise KeyError(f"unknown opera module '{step_name}'")
        updated = module.execute(sample.observables, sample.info)
        sample.observables.update(updated)
        if isinstance(sample.info, dict):
            sample.info["observables"] = dict(sample.observables)

    def _run_likelihood(self, sample: Sample) -> None:
        if self._likelihood is None:
            return
        if not isinstance(sample.info, dict):
            sample.info = {}
        value = self._likelihood.calculate(sample.info)
        sample.observables = dict(sample.info.get("observables", sample.observables))
        sample._likelihood = value

    def _collect_cleanup_paths(self, sample: Sample) -> list[str]:
        if not isinstance(sample.info, dict):
            return []
        paths: list[str] = []
        for key in ("cleanup_paths", "staging_paths"):
            raw = sample.info.get(key)
            if isinstance(raw, str) and raw.strip():
                paths.append(raw.strip())
            elif isinstance(raw, list):
                paths.extend(str(item).strip() for item in raw if str(item).strip())
        return paths

    def _cleanup_transient_paths(self, sample: Sample) -> None:
        paths = self._collect_cleanup_paths(sample)
        if paths:
            delete_paths(paths, method=self._delete_method, missing_ok=True)

    def _handoff_sample_to_staging(self, sample: Sample) -> None:
        """Fast metadata-only move of materialized work dirs into staging (WP-D4.1)."""
        if not self._handoff_to_staging or not self._staging_dir:
            return
        if not isinstance(sample.info, dict):
            return
        save_dir = str(sample.info.get("save_dir") or "").strip()
        if not save_dir or not os.path.isdir(save_dir):
            return

        sample.close_logger()
        staging_path = stage_sample_dir(save_dir, self._staging_dir, sample.uuid)
        sample.info["staging_path"] = staging_path
        sample.info["product_list"] = list_product_names(staging_path)
        sample.info["save_dir"] = None

    def _stage_and_submit(self, sample: Sample) -> None:
        if self._redis is None:
            return
        self._redis.submit_result(sample.to_info_dict())

    def process_task(self, task: Mapping[str, Any]) -> None:
        """Core pipeline: rebuild Sample, execute workflow, submit result."""
        sample = Sample.from_task_dict(task)
        self._current_sample_uuid = sample.uuid
        top = get_jarvis_logger("worker").bind(
            worker_id=self.worker_id,
            sample_uuid=sample.uuid,
        )
        sample_config = dict(self.worker_config.get("sample_config") or {})
        sample.set_config(sample_config)
        sample.start()
        try:
            if self._mapper is not None:
                sample.bind_params(self._mapper)
            sample.materialize(worker_id=str(self.worker_id))
            layers = group_by_layer(sample.execution_plan)
            _profile_path = os.environ.get("JARVIS2_WORKER_PROFILE")
            _profile_started = time.monotonic()
            for layer in layers:
                self._run_layer(layer, sample)
            if _profile_path:
                with open(_profile_path, "a", encoding="utf-8") as _profile_handle:
                    _profile_handle.write(
                        json.dumps(
                            {
                                "worker_id": self.worker_id,
                                "layers": [
                                    {
                                        "calc_steps": [
                                            step.name
                                            for step in layer
                                            if step.type == "calculator"
                                        ],
                                        "width": len(layer),
                                    }
                                    for layer in layers
                                ],
                                "elapsed_sec": time.monotonic() - _profile_started,
                                "scheduler": (
                                    self._scheduler.snapshot()
                                    if self._scheduler is not None
                                    else {}
                                ),
                            }
                        )
                        + "\n"
                    )
            sample.status = "Completed"
            if isinstance(sample.info, dict):
                sample.info["status"] = "Completed"
        except Exception as exc:
            sample.status = "Failed"
            if isinstance(sample.info, dict):
                sample.info["status"] = "Failed"
            materialize_failure_artifacts(sample.info, error=exc)
            top.error("sample failed; see sample log -> %s", exc)
        finally:
            self._handoff_sample_to_staging(sample)
            self._stage_and_submit(sample)
            self._cleanup_transient_paths(sample)
            sample.close()
            self._current_sample_uuid = None

    def _main_loop(self) -> None:
        assert self._redis is not None
        pull_timeout = int(self.worker_config.get("pull_timeout", 5))
        while self._is_running:
            task = self._redis.pull_task(timeout=pull_timeout)
            if task is None:
                self._heartbeat("idle")
                continue
            self._heartbeat("busy")
            self.process_task(task)
            self._heartbeat("idle")

    def run(self) -> None:
        setup_jarvis_logging(role="worker")
        signal.signal(signal.SIGTERM, self._handle_signal)
        signal.signal(signal.SIGINT, self._handle_signal)
        self._init_redis()
        self._init_runtime()
        self._heartbeat("starting")
        try:
            self._main_loop()
        finally:
            self._heartbeat("stopped")
            if self._scheduler is not None:
                self._scheduler.shutdown(wait=True)
                self._scheduler = None
            if self._redis is not None:
                self._redis.close()
                self._redis = None


__all__ = ["Worker"]