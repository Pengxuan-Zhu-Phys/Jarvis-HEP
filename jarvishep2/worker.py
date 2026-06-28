#!/usr/bin/env python3
"""Redis-backed Worker process for Jarvis-HEP V2."""

from __future__ import annotations

import signal
import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.likelihood import LogLikelihoodEvaluator
from jarvishep2.Module.calculator import CalculatorModule, mint_pack_id
from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging
from jarvishep2.mapper import build_mapper
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.operas import preload_operas
from jarvishep2.redis_queue import RedisQueue
from jarvishep2.sample import Sample, materialize_failure_artifacts
from jarvishep2.workflow import group_by_layer

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
        pack_id = mint_pack_id()
        if isinstance(sample.info, dict):
            sample.info["pack_id"] = pack_id
        module.acquire_pack_id(pack_id)
        updated = module.execute(sample.info)
        sample.observables.update(updated)
        if isinstance(sample.info, dict):
            sample.info["observables"] = dict(sample.observables)

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

    def _stage_and_submit(self, sample: Sample) -> None:
        if self._redis is None:
            return
        info = sample.to_info_dict()
        self._redis.submit_result(info)
        self._redis.incr_op("sample")

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
            for layer in group_by_layer(sample.execution_plan):
                for step in layer:
                    if step.type == "calculator":
                        self._run_calculator_step(step.name, sample)
                    elif step.type == "opera":
                        self._run_opera_step(step.name, sample)
                    elif step.type == "likelihood":
                        self._run_likelihood(sample)
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
            self._stage_and_submit(sample)
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
            if self._redis is not None:
                self._redis.close()
                self._redis = None


__all__ = ["Worker"]