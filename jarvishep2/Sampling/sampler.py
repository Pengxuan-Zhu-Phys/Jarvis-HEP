#!/usr/bin/env python3
"""Sampler base with Redis submission path for Jarvis-HEP V2."""

from __future__ import annotations

from typing import Any, Mapping
from uuid import uuid4

import numpy as np

from jarvishep2.redis_queue import RedisQueue
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample
from jarvishep2.workflow import execution_plan_template


class SamplingVirtial:
    """Minimal sampler base: proposes Samples and submits to Redis in V2."""

    def __init__(self) -> None:
        self.config: dict[str, Any] = {}
        self.info: dict[str, Any] = {}
        self.redis: RedisQueue | None = None
        self.runtime_mode = "auto"
        self._execution_plan_template: list[dict[str, Any]] = []
        self._opera_modules: list[dict[str, Any]] = []
        self._calculator_modules: list[dict[str, Any]] = []
        self._sample_artifacts = "auto"

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        self.config = dict(config_info)
        runtime = get_runtime_block(self.config)
        self.runtime_mode = str(runtime.get("mode", "auto"))
        self._sample_artifacts = str(runtime.get("sample_artifacts", "auto"))

    def set_redis(self, redis: RedisQueue) -> None:
        self.redis = redis

    def set_execution_plan_template(
        self,
        opera_modules: list[dict[str, Any]] | None = None,
        *,
        calculator_modules: list[dict[str, Any]] | None = None,
        include_likelihood: bool = True,
    ) -> None:
        self._opera_modules = list(opera_modules or [])
        self._calculator_modules = list(calculator_modules or [])
        self._execution_plan_template = execution_plan_template(
            opera_modules=self._opera_modules,
            calculator_modules=self._calculator_modules,
            include_likelihood=include_likelihood,
        )

    def _build_sample(self, u_coords: np.ndarray | list[float] | None = None) -> Sample:
        coords = np.asarray(u_coords if u_coords is not None else [], dtype=np.float64)
        sample = Sample(
            uuid=str(uuid4()),
            u_coords=coords,
            execution_plan=[],
            sample_artifacts=self._sample_artifacts,
        )
        if self._execution_plan_template:
            from jarvishep2.sample import ExecutionStep

            sample.execution_plan = [
                ExecutionStep.from_dict(step) for step in self._execution_plan_template
            ]
        return sample

    def _submit(self, sample: Sample) -> None:
        if self.runtime_mode != "redis":
            raise RuntimeError("redis submission requires Runtime.mode == 'redis'")
        if self.redis is None:
            raise RuntimeError("sampler redis queue is not configured")
        self.redis.push_task(sample.to_task_dict())

    def _submit_group(self, samples: list[Sample]) -> None:
        if self.runtime_mode != "redis":
            raise RuntimeError("redis submission requires Runtime.mode == 'redis'")
        if self.redis is None:
            raise RuntimeError("sampler redis queue is not configured")
        self.redis.push_many_tasks([sample.to_task_dict() for sample in samples])


__all__ = ["SamplingVirtial"]