#!/usr/bin/env python3
from __future__ import annotations

import os
import threading
from copy import deepcopy
from typing import Any, Dict, Sequence

import numpy as np

from jarvishep.Sampling.bucketallocator import BucketAllocator
from jarvishep.sample import Sample


class NestedLikelihoodBridge:
    """Pickle-friendly runtime bridge for dynesty-style nested samplers."""

    def __init__(
        self,
        *,
        sampler_name: str,
        variables: Sequence[Any],
        base_sample_cfg: Dict[str, Any],
        sample_cls=Sample,
        bucket_state: Dict[str, Any] | None = None,
        bucket_limit: int = 200,
        bucket_width: int = 6,
        submit_limit: int = 1,
    ) -> None:
        self.sampler_name = str(sampler_name)
        self.variables = tuple(variables or ())
        self.base_sample_cfg = deepcopy(base_sample_cfg or {})
        self.sample_cls = sample_cls or Sample
        self.bucket_state = deepcopy(bucket_state or {})
        self.bucket_limit = max(1, int(bucket_limit))
        self.bucket_width = max(1, int(bucket_width))
        self.submit_limit = max(1, int(submit_limit))
        self._factory = None
        self._logger = None
        self._bucket_alloc: BucketAllocator | None = None
        self._submit_gate = threading.BoundedSemaphore(value=self.submit_limit)

    def attach_runtime(self, *, factory=None, logger=None) -> None:
        if factory is not None:
            self._factory = factory
        if logger is not None:
            self._logger = logger

    def _log(self, level: str, message: str) -> None:
        logger = self._logger
        if logger is None:
            return
        method = getattr(logger, str(level).lower(), None)
        if callable(method):
            try:
                method(message)
            except Exception:
                pass

    def _ensure_bucket_allocator(self) -> BucketAllocator:
        if self._bucket_alloc is not None:
            return self._bucket_alloc

        base_path = self.bucket_state.get("base_path")
        if not isinstance(base_path, str) or not base_path:
            base_path = self.base_sample_cfg.get("sample_dirs")
        if not isinstance(base_path, str) or not base_path:
            task_root = self.base_sample_cfg.get("task_result_dir", os.getcwd())
            base_path = os.path.join(str(task_root), "SAMPLE")

        self._bucket_alloc = BucketAllocator(
            base_path=base_path,
            limit=self.bucket_limit,
            width=self.bucket_width,
            start_bucket=1,
            on_bucket_sealed=None,
        )
        if self.bucket_state:
            self._bucket_alloc.set_state(self.bucket_state)
        else:
            self._bucket_alloc.check_and_update()
        return self._bucket_alloc

    def _next_bucket_dir(self) -> str | None:
        bucket_alloc = self._ensure_bucket_allocator()
        return bucket_alloc.next_bucket_dir()

    def _mark_sample_finished(self, sample_info: Dict[str, Any] | None) -> None:
        if not isinstance(sample_info, dict):
            return
        save_dir = sample_info.get("save_dir")
        if not save_dir:
            return
        bucket_dir = os.path.dirname(str(save_dir).rstrip(os.sep))
        bucket_alloc = self._ensure_bucket_allocator()
        try:
            bucket_alloc.mark_sample_finished(bucket_dir)
        except Exception:
            pass

    def _build_sample_config(self, save_dir: str | None) -> Dict[str, Any]:
        cfg = deepcopy(self.base_sample_cfg)
        if save_dir is not None:
            cfg["save_dir"] = save_dir
        nuisance = cfg.get("nuisance")
        if isinstance(nuisance, dict):
            cfg["nuisance"] = deepcopy(nuisance)
        return cfg

    def __call__(self, params):
        if self._factory is None:
            raise RuntimeError(f"{self.sampler_name} runtime bridge has no attached factory")

        arr = np.asarray(params, dtype=object)
        if arr.size == 0:
            raise ValueError(f"{self.sampler_name} received empty parameter vector")

        unit = np.asarray(arr[:-1], dtype=np.float64)
        uuid = arr[-1]
        values = {}
        for ii, variable in enumerate(self.variables):
            values[getattr(variable, "name", f"x{ii}")] = variable.map_standard_random_to_distribution(unit[ii])

        sample = self.sample_cls(values)
        sample.update_uuid(str(uuid))
        sample_cfg = self._build_sample_config(self._next_bucket_dir())
        sample.set_config(sample_cfg)

        try:
            if not self._submit_gate.acquire(blocking=False):
                self._submit_gate.acquire()
            try:
                future = self._factory.submit_task(sample.info)
                result = future.result()
            finally:
                self._submit_gate.release()
            return result
        except Exception as exc:
            self._log("error", f"{self.sampler_name} likelihood bridge failure -> {exc}")
            raise
        finally:
            self._mark_sample_finished(sample.info)
            try:
                sample.close()
            except Exception:
                pass

    def __getstate__(self) -> Dict[str, Any]:
        state = dict(self.__dict__)
        bucket_alloc = state.pop("_bucket_alloc", None)
        if bucket_alloc is not None:
            state["bucket_state"] = bucket_alloc.get_state()
        state["_factory"] = None
        state["_logger"] = None
        state["_submit_gate"] = None
        return state

    def __setstate__(self, state: Dict[str, Any]) -> None:
        self.__dict__.update(state)
        self._factory = None
        self._logger = None
        self._bucket_alloc = None
        self._submit_gate = threading.BoundedSemaphore(value=self.submit_limit)
