#!/usr/bin/env python3
"""Random uniform sampler for Jarvis-HEP V2."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Mapping

import numpy as np

from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.sampling_utils import (
    BoolConversionError,
    evaluate_selection,
    map_u_to_physical,
)
from jarvishep2.Sampling.stateless_batch import (
    deterministic_sampler_uuid,
    run_stateless_distributed,
)
from jarvishep2.Sampling.variables import Variable, load_variables
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample


class RandomS(CheckpointedSampler):
    method = "Random"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.random")
        self.vars: list[Variable] = []
        self._index = 0
        self._accepted_index = 0
        self._maxp = 0
        self._dimensions = 0
        self._selectionexp: str | None = None
        self._seed = 0
        self._batch_size = 16
        self._uuid_by_accepted_index: dict[int, str] = {}
        self._u_by_uuid: dict[str, np.ndarray] = {}
        self._generator_ready = False

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._dimensions = len(self.vars)
        point_number = sampling.get("Point number", sampling.get("point_number"))
        if point_number is None:
            raise ValueError("Random sampler requires Sampling['Point number']")
        self._maxp = int(point_number)
        self._selectionexp = sampling.get("selection")
        self._seed = int(sampling.get("Seed", sampling.get("seed", 0)) or 0)
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        self._index = 0
        self._accepted_index = 0
        self._generator_ready = False

    def initialize(self) -> None:
        if self._seed:
            np.random.seed(self._seed)
        if self._selectionexp:
            probe = map_u_to_physical(np.random.rand(self._dimensions), self.vars)
            try:
                evaluate_selection(self._selectionexp, probe)
            except BoolConversionError as exc:
                raise ValueError(f"Invalid selection expression: {self._selectionexp}") from exc
        self._generator_ready = True

    def _ensure_ready(self) -> None:
        if not self._generator_ready:
            self.initialize()

    def propose_next(self) -> Sample | None:
        self._ensure_ready()
        while self._index < self._maxp:
            u_coords = np.random.random(self._dimensions).astype(np.float64)
            self._index += 1
            physical = map_u_to_physical(u_coords, self.vars)
            if self._selectionexp and not evaluate_selection(self._selectionexp, physical):
                continue
            accepted_index = self._accepted_index
            self._accepted_index += 1
            sample = self._build_sample(u_coords)
            sample.uuid = self._uuid_for_accepted_index(accepted_index)
            self._u_by_uuid[sample.uuid] = np.asarray(u_coords, dtype=np.float64)
            return sample
        return None

    def _uuid_for_accepted_index(self, accepted_index: int) -> str:
        if accepted_index in self._uuid_by_accepted_index:
            return self._uuid_by_accepted_index[accepted_index]
        uuid = deterministic_sampler_uuid(
            prefix="random",
            seed=self._seed,
            sample_index=accepted_index,
        )
        self._uuid_by_accepted_index[accepted_index] = uuid
        return uuid

    def repropose_unfinished(self) -> list[str]:
        if not self._repropose_after_resume:
            return []
        pending = [uuid for uuid in self._submitted_uuids if uuid not in self._completed_uuids]
        requeued: list[str] = []
        for uuid in pending:
            u_coords = self._u_by_uuid.get(uuid)
            if u_coords is None:
                continue
            sample = self._build_sample(u_coords)
            sample.uuid = uuid
            self._submit(sample)
            requeued.append(uuid)
        return requeued

    def run_distributed(self) -> int:
        return run_stateless_distributed(self, propose_next=self.propose_next)

    def at_safe_barrier(self) -> bool:
        if self._index < self._maxp:
            return False
        if not self._submitted_uuids:
            return True
        return set(self._submitted_uuids) <= self._completed_uuids

    def export_runtime_state(self) -> dict[str, Any]:
        return {
            "index": int(self._index),
            "maxp": int(self._maxp),
            "selectionexp": self._selectionexp,
            "seed": int(self._seed),
            "uuid_by_accepted_index": dict(self._uuid_by_accepted_index),
            "u_by_uuid": {key: value.tolist() for key, value in self._u_by_uuid.items()},
            "submitted_uuids": list(self._submitted_uuids),
            "completed_uuids": sorted(self._completed_uuids),
            "chains": [],
            "ready_queue": [],
            "control_state": self._checkpoint_control_state(),
            "numpy_random_state": np.random.get_state(),
        }

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        self._index = int(state.get("index", self._index) or 0)
        self._maxp = int(state.get("maxp", self._maxp) or self._maxp)
        self._selectionexp = state.get("selectionexp", self._selectionexp)
        self._seed = int(state.get("seed", self._seed) or self._seed)
        raw_uuid_map = state.get("uuid_by_accepted_index") or {}
        self._uuid_by_accepted_index = {int(k): str(v) for k, v in raw_uuid_map.items()}
        raw_u_by_uuid = state.get("u_by_uuid") or {}
        self._u_by_uuid = {
            str(key): np.asarray(value, dtype=np.float64) for key, value in raw_u_by_uuid.items()
        }
        np_state = state.get("numpy_random_state")
        if np_state is not None:
            np.random.set_state(np_state)
        self._import_checkpoint_control_state(state)


__all__ = ["RandomS"]