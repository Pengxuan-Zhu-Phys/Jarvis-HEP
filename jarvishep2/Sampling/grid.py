#!/usr/bin/env python3
"""Cartesian grid sampler for Jarvis-HEP V2."""

from __future__ import annotations

import itertools
import time
from copy import deepcopy
from typing import Any, Mapping

import numpy as np

from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.sampling_utils import evaluate_selection, map_u_to_physical
from jarvishep2.Sampling.stateless_batch import (
    deterministic_sampler_uuid,
    run_stateless_distributed,
)
from jarvishep2.Sampling.variables import Variable, load_variables
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample


def grid_sampling(dimensions: np.ndarray) -> np.ndarray:
    eps = np.finfo(np.float64).eps
    grid_ranges: list[np.ndarray] = []
    for steps in dimensions:
        nsteps = int(steps)
        if nsteps <= 1:
            grid_ranges.append(np.array([0.5], dtype=np.float64))
            continue
        vals = np.linspace(0.0, 1.0, nsteps, dtype=np.float64)
        grid_ranges.append(np.clip(vals, eps, 1.0 - eps))
    return np.array(list(itertools.product(*grid_ranges)), dtype=np.float64)


class Grid(CheckpointedSampler):
    method = "Grid"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.grid")
        self.vars: list[Variable] = []
        self._P: np.ndarray | None = None
        self._index = 0
        self._accepted_index = 0
        self._selectionexp: str | None = None
        self._seed = 0
        self._batch_size = 16
        self._uuid_by_accepted_index: dict[int, str] = {}
        self._u_by_uuid: dict[str, np.ndarray] = {}

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._selectionexp = sampling.get("selection")
        self._seed = int(sampling.get("Seed", sampling.get("seed", 0)) or 0)
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))

    def initialize(self) -> None:
        t0 = time.time()
        dims = []
        for var in self.vars:
            if "num" not in var.parameters:
                raise ValueError(
                    f"Grid sampler requires distribution.parameters.num for variable '{var.name}'"
                )
            dims.append(int(var.parameters["num"]))
        dims_arr = np.array(dims, dtype=np.int64)
        self._P = grid_sampling(dims_arr)
        self.info["NSamples"] = int(self._P.shape[0])
        self.info["t0"] = time.time() - t0
        self._index = 0
        self._accepted_index = 0
        self._logger.info(
            "Grid generated %d points in %.2f s",
            self.info["NSamples"],
            self.info["t0"],
        )

    def _ensure_grid(self) -> None:
        if self._P is None:
            self.initialize()

    def propose_next(self) -> Sample | None:
        self._ensure_grid()
        assert self._P is not None
        while self._index < len(self._P):
            row = self._P[self._index]
            self._index += 1
            physical = map_u_to_physical(row, self.vars)
            if self._selectionexp and not evaluate_selection(self._selectionexp, physical):
                continue
            accepted_index = self._accepted_index
            self._accepted_index += 1
            u_coords = np.asarray(row, dtype=np.float64)
            sample = self._build_sample(u_coords)
            sample.uuid = self._uuid_for_accepted_index(accepted_index)
            self._u_by_uuid[sample.uuid] = u_coords
            return sample
        return None

    def _uuid_for_accepted_index(self, accepted_index: int) -> str:
        if accepted_index in self._uuid_by_accepted_index:
            return self._uuid_by_accepted_index[accepted_index]
        uuid = deterministic_sampler_uuid(
            prefix="grid",
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
                accepted_index = next(
                    (idx for idx, mapped in self._uuid_by_accepted_index.items() if mapped == uuid),
                    None,
                )
                if accepted_index is not None:
                    row = self._find_row_for_accepted_index(accepted_index)
                    if row is not None:
                        u_coords = np.asarray(row, dtype=np.float64)
            if u_coords is None:
                continue
            sample = self._build_sample(u_coords)
            sample.uuid = uuid
            self._submit(sample)
            requeued.append(uuid)
        return requeued

    def _find_row_for_accepted_index(self, target: int) -> np.ndarray | None:
        if self._P is None:
            return None
        accepted = 0
        for scan_index in range(len(self._P)):
            row = self._P[scan_index]
            physical = map_u_to_physical(row, self.vars)
            if self._selectionexp and not evaluate_selection(self._selectionexp, physical):
                continue
            if accepted == target:
                return np.asarray(row, dtype=np.float64)
            accepted += 1
        return None

    def run_distributed(self) -> int:
        return run_stateless_distributed(self, propose_next=self.propose_next)

    def at_safe_barrier(self) -> bool:
        if self._P is None or self._index < len(self._P):
            return False
        if not self._submitted_uuids:
            return True
        return set(self._submitted_uuids) <= self._completed_uuids

    def export_runtime_state(self) -> dict[str, Any]:
        return {
            "grid_points": self._P,
            "index": int(self._index),
            "accepted_index": int(self._accepted_index),
            "selectionexp": self._selectionexp,
            "seed": int(self._seed),
            "uuid_by_accepted_index": dict(self._uuid_by_accepted_index),
            "u_by_uuid": {key: value.tolist() for key, value in self._u_by_uuid.items()},
            "submitted_uuids": list(self._submitted_uuids),
            "completed_uuids": sorted(self._completed_uuids),
            "chains": [],
            "ready_queue": [],
            "control_state": self._checkpoint_control_state(),
            "numpy_random_state": None,
        }

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        grid = state.get("grid_points")
        if grid is not None:
            self._P = np.asarray(grid)
        self._index = int(state.get("index", self._index) or 0)
        self._accepted_index = int(state.get("accepted_index", self._accepted_index) or 0)
        self._selectionexp = state.get("selectionexp", self._selectionexp)
        self._seed = int(state.get("seed", self._seed) or self._seed)
        raw_uuid_map = state.get("uuid_by_accepted_index") or {}
        self._uuid_by_accepted_index = {int(k): str(v) for k, v in raw_uuid_map.items()}
        raw_u_by_uuid = state.get("u_by_uuid") or {}
        self._u_by_uuid = {
            str(key): np.asarray(value, dtype=np.float64) for key, value in raw_u_by_uuid.items()
        }
        self._import_checkpoint_control_state(state)
        if self._P is not None:
            self.info["NSamples"] = int(self._P.shape[0])


__all__ = ["Grid", "grid_sampling"]