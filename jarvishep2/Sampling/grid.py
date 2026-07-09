#!/usr/bin/env python3
"""Cartesian grid sampler for Jarvis-HEP V2."""

from __future__ import annotations

import itertools
import time
from typing import Any, Mapping

import numpy as np

from jarvishep2.Sampling.fixed_set_sampler import FixedSetSampler
from jarvishep2.Sampling.sampling_utils import evaluate_selection, map_u_to_physical
from jarvishep2.logging import get_jarvis_logger
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


class Grid(FixedSetSampler):
    method = "Grid"
    uuid_prefix = "grid"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.grid")
        self._P: np.ndarray | None = None

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

    def _stream_exhausted(self) -> bool:
        return self._P is not None and self._index >= len(self._P)

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

    def export_runtime_state(self) -> dict[str, Any]:
        payload = self._common_export_fields()
        payload.update({"grid_points": self._P, "numpy_random_state": None})
        return payload

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        grid = state.get("grid_points")
        if grid is not None:
            self._P = np.asarray(grid)
        self._import_common_fields(state)
        if self._P is not None:
            self.info["NSamples"] = int(self._P.shape[0])


__all__ = ["Grid", "grid_sampling"]
