#!/usr/bin/env python3
"""Cartesian grid sampler for Jarvis-HEP V2."""

from __future__ import annotations

import itertools
import time
from typing import Any, Mapping

import numpy as np

from jarvishep2.sampling.fixed_set_sampler import FixedSetSampler
from jarvishep2.sampling.sampling_utils import evaluate_selection, physical_from_u
from jarvishep2.log_kv import PermilleProgress
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
    _checkpoint_excluded_attributes = frozenset({"_logger", "_P", "_submit_progress"})

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.grid")
        self._P: np.ndarray | None = None
        self._submit_progress: PermilleProgress | None = None

    def initialize(self, *, reset: bool = True) -> None:
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
        if reset:
            self._index = 0
            self._accepted_index = 0
        if reset or self._submit_progress is None:
            self._submit_progress = PermilleProgress(
                self._logger,
                total=int(self.info["NSamples"]),
                label="grid points sampled",
            )
            if self._index == 0 and not bool(
                getattr(self, "_suppress_submit_progress", False)
            ):
                self._submit_progress.update(0, force=True)
        self._logger.info(
            "Grid generated %d points in %.2f s",
            self.info["NSamples"],
            self.info["t0"],
        )

    def _ensure_grid(self) -> None:
        if self._P is None:
            self.initialize(reset=False)

    def propose_next(self) -> Sample | None:
        self._ensure_grid()
        assert self._P is not None
        while self._index < len(self._P):
            row = self._P[self._index]
            self._index += 1
            self._emit_progress()
            physical = physical_from_u(row, self.vars, self._mapper_pipeline)
            if self._selectionexp and not evaluate_selection(
                self._selectionexp,
                physical,
                context=self._expression_context,
            ):
                continue
            accepted_index = self._accepted_index
            self._accepted_index += 1
            u_coords = np.asarray(row, dtype=np.float64)
            sample = self._build_sample(u_coords)
            sample.sample_index = accepted_index
            sample.uuid = self._uuid_for_accepted_index(accepted_index)
            return sample
        return None

    def _stream_exhausted(self) -> bool:
        return self._P is not None and self._index >= len(self._P)

    def _emit_progress(self) -> None:
        """Emit one progress heartbeat for each newly traversed grid permille."""
        if bool(getattr(self, "_suppress_submit_progress", False)):
            return
        total = int(self.info.get("NSamples", 0) or 0)
        if self._P is not None:
            total = max(total, int(len(self._P)))
        if total <= 0:
            return
        if self._submit_progress is None or self._submit_progress.total != total:
            self._submit_progress = PermilleProgress(
                self._logger,
                total=total,
                label="grid points sampled",
            )
        self._submit_progress.update(self._index)

    def export_runtime_state(self) -> dict[str, Any]:
        payload = self._common_export_fields()
        payload.update({"numpy_random_state": None})
        return payload

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        self._import_common_fields(state)
        self._submit_progress = None


__all__ = ["Grid", "grid_sampling"]
