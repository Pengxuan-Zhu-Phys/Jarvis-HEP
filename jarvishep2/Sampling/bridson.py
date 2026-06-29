#!/usr/bin/env python3
"""Bridson (Poisson-disk / blue-noise) sampler for Jarvis-HEP V2."""

from __future__ import annotations

import hashlib
import time
from copy import deepcopy
from typing import Any, Mapping

import numpy as np
from scipy.special import gammainc

from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.sampling_utils import (
    evaluate_selection,
    map_row_to_physical,
    row_to_u_coords,
)
from jarvishep2.Sampling.variables import load_variables
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample


def hypersphere_volume_sample(center: np.ndarray, radius: float, k: int = 1) -> np.ndarray:
    ndim = center.size
    x = np.random.normal(size=(k, ndim))
    ssq = np.sum(x**2, axis=1)
    fr = radius * gammainc(ndim / 2, ssq / 2) ** (1 / ndim) / np.sqrt(ssq)
    frtiled = np.tile(fr.reshape(k, 1), (1, ndim))
    return center + np.multiply(x, frtiled)


def hypersphere_surface_sample(center: np.ndarray, radius: float, k: int = 1) -> np.ndarray:
    ndim = center.size
    vec = np.random.standard_normal(size=(k, ndim))
    vec /= np.linalg.norm(vec, axis=1)[:, None]
    return center + np.multiply(vec, radius)


def squared_distance(p0: np.ndarray, p1: np.ndarray) -> float:
    return float(np.sum(np.square(p0 - p1)))


def Bridson_sampling(
    *,
    dims: np.ndarray,
    radius: float,
    k: int = 30,
    hypersphere_sample=hypersphere_volume_sample,
) -> np.ndarray:
    """Robert Bridson 2007 Poisson-disk sampling in arbitrary dimensions."""
    ndim = dims.size
    sample_factor = 2
    if hypersphere_sample == hypersphere_volume_sample:
        sample_factor = 2
    if hypersphere_sample == hypersphere_surface_sample:
        eps = 0.001
        sample_factor = 1 + eps

    def in_limits(p: np.ndarray) -> bool:
        return bool(np.all(np.zeros(ndim) <= p) and np.all(p < dims))

    cellsize = radius / np.sqrt(ndim)
    gridsize = np.ceil(dims / cellsize).astype(int)
    squared_radius = radius * radius
    p_grid = np.empty(np.append(gridsize, ndim), dtype=np.float32)
    p_grid.fill(np.nan)
    points: list[np.ndarray] = []

    def in_neighborhood(p: np.ndarray, n: int = 2) -> bool:
        indices = (p / cellsize).astype(int)
        indmin = np.maximum(indices - n, np.zeros(ndim, dtype=int))
        indmax = np.minimum(indices + n + 1, gridsize)
        if not np.isnan(p_grid[tuple(indices)][0]):
            return True
        slices = tuple(slice(indmin[i], indmax[i]) for i in range(ndim))
        with np.errstate(invalid="ignore"):
            if np.any(np.sum(np.square(p - p_grid[slices]), axis=-1) < squared_radius):
                return True
        return False

    def add_point(p: np.ndarray) -> None:
        points.append(p)
        indices = (p / cellsize).astype(int)
        p_grid[tuple(indices)] = p

    add_point(np.random.uniform(np.zeros(ndim), dims))
    while points:
        index = np.random.randint(len(points))
        p = points[index]
        del points[index]
        candidates = hypersphere_sample(np.array(p), radius * sample_factor, k)
        for q in candidates:
            if in_limits(q) and not in_neighborhood(q):
                add_point(q)
    return p_grid[~np.isnan(p_grid).any(axis=-1)]


def deterministic_bridson_uuid(*, seed: int, sample_index: int) -> str:
    digest = hashlib.sha256(f"bridson:{int(seed)}:{int(sample_index)}".encode("utf-8")).hexdigest()
    return (
        f"{digest[0:8]}-{digest[8:12]}-{digest[12:16]}-{digest[16:20]}-{digest[20:32]}"
    )


class Bridson(CheckpointedSampler):
    """Stateless batch proposer: precompute Poisson-disk set, submit via Redis."""

    method = "Bridson"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.bridson")
        self.vars: list = []
        self._P: np.ndarray | None = None
        self._index = 0
        self._accepted_index = 0
        self._radius = 0.0
        self._k = 30
        self._selectionexp: str | None = None
        self._seed = 0
        self._batch_size = 16
        self._max_inflight = 1
        self.barinfo: dict[str, Any] = {}
        self._uuid_by_accepted_index: dict[int, str] = {}

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._radius = float(sampling["Radius"])
        self._k = int(sampling["MaxAttempt"])
        self._selectionexp = sampling.get("selection")
        self._seed = int(sampling.get("Seed", sampling.get("seed", 0)) or 0)
        workers = int(runtime.get("workers", 1) or 1)
        self._max_inflight = max(1, int(sampling.get("MaxWorker", workers) or workers))
        self._batch_size = max(1, int(runtime.get("batch_size", self._max_inflight) or self._max_inflight))

    def initialize(self) -> None:
        ndim = len(self.vars)
        if ndim < 2 or ndim >= 5:
            raise ValueError("Bridson supports 2d to 4d parameter spaces only")
        if self._seed:
            np.random.seed(self._seed)
        t0 = time.time()
        dims = np.array([float(var.parameters["length"]) for var in self.vars], dtype=np.float64)
        self._P = Bridson_sampling(
            dims=dims,
            radius=self._radius,
            k=self._k,
            hypersphere_sample=hypersphere_surface_sample,
        )
        self.info["NSamples"] = int(self._P.shape[0])
        self.info["t0"] = time.time() - t0
        self._index = 0
        self._accepted_index = 0
        self._logger.info(
            "Bridson generated %d accepted grid cells in %.2f s",
            self.info["NSamples"],
            self.info["t0"],
        )

    def _ensure_grid(self) -> None:
        if self._P is None:
            self.initialize()

    def _uuid_for_accepted_index(self, accepted_index: int) -> str:
        if accepted_index in self._uuid_by_accepted_index:
            return self._uuid_by_accepted_index[accepted_index]
        uuid = deterministic_bridson_uuid(seed=self._seed, sample_index=accepted_index)
        self._uuid_by_accepted_index[accepted_index] = uuid
        return uuid

    def propose_next(self) -> Sample | None:
        self._ensure_grid()
        assert self._P is not None
        while self._index < len(self._P):
            row = self._P[self._index]
            self._index += 1
            physical = map_row_to_physical(row, self.vars)
            if self._selectionexp and not evaluate_selection(self._selectionexp, physical):
                continue
            accepted_index = self._accepted_index
            self._accepted_index += 1
            sample = self._build_sample(row_to_u_coords(row, self.vars))
            sample.uuid = self._uuid_for_accepted_index(accepted_index)
            self._emit_progress()
            return sample
        return None

    def _emit_progress(self) -> None:
        if not self.barinfo:
            self.barinfo = {"total": int(self.info.get("NSamples", 0)), "permille": 0}
        if self.barinfo["total"] <= 0:
            return
        permille = int(self._index / self.barinfo["total"] * 1000)
        if permille != self.barinfo.get("permille"):
            self.barinfo["permille"] = permille

    def repropose_unfinished(self) -> list[str]:
        if not self._repropose_after_resume:
            return []
        pending = [uuid for uuid in self._submitted_uuids if uuid not in self._completed_uuids]
        if not pending or self._P is None:
            return []
        requeued: list[str] = []
        index_by_uuid = {uuid: idx for idx, uuid in self._uuid_by_accepted_index.items()}
        for uuid in pending:
            accepted_index = index_by_uuid.get(uuid)
            if accepted_index is None:
                continue
            row = self._find_row_for_accepted_index(accepted_index)
            if row is None:
                continue
            sample = self._build_sample(row_to_u_coords(row, self.vars))
            sample.uuid = uuid
            self._submit(sample)
            requeued.append(uuid)
        return requeued

    def _find_row_for_accepted_index(self, target: int) -> np.ndarray | None:
        assert self._P is not None
        accepted = 0
        scan_index = 0
        while scan_index < len(self._P):
            row = self._P[scan_index]
            scan_index += 1
            physical = map_row_to_physical(row, self.vars)
            if self._selectionexp and not evaluate_selection(self._selectionexp, physical):
                continue
            if accepted == target:
                return row
            accepted += 1
        return None

    def run_distributed(self) -> int:
        """Submit all accepted Bridson points in Redis pipeline batches."""
        self._ensure_grid()
        pushed = 0
        batch: list[Sample] = []
        while True:
            sample = self.propose_next()
            if sample is None:
                break
            batch.append(sample)
            if len(batch) >= self._batch_size:
                pushed += self._flush_batch(batch)
                batch = []
        if batch:
            pushed += self._flush_batch(batch)
        return pushed

    def _flush_batch(self, batch: list[Sample]) -> int:
        if not batch:
            return 0
        if len(batch) == 1:
            self._submit(batch[0])
        else:
            self._submit_group(batch)
        self._submitted_uuids.extend(sample.uuid for sample in batch)
        return len(batch)

    def at_safe_barrier(self) -> bool:
        if self._P is None:
            return False
        if self._index < len(self._P):
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
            "barinfo": deepcopy(self.barinfo),
            "uuid_by_accepted_index": dict(self._uuid_by_accepted_index),
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
        self.barinfo = deepcopy(state.get("barinfo", self.barinfo))
        raw_uuid_map = state.get("uuid_by_accepted_index") or {}
        self._uuid_by_accepted_index = {int(k): str(v) for k, v in raw_uuid_map.items()}
        self._import_checkpoint_control_state(state)
        if self._P is not None:
            self.info["NSamples"] = int(self._P.shape[0])


__all__ = [
    "Bridson",
    "Bridson_sampling",
    "deterministic_bridson_uuid",
    "hypersphere_surface_sample",
    "hypersphere_volume_sample",
]