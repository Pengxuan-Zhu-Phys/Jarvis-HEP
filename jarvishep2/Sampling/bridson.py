#!/usr/bin/env python3
"""Bridson (Poisson-disk / blue-noise) sampler for Jarvis-HEP V2."""

from __future__ import annotations

import hashlib
import time
from copy import deepcopy
from typing import Any, Mapping

import numpy as np
from scipy.special import gammainc

from jarvishep2.Sampling.fixed_set_sampler import FixedSetSampler
from jarvishep2.Sampling.sampling_utils import (
    evaluate_selection,
    map_row_to_physical,
    row_to_u_coords,
)
from jarvishep2.log_kv import format_duration
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


class Bridson(FixedSetSampler):
    """Stateless batch proposer: precompute Poisson-disk set, submit via Redis."""

    method = "Bridson"
    uuid_prefix = "bridson"
    _checkpoint_saved_attributes = frozenset({"barinfo"})
    _checkpoint_excluded_attributes = frozenset(
        {"_logger", "_P", "_radius", "_k", "_max_inflight"}
    )

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.bridson")
        self._P: np.ndarray | None = None
        self._radius = 0.0
        self._k = 30
        self._max_inflight = 1
        self.barinfo: dict[str, Any] = {}

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        # FixedSetSampler already set _batch_size from get_runtime_block (default 256).
        # Do not re-read raw Runtime here — missing batch_size would wrongly fall back
        # to MaxWorker/workers and change Redis pipeline sizing vs Random/Grid.
        sampling = dict(self.config.get("Sampling") or {})
        bounds = sampling.get("Bounds") if isinstance(sampling.get("Bounds"), Mapping) else {}
        if not isinstance(bounds, Mapping) or not bounds:
            raise ValueError(
                "Bridson requires Sampling.Bounds with Radius and MaxAttempt"
            )
        runtime = get_runtime_block(self.config)
        workers = int(runtime.get("workers", 1) or 1)
        if "Radius" not in bounds:
            raise ValueError("Bridson requires Sampling.Bounds.Radius")
        if "MaxAttempt" not in bounds:
            raise ValueError("Bridson requires Sampling.Bounds.MaxAttempt")
        self._radius = float(bounds["Radius"])
        self._k = int(bounds["MaxAttempt"])
        # V1-like backpressure: at most Bounds.MaxWorker (default = Runtime.workers).
        self._max_inflight = max(1, int(bounds.get("MaxWorker", workers) or workers))

    def initialize(self, *, reset: bool = True) -> None:
        ndim = len(self.vars)
        if ndim < 2 or ndim >= 5:
            raise ValueError("Bridson supports 2d to 4d parameter spaces only")
        self._logger.warning("Initializing the Bridson Sampling")
        # Seed 0 is a valid deterministic seed — do not treat it as "unset"
        # (``if self._seed`` would skip seeding and break resume: same uuid
        # stream, different coordinates).
        np.random.seed(int(self._seed))
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
        if reset:
            self._index = 0
            self._accepted_index = 0
            # Fresh submit-progress bar for this grid generation.
            self.barinfo = {}
        self._logger.info(
            "Bridson Sampler obtains %d samples in %.2f sec",
            self.info["NSamples"],
            self.info["t0"],
        )

    def _ensure_grid(self) -> None:
        if self._P is None:
            self.initialize(reset=False)

    def _uuid_for_accepted_index(self, accepted_index: int) -> str:
        # Keep historical Bridson UUID formula (prefix "bridson:").
        return deterministic_bridson_uuid(seed=self._seed, sample_index=accepted_index)

    def propose_next(self) -> Sample | None:
        self._ensure_grid()
        assert self._P is not None
        while self._index < len(self._P):
            row = self._P[self._index]
            self._index += 1
            # V1 progress heartbeat: advance on every grid step (even if selection rejects).
            self._emit_progress()
            physical = map_row_to_physical(row, self.vars)
            if self._selectionexp and not evaluate_selection(
                self._selectionexp,
                physical,
                context=self._expression_context,
            ):
                continue
            accepted_index = self._accepted_index
            self._accepted_index += 1
            sample = self._build_sample(row_to_u_coords(row, self.vars))
            sample.sample_index = accepted_index
            sample.uuid = self._uuid_for_accepted_index(accepted_index)
            return sample
        return None

    def _stream_exhausted(self) -> bool:
        return self._P is not None and self._index >= len(self._P)

    def _emit_progress(self) -> None:
        """V1 Bridson submit heartbeat: one log line per new ‰ of the grid."""
        # Local DATABASE prefix replay must not look like real submissions.
        if bool(getattr(self, "_suppress_submit_progress", False)):
            return
        total = int(self.info.get("NSamples", 0) or 0)
        if self._P is not None:
            total = max(total, int(len(self._P)))
        if total <= 0:
            return

        def _log_progress(permille: int) -> None:
            # Keep V1 wording (including the historical "submited" spelling).
            msg = "{}‰ of {}/{} samples submited in {}".format(
                permille,
                int(self._index),
                total,
                format_duration(time.time() - float(self.barinfo["t0"])),
            )
            # INFO for 1‰ steps; WARNING only for exact 1%, 2%, … milestones.
            if permille > 0 and permille % 10 == 0:
                self._logger.warning(msg)
            else:
                self._logger.info(msg)

        if not self.barinfo:
            self.barinfo = {
                "total": total,
                "t0": time.time(),
                "permille": 0,
            }
            _log_progress(0)
            return

        self.barinfo.setdefault("total", total)
        self.barinfo.setdefault("t0", time.time())
        permille = int(self._index / float(self.barinfo["total"]) * 1000)
        if permille != int(self.barinfo.get("permille", -1)):
            self.barinfo["permille"] = permille
            _log_progress(permille)

    def export_runtime_state(self) -> dict[str, Any]:
        payload = self._common_export_fields()
        payload.update(
            {
                "barinfo": deepcopy(self.barinfo),
                "numpy_random_state": None,
            }
        )
        return payload

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        self._import_common_fields(state)
        self.barinfo = deepcopy(state.get("barinfo", self.barinfo))


__all__ = [
    "Bridson",
    "Bridson_sampling",
    "deterministic_bridson_uuid",
    "hypersphere_surface_sample",
    "hypersphere_volume_sample",
]
