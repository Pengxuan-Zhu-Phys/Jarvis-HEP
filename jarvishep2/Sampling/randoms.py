#!/usr/bin/env python3
"""Random uniform sampler for Jarvis-HEP V2."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from jarvishep2.Sampling.fixed_set_sampler import FixedSetSampler
from jarvishep2.Sampling.sampling_utils import (
    BoolConversionError,
    evaluate_selection,
    map_u_to_physical,
)
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.sample import Sample


class RandomS(FixedSetSampler):
    method = "Random"
    uuid_prefix = "random"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.random")
        self._maxp = 0
        self._dimensions = 0
        self._generator_ready = False

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        self._dimensions = len(self.vars)
        point_number = sampling.get("Point number", sampling.get("point_number"))
        if point_number is None:
            raise ValueError("Random sampler requires Sampling['Point number']")
        self._maxp = int(point_number)
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

    def _stream_exhausted(self) -> bool:
        return self._index >= self._maxp

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

    def export_runtime_state(self) -> dict[str, Any]:
        payload = self._common_export_fields()
        payload.update(
            {
                "maxp": int(self._maxp),
                "numpy_random_state": np.random.get_state(),
            }
        )
        return payload

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        self._import_common_fields(state)
        self._maxp = int(state.get("maxp", self._maxp) or self._maxp)
        np_state = state.get("numpy_random_state")
        if np_state is not None:
            np.random.set_state(np_state)


__all__ = ["RandomS"]
