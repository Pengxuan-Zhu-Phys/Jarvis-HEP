#!/usr/bin/env python3
"""Shared base for finite, stateless batch samplers (D9.7)."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.stateless_batch import deterministic_sampler_uuid
from jarvishep2.Sampling.variables import Variable, load_variables
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample


class FixedSetSampler(CheckpointedSampler):
    """Template base: index/accepted maps + batch flush + common set_config fields.

    Subclasses implement generation via ``propose_next`` and any sampler-specific
    initialize/export state. Shared batch submission lives here so modules no
    longer poke ``_submit`` private members from the outside.
    """

    uuid_prefix: str = "fixed"

    def __init__(self) -> None:
        super().__init__()
        self.vars: list[Variable] = []
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

    def _uuid_for_accepted_index(self, accepted_index: int) -> str:
        if accepted_index in self._uuid_by_accepted_index:
            return self._uuid_by_accepted_index[accepted_index]
        uuid = deterministic_sampler_uuid(
            prefix=self.uuid_prefix,
            seed=self._seed,
            sample_index=accepted_index,
        )
        self._uuid_by_accepted_index[accepted_index] = uuid
        return uuid

    def flush_batch(self, batch: list[Sample]) -> int:
        if not batch:
            return 0
        if len(batch) == 1:
            self._submit(batch[0])
        else:
            self._submit_group(batch)
        self._submitted_uuids.extend(sample.uuid for sample in batch)
        return len(batch)

    def run_distributed(self) -> int:
        pushed = 0
        batch: list[Sample] = []
        while True:
            sample = self.propose_next()
            if sample is None:
                break
            batch.append(sample)
            if len(batch) >= self._batch_size:
                pushed += self.flush_batch(batch)
                batch = []
        if batch:
            pushed += self.flush_batch(batch)
        return pushed

    def propose_next(self) -> Sample | None:
        raise NotImplementedError

    def _stream_exhausted(self) -> bool:
        """Return True when the sampler has no more candidates to propose."""
        raise NotImplementedError

    def at_safe_barrier(self) -> bool:
        if not self._stream_exhausted():
            return False
        if not self._submitted_uuids:
            return True
        return set(self._submitted_uuids) <= self._completed_uuids

    def _common_export_fields(self) -> dict[str, Any]:
        return {
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
        }

    def _import_common_fields(self, state: Mapping[str, Any]) -> None:
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


__all__ = ["FixedSetSampler"]
