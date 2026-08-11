#!/usr/bin/env python3
"""Random uniform sampler for Jarvis-HEP V2."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from jarvishep2.Sampling.fixed_set_sampler import FixedSetSampler
from jarvishep2.Sampling.sampling_utils import (
    BoolConversionError,
    evaluate_selection,
    physical_from_u,
)
from jarvishep2.log_kv import PermilleProgress
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.sample import Sample


class RandomS(FixedSetSampler):
    method = "Random"
    uuid_prefix = "random"
    _checkpoint_saved_attributes = frozenset({"_maxp"})
    _checkpoint_excluded_attributes = frozenset(
        {"_logger", "_dimensions", "_generator_ready", "_submit_progress"}
    )

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.random")
        self._maxp = 0
        self._dimensions = 0
        self._generator_ready = False
        self._submit_progress: PermilleProgress | None = None

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        bounds = sampling.get("Bounds") if isinstance(sampling.get("Bounds"), Mapping) else {}
        self._dimensions = len(self.vars)
        point_number = bounds.get("point_number")
        if point_number is None:
            raise ValueError("Random sampler requires Sampling.Bounds.point_number")
        self._maxp = int(point_number)
        self._index = 0
        self._accepted_index = 0
        self._generator_ready = False
        self._submit_progress = None

    def initialize(self) -> None:
        # Seed 0 is valid and must reseed for reproducible resume trajectories.
        np.random.seed(int(self._seed))
        if self._selectionexp:
            probe = physical_from_u(
                np.random.rand(self._dimensions),
                self.vars,
                self._mapper_pipeline,
            )
            try:
                evaluate_selection(
                    self._selectionexp,
                    probe,
                    context=self._expression_context,
                )
            except BoolConversionError as exc:
                raise ValueError(f"Invalid selection expression: {self._selectionexp}") from exc
        self._generator_ready = True
        self._submit_progress = PermilleProgress(
            self._logger,
            total=self._maxp,
            label="random candidates sampled",
        )
        if not bool(getattr(self, "_suppress_submit_progress", False)):
            self._submit_progress.update(0, force=True)

    def _ensure_ready(self) -> None:
        if not self._generator_ready:
            self.initialize()

    def propose_next(self) -> Sample | None:
        self._ensure_ready()
        while self._index < self._maxp:
            u_coords = np.random.random(self._dimensions).astype(np.float64)
            self._index += 1
            self._emit_progress()
            physical = physical_from_u(u_coords, self.vars, self._mapper_pipeline)
            if self._selectionexp and not evaluate_selection(
                self._selectionexp,
                physical,
                context=self._expression_context,
            ):
                continue
            accepted_index = self._accepted_index
            self._accepted_index += 1
            sample = self._build_sample(u_coords)
            sample.sample_index = accepted_index
            sample.uuid = self._uuid_for_accepted_index(accepted_index)
            return sample
        return None

    def _emit_progress(self) -> None:
        """Log every new 1‰ of the configured draws; warn at whole percents."""
        # DATABASE prefix replay is local bookkeeping, not fresh sampling work.
        if bool(getattr(self, "_suppress_submit_progress", False)):
            return
        if self._submit_progress is None:
            self._submit_progress = PermilleProgress(
                self._logger,
                total=self._maxp,
                label="random candidates sampled",
            )
        self._submit_progress.update(self._index)

    def _stream_exhausted(self) -> bool:
        return self._index >= self._maxp

    def export_runtime_state(self) -> dict[str, Any]:
        # The initial checkpoint must capture the seeded generator state, not
        # whichever global NumPy state happened to precede the first proposal.
        self._ensure_ready()
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
            self._generator_ready = True
        # Logging state is process-local. The first new proposal reports the
        # restored cursor without adding mutable progress data to checkpoints.
        self._submit_progress = None


__all__ = ["RandomS"]
