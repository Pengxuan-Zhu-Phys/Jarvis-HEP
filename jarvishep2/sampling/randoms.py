#!/usr/bin/env python3
"""Random uniform sampler for Jarvis-HEP V2."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from jarvishep2.sampling.completion_sampler import CompletionSampler
from jarvishep2.sampling.sampling_utils import (
    BoolConversionError,
    evaluate_selection,
    physical_from_u,
)
from jarvishep2.log_kv import PermilleProgress
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.sample import Sample


class RandomS(CompletionSampler):
    method = "Random"
    uuid_prefix = "random"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.random")
        self._maxp = 0
        self._dimensions = 0
        self._generator_ready = False
        self._rng = np.random.default_rng(0)
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
        self._rng = np.random.default_rng(self._seed)
        self._submit_progress = None

    def log_settings_rows(self) -> list[tuple[str, Any]]:
        return [("point number", self._maxp), *super().log_settings_rows()]

    def initialize(self) -> None:
        self._rng = np.random.default_rng(self._seed)
        if self._selectionexp:
            probe = physical_from_u(
                self._rng.random(self._dimensions),
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
            u_coords = self._rng.random(self._dimensions)
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
            return self._remember(sample)
        return None

    def _emit_progress(self) -> None:
        """Log every new 1‰ of the configured draws; warn at whole percents."""
        if self._submit_progress is None:
            self._submit_progress = PermilleProgress(
                self._logger,
                total=self._maxp,
                label="random candidates sampled",
            )
        self._submit_progress.update(self._index)

    def _stream_exhausted(self) -> bool:
        return self._index >= self._maxp

    def _source_state(self) -> dict[str, Any]:
        return {"index": self._index, "rng": self._rng.bit_generator.state,
                "initialized": self._generator_ready}

    def _restore_source_state(self, state: Mapping[str, Any]) -> None:
        index = state["index"]
        if type(index) is not int or not 0 <= index <= self._maxp:
            raise ValueError("Random checkpoint cursor is outside point_number")
        self._index = index
        self._rng.bit_generator.state = state["rng"]
        if type(state["initialized"]) is not bool or (index and not state["initialized"]):
            raise ValueError("invalid Random generator initialization state")
        self._generator_ready = state["initialized"]
        self._submit_progress = None


__all__ = ["RandomS"]
