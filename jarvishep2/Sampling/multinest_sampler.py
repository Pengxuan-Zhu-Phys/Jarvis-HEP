#!/usr/bin/env python3
"""MultiNest nested sampling on the V2 Redis runtime.

V1 ``Sampling.Method: MultiNest`` is **not** Fortran MultiNest / pymultinest —
it is a thin static wrapper around vendored dynesty ``NestedSampler`` (same UUID
channel + factory pool as Dynesty). V2 keeps that contract:

* engine: vendored ``NestedSampler`` only (never DynamicNestedSampler)
* transport: :class:`RedisEvaluationPool` + ``hep:feedback``
* diagnostics: ``DATABASE/multinest_result.csv`` (same column schema as dynesty)
* plot: Jarvis-PLOT ``dynesty_runplot`` with ``DataSet.name: dynesty`` and path
  pointing at the MultiNest CSV (V1 plot.py parity)
"""

from __future__ import annotations

import os
from collections.abc import Mapping
from typing import Any

from jarvishep2.Sampling.dynesty_sampler import (
    DynestySampler,
    export_dynesty_results_csv,
)
from jarvishep2.logging import get_jarvis_logger


class MultiNestSampler(DynestySampler):
    """YAML ``Sampling.Method: MultiNest`` — static nested sampling (V1 parity)."""

    method = "MultiNest"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.multinest")
        # V1 MultiNest always uses NestedSampler, never DynamicNestedSampler.
        self._use_dynamic = False
        self._multinest_csv_path: str | None = None

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        # Ignore Bounds.dynamic if a Dynesty card was copy-pasted — MultiNest is static.
        self._use_dynamic = False

    def multinest_result_csv_path(self, task_result_dir: str | None = None) -> str:
        """Canonical path: ``<task_result_dir>/DATABASE/multinest_result.csv``."""
        root = os.path.abspath(str(task_result_dir or self._task_result_dir()))
        return os.path.join(root, "DATABASE", "multinest_result.csv")

    def dynesty_result_csv_path(self, task_result_dir: str | None = None) -> str:
        # DynestySampler base calls this for the nested CSV; redirect to MultiNest path.
        return self.multinest_result_csv_path(task_result_dir)

    def save_dynesty_results_to_csv(
        self,
        output_csv: str | None = None,
    ) -> str | None:
        """Export nested results to V1 ``multinest_result.csv`` (Jarvis-PLOT schema)."""
        if self._sampler is None:
            self._logger.warning("save_multinest_results_to_csv: no sampler results")
            return None
        try:
            results = self._sampler.results
        except Exception as exc:
            self._logger.warning(
                "save_multinest_results_to_csv: cannot read results: %s", exc
            )
            return None
        path = output_csv or self.multinest_result_csv_path()
        written = export_dynesty_results_csv(
            results, path, fallback_nlive=self._nlive
        )
        self._dynesty_csv_path = written
        self._multinest_csv_path = written
        self._logger.info("MultiNest results saved → %s", written)
        try:
            summary_text = results.summary()
            if summary_text:
                self._logger.info("MultiNest summary → %s", summary_text)
        except Exception:
            pass
        return written

    # V1 name alias
    save_multinest_results_to_csv = save_dynesty_results_to_csv

    def run_adaptive(self, *, timeout: float = 3600.0) -> int:
        # Force static NestedSampler even if a parent left dynamic=True.
        self._use_dynamic = False
        ncall = super().run_adaptive(timeout=timeout)
        if self._summary is not None:
            self._summary["method"] = self.method
            csv_path = self._multinest_csv_path or self._dynesty_csv_path
            if csv_path:
                self._summary["multinest_result_csv"] = csv_path
                self._summary.pop("dynesty_result_csv", None)
        return ncall

    def export_runtime_state(self) -> dict[str, Any]:
        state = super().export_runtime_state()
        state["method"] = self.method
        state["use_dynamic"] = False
        state["multinest_csv_path"] = self._multinest_csv_path or self._dynesty_csv_path
        return state

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        super().import_runtime_state(state)
        self._use_dynamic = False
        path = state.get("multinest_csv_path") or state.get("dynesty_csv_path")
        self._multinest_csv_path = str(path) if path else None


def create_multinest() -> MultiNestSampler:
    return MultiNestSampler()


__all__ = [
    "MultiNestSampler",
    "create_multinest",
]
