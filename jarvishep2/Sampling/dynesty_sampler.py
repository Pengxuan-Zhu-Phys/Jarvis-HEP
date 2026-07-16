#!/usr/bin/env python3
"""Dynesty nested sampling on the V2 Redis runtime (D13.5, route A).

Uses the **vendored** dynesty package under
``jarvishep2.Sampling.Source.Dynesty`` (stock 3.0.0 + Jarvis UUID channel),
with likelihood evaluations via :class:`RedisEvaluationPool`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any
from uuid import uuid4

import numpy as np

from jarvishep2.Sampling.feedback_sampler import FeedbackSampler
from jarvishep2.Sampling.redis_evaluation_pool import RedisEvaluationPool
from jarvishep2.Sampling.sampling_utils import evaluate_selection, map_u_to_physical
from jarvishep2.Sampling.variables import load_variables
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample


def _jarvis_prior_transform(u: np.ndarray) -> np.ndarray:
    """V1-style prior_transform: append a uuid to the unit-cube vector.

    Physical mapping is done in the loglikelihood bridge (Worker Sample path
    still maps u→x via the UMapper). Dynesty stores coords in ``live_v`` and
    uuid in ``live_uid`` after the vendored split.
    """
    u = np.asarray(u, dtype=float).reshape(-1)
    out = np.empty(u.size + 1, dtype=object)
    out[:-1] = u
    out[-1] = str(uuid4())
    return out


class DynestySampler(FeedbackSampler):
    """YAML ``Sampling.Method: Dynesty``."""

    method = "Dynesty"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.dynesty")
        self.vars: list = []
        self._dim = 0
        self._nlive = 100
        self._dlogz = 0.5
        self._rseed = 0
        self._run_nested_kwargs: dict[str, Any] = {}
        self._selectionexp: str | None = None
        self._sampler = None  # dynesty Nested/Dynamic sampler
        self._finished = False
        self._summary: dict[str, Any] | None = None
        self._use_dynamic = True

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._dim = len(self.vars)
        if self._dim < 1:
            raise ValueError("Dynesty requires at least one Sampling.Variable")
        bounds = sampling.get("Bounds") or {}
        if not isinstance(bounds, Mapping):
            bounds = {}
        self._nlive = max(2, int(bounds.get("nlive", bounds.get("n_live", 100)) or 100))
        self._dlogz = float(bounds.get("dlogz", 0.5) or 0.5)
        self._rseed = int(
            bounds.get("rseed", sampling.get("Seed", sampling.get("seed", 0))) or 0
        )
        self._init_seed_sequence(self._rseed)
        self._selectionexp = sampling.get("selection")
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        run_nested = bounds.get("run_nested") or {}
        if isinstance(run_nested, Mapping):
            self._run_nested_kwargs = dict(run_nested)
        else:
            self._run_nested_kwargs = {}
        self._run_nested_kwargs.setdefault("dlogz", self._dlogz)
        # print_progress can be True in YAML; keep as bool
        dynamic = bounds.get("dynamic", bounds.get("Dynamic", True))
        self._use_dynamic = bool(dynamic)

    def propose_generation(self) -> Sequence[Sample] | None:
        # Dynesty drives proposals via its own state machine + pool.map.
        return None

    def absorb_generation(self, results: Sequence[Mapping[str, Any]]) -> None:
        return None

    def _build_sample_for_pool(self, payload: np.ndarray, uuid: str) -> Sample:
        """Build a Sample from unit-cube coords (uuid already assigned)."""
        u = np.asarray(payload, dtype=np.float64).reshape(-1)
        if u.size != self._dim:
            # If physical coords slipped through, clip to dim
            u = u[: self._dim]
        sample = self._build_sample(u)
        sample.uuid = uuid
        if self._selectionexp:
            physical = map_u_to_physical(u, self.vars)
            if not evaluate_selection(
                self._selectionexp, physical, context=self._expression_context
            ):
                # Mark for -inf via a side channel — Worker still runs; bridge uses selection.
                sample.observables["_dynesty_selection_reject"] = True
        return sample

    def _local_loglike(self, params: Any) -> float:
        """Synchronous path when pool is unavailable (should not be used in production)."""
        arr = np.asarray(params, dtype=object).reshape(-1)
        if arr.size >= 2 and isinstance(arr[-1], (str, np.str_)):
            u = np.asarray(arr[:-1], dtype=float)
        else:
            u = np.asarray(arr, dtype=float)
        # Without Workers we cannot evaluate real logL — return a toy Gaussian.
        return float(-0.5 * np.sum((u - 0.5) ** 2))

    def run_adaptive(self, *, timeout: float = 3600.0) -> int:
        self._require_redis("Dynesty.run_adaptive")
        self._ensure_seed_sequence()
        if self._finished and self._sampler is not None:
            return 0

        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import (
            DynamicNestedSampler,
            NestedSampler,
        )
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty.jarvis_logging import (
            bind_inner,
            install_warnings_bridge,
            set_dynesty_logger,
        )

        # All dynesty progress/warnings use this sampler's Jarvis logger.
        inner_logger = bind_inner(self._logger)
        set_dynesty_logger(inner_logger)
        install_warnings_bridge()

        pool = RedisEvaluationPool(
            self.redis,
            build_sample=self._build_sample_for_pool,
            batch_size=self._batch_size,
            timeout=timeout,
            seed=self._seed,
            method=self.method,
            logger=inner_logger,
        )

        rstate = np.random.default_rng(self._rseed)
        sampler_cls = DynamicNestedSampler if self._use_dynamic else NestedSampler
        self._logger.info(
            "Starting %s nlive=%d ndim=%d dynamic=%s",
            sampler_cls.__name__,
            self._nlive,
            self._dim,
            self._use_dynamic,
        )
        # ndim = unit-cube dim; prior_transform appends uuid (vendored split).
        self._sampler = sampler_cls(
            loglikelihood=self._local_loglike,
            prior_transform=_jarvis_prior_transform,
            ndim=self._dim,
            nlive=self._nlive,
            pool=pool,
            rstate=rstate,
            queue_size=max(1, self._batch_size),
        )
        # Attach logger onto nested sampler instance (used by run_nested progress).
        try:
            self._sampler.logger = inner_logger
            # DynamicNestedSampler wraps an inner sampler object
            for attr in ("sampler", "base_sampler"):
                inner = getattr(self._sampler, attr, None)
                if inner is not None and hasattr(inner, "logger"):
                    inner.logger = inner_logger
                elif inner is not None:
                    setattr(inner, "logger", inner_logger)
        except Exception:
            pass
        # Dynesty LogLikelihood uses pool.map(loglikelihood, pars) — our pool
        # evaluates remotely and returns floats; wrap is still required by API.
        run_kwargs = dict(self._run_nested_kwargs)
        # Ensure dlogz present
        run_kwargs.setdefault("dlogz", self._dlogz)
        # Progress always on → routed to Jarvis logger (not bare stderr).
        run_kwargs["print_progress"] = True

        self._sampler.run_nested(**run_kwargs)
        self._finished = True
        self._summary = self._build_summary()
        self.checkpoint_at_barrier(reason="dynesty_finished")
        self._logger.info(
            "Dynesty finished logZ=%.4f ± %.4f niter=%s",
            self._summary.get("logz", float("nan")),
            self._summary.get("logzerr", float("nan")),
            self._summary.get("niter"),
        )
        # Approximate submitted = ncall
        return int(self._summary.get("ncall") or 0)

    def _build_summary(self) -> dict[str, Any]:
        if self._sampler is None:
            return {"method": self.method, "finished": self._finished}
        try:
            res = self._sampler.results
            logz = float(np.asarray(res["logz"])[-1]) if len(res["logz"]) else float("nan")
            logzerr = (
                float(np.asarray(res["logzerr"])[-1]) if len(res["logzerr"]) else float("nan")
            )
            ncall = int(np.sum(res["ncall"])) if "ncall" in res.keys() else 0
            # ncall may already be scalar
            try:
                ncall = int(res["ncall"]) if np.ndim(res["ncall"]) == 0 else ncall
            except Exception:
                pass
            niter = int(res["niter"]) if "niter" in res.keys() else 0
            samples_uid = None
            if "samples_uid" in res.keys():
                samples_uid = list(res["samples_uid"])
            return {
                "method": self.method,
                "finished": True,
                "logz": logz,
                "logzerr": logzerr,
                "ncall": ncall,
                "niter": niter,
                "nlive": self._nlive,
                "samples_uid": samples_uid,
            }
        except Exception as exc:
            self._logger.warning("failed to build dynesty summary: %s", exc)
            return {"method": self.method, "finished": self._finished, "error": str(exc)}

    def summary(self) -> dict[str, Any]:
        if self._summary is None:
            self._summary = self._build_summary()
        return dict(self._summary)

    def export_runtime_state(self) -> dict[str, Any]:
        state = self._feedback_export_state()
        state.update(
            {
                "method": self.method,
                "finished": self._finished,
                "summary": self._summary,
                "nlive": self._nlive,
                "dlogz": self._dlogz,
                "rseed": self._rseed,
                "call_index": getattr(
                    getattr(self, "_pool_ref", None), "_call_index", 0
                ),
            }
        )
        # Full dynesty pickle resume is D13.6-progressive; store summary + flags.
        return state

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        self._feedback_import_state(state)
        self._finished = bool(state.get("finished", False))
        self._summary = state.get("summary")
        self._nlive = int(state.get("nlive", self._nlive) or self._nlive)
        self._dlogz = float(state.get("dlogz", self._dlogz) or self._dlogz)
        self._rseed = int(state.get("rseed", self._rseed) or self._rseed)


def create_dynesty() -> DynestySampler:
    return DynestySampler()


__all__ = ["DynestySampler", "create_dynesty", "_jarvis_prior_transform"]
