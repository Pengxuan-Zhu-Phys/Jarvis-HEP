#!/usr/bin/env python3
"""Dynesty nested sampling on the V2 Redis runtime (D13.5, route A).

Uses the **vendored** dynesty package under
``jarvishep2.Sampling.Source.Dynesty`` (stock 3.1.0 + Jarvis UUID channel),
with likelihood evaluations via :class:`RedisEvaluationPool`.

YAML ``Sampling.Bounds`` is aligned with the **official dynesty 3.x API**:

* constructor kwargs → ``NestedSampler`` / ``DynamicNestedSampler``
  (via ``Bounds.sampler`` and/or flat Bounds keys)
* run kwargs → ``run_nested`` (via ``Bounds.run_nested``)
* ``dlogz`` is mapped to static ``dlogz`` or dynamic ``dlogz_init``

Post-run diagnostics follow the V1 Jarvis-PLOT path: write
``DATABASE/dynesty_result.csv`` (column schema for
``dynesty_runplot``), then let ``plot_scene`` emit a stock jplot YAML.
Rendering stays in Jarvis-PLOT — HEP only maps results → CSV → scene.
"""

from __future__ import annotations

import csv
import inspect
import os
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


# Official NestedSampler / DynamicNestedSampler user-facing constructor keys
# (excluding HEP-injected: loglikelihood, prior_transform, ndim, pool, rstate).
NESTED_CONSTRUCTOR_USER_KEYS: frozenset[str] = frozenset(
    {
        "nlive",
        "bound",
        "sample",
        "periodic",
        "reflective",
        "update_interval",
        "first_update",
        "queue_size",
        "use_pool",
        "live_points",
        "logl_args",
        "logl_kwargs",
        "ptform_args",
        "ptform_kwargs",
        "enlarge",
        "bootstrap",
        "walks",
        "facc",
        "slices",
        "ncdim",
        "blob",
        "save_evaluation_history",
        "history_filename",
    }
)

# NestedSampler.run_nested (static)
STATIC_RUN_NESTED_KEYS: frozenset[str] = frozenset(
    {
        "maxiter",
        "maxcall",
        "dlogz",
        "logl_max",
        "add_live",
        "print_progress",
        "print_func",
        "save_bounds",
        "checkpoint_file",
        "checkpoint_every",
        "resume",
    }
)

# DynamicNestedSampler.run_nested
DYNAMIC_RUN_NESTED_KEYS: frozenset[str] = frozenset(
    {
        "nlive_init",
        "maxiter_init",
        "maxcall_init",
        "dlogz_init",
        "logl_max_init",
        "nlive_batch",
        "wt_function",
        "wt_kwargs",
        "maxiter_batch",
        "maxcall_batch",
        "maxiter",
        "maxcall",
        "maxbatch",
        "n_effective",
        "stop_function",
        "stop_kwargs",
        "use_stop",
        "save_bounds",
        "print_progress",
        "print_func",
        "live_points",
        "resume",
        "checkpoint_file",
        "checkpoint_every",
    }
)

# Bounds keys owned by HEP / meta (never forwarded as constructor kwargs).
_BOUNDS_META_KEYS: frozenset[str] = frozenset(
    {
        "nlive",
        "n_live",
        "rseed",
        "seed",
        "Seed",
        "dynamic",
        "Dynamic",
        "dlogz",
        "dlogz_init",
        "run_nested",
        "sampler",
        "Sampler",
        "constructor",
        "Constructor",
    }
)

# Always injected by HEP — strip if a user pastes them into sampler: block.
_HEP_INJECTED_CONSTRUCTOR_KEYS: frozenset[str] = frozenset(
    {
        "loglikelihood",
        "prior_transform",
        "ndim",
        "pool",
        "rstate",
    }
)


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


def _filter_known_kwargs(
    raw: Mapping[str, Any],
    allowed: frozenset[str],
    *,
    context: str,
    logger: Any | None = None,
) -> dict[str, Any]:
    """Keep only keys in *allowed*; warn on unknown (no crash)."""
    out: dict[str, Any] = {}
    for key, value in raw.items():
        name = str(key)
        if name in allowed:
            out[name] = value
        elif logger is not None:
            logger.warning(
                "nested sampler: ignoring unknown %s key %r (not in dynesty API)",
                context,
                name,
            )
    return out


def extract_nested_constructor_kwargs(
    bounds: Mapping[str, Any],
    *,
    logger: Any | None = None,
) -> dict[str, Any]:
    """Build NestedSampler/DynamicNestedSampler kwargs from Bounds.

    Sources (later wins):

    1. Flat ``Bounds.<constructor_key>`` (official names: bound, sample, walks, …)
    2. Nested ``Bounds.sampler`` / ``Bounds.constructor`` block (preferred)

    HEP-injected keys (loglikelihood, pool, …) are stripped. ``nlive`` is left
    in so callers can still override via sampler block if desired.
    """
    raw: dict[str, Any] = {}
    for key, value in bounds.items():
        name = str(key)
        if name in _BOUNDS_META_KEYS:
            continue
        if name in NESTED_CONSTRUCTOR_USER_KEYS:
            raw[name] = value
    for block_key in ("sampler", "Sampler", "constructor", "Constructor"):
        block = bounds.get(block_key)
        if isinstance(block, Mapping):
            for key, value in block.items():
                raw[str(key)] = value
    for banned in _HEP_INJECTED_CONSTRUCTOR_KEYS:
        raw.pop(banned, None)
    return _filter_known_kwargs(
        raw, NESTED_CONSTRUCTOR_USER_KEYS, context="constructor", logger=logger
    )


def extract_run_nested_kwargs(
    bounds: Mapping[str, Any],
    *,
    dynamic: bool,
    dlogz_default: float = 0.5,
    logger: Any | None = None,
) -> dict[str, Any]:
    """Build ``run_nested`` kwargs with static/dynamic evidence-threshold aliases.

    * Static: uses ``dlogz`` (from ``run_nested.dlogz`` or ``Bounds.dlogz``).
    * Dynamic: uses ``dlogz_init``; a user ``dlogz`` is remapped to ``dlogz_init``.
    """
    run_raw = bounds.get("run_nested") or {}
    if not isinstance(run_raw, Mapping):
        run_raw = {}
    run: dict[str, Any] = {str(k): v for k, v in run_raw.items()}

    if dynamic:
        # Official Dynamic API: dlogz_init (not dlogz).
        if "dlogz" in run and "dlogz_init" not in run:
            run["dlogz_init"] = run.pop("dlogz")
            if logger is not None:
                logger.info(
                    "nested sampler: mapped run_nested.dlogz → dlogz_init (Dynamic)"
                )
        elif "dlogz" in run:
            run.pop("dlogz")
            if logger is not None:
                logger.warning(
                    "nested sampler: dropped run_nested.dlogz "
                    "(Dynamic uses dlogz_init; already set)"
                )
        if "dlogz_init" not in run:
            if bounds.get("dlogz_init") is not None:
                run["dlogz_init"] = bounds["dlogz_init"]
            elif bounds.get("dlogz") is not None:
                run["dlogz_init"] = bounds["dlogz"]
            else:
                run["dlogz_init"] = dlogz_default
        allowed = DYNAMIC_RUN_NESTED_KEYS
        context = "run_nested(Dynamic)"
    else:
        # Official NestedSampler API: dlogz.
        if "dlogz_init" in run and "dlogz" not in run:
            run["dlogz"] = run.pop("dlogz_init")
            if logger is not None:
                logger.info(
                    "nested sampler: mapped run_nested.dlogz_init → dlogz (static)"
                )
        elif "dlogz_init" in run:
            run.pop("dlogz_init")
        if "dlogz" not in run:
            if bounds.get("dlogz") is not None:
                run["dlogz"] = bounds["dlogz"]
            else:
                run["dlogz"] = dlogz_default
        allowed = STATIC_RUN_NESTED_KEYS
        context = "run_nested(static)"

    return _filter_known_kwargs(run, allowed, context=context, logger=logger)


def _signature_param_names(func: Any) -> frozenset[str]:
    try:
        return frozenset(
            name
            for name in inspect.signature(func).parameters
            if name not in {"self", "cls"}
        )
    except (TypeError, ValueError):
        return frozenset()


def filter_kwargs_to_callable(
    kwargs: Mapping[str, Any],
    func: Any,
    *,
    context: str,
    logger: Any | None = None,
) -> dict[str, Any]:
    """Final safety filter against the live callable signature."""
    allowed = _signature_param_names(func)
    if not allowed:
        return dict(kwargs)
    return _filter_known_kwargs(kwargs, allowed, context=context, logger=logger)


def _results_get(results: Any, key: str, default: Any = None) -> Any:
    """Read a Results key whether Results supports ``.keys()`` / ``[]``."""
    try:
        if hasattr(results, "keys") and key in results.keys():
            return results[key]
    except Exception:
        pass
    try:
        return results[key]
    except Exception:
        return getattr(results, key, default)


def _samples_nlive_array(results: Any, n_rows: int, fallback_nlive: int) -> np.ndarray:
    """Per-sample live-point count (V1 ``samples_nlive`` / dynesty ``samples_n``)."""
    raw = _results_get(results, "samples_n", None)
    if raw is not None:
        arr = np.asarray(raw, dtype=float).reshape(-1)
        if arr.size == n_rows:
            return arr
        if arr.size > 0:
            # Pad / trim defensively
            out = np.empty(n_rows, dtype=float)
            n = min(arr.size, n_rows)
            out[:n] = arr[:n]
            out[n:] = arr[-1] if arr.size else float(fallback_nlive)
            return out
    try:
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty.utils import (
            _get_nsamps_samples_n,
        )

        _nsamps, samples_n = _get_nsamps_samples_n(results)
        arr = np.asarray(samples_n, dtype=float).reshape(-1)
        if arr.size == n_rows:
            return arr
    except Exception:
        pass
    nlive = int(_results_get(results, "nlive", fallback_nlive) or fallback_nlive)
    return np.full(n_rows, float(nlive), dtype=float)


def export_dynesty_results_csv(
    results: Any,
    output_csv: str,
    *,
    fallback_nlive: int = 100,
    param_names: Sequence[str] | None = None,
) -> str:
    """Write clean nested-result CSV (Dynesty dead points only).

    This is **not** the full likelihood-call archive (``samples.csv`` /
    ``samples.hdf5``). Rows equal ``niter`` nested samples with weights,
    evidence trajectory, uuids, and parameter coordinates.

    Column schema matches V1 ``Sampling.dynesty.save_dynesty_results_to_csv`` and
    Jarvis-PLOT ``dynesty_runplot`` (``samples_v[i]`` / ``samples_u[i]``).
    When *param_names* is provided, physical parameter columns with those names
    are added alongside ``samples_v[i]`` for readability.
    """
    logl = np.asarray(_results_get(results, "logl", []), dtype=float).reshape(-1)
    n_rows = int(logl.size)
    if n_rows < 1:
        raise ValueError("dynesty results have no logl rows to export")

    def _col(key: str, *, dtype=float, default: float | str = 0.0) -> np.ndarray:
        raw = _results_get(results, key, None)
        if raw is None:
            if dtype is str or dtype is object:
                return np.array([""] * n_rows, dtype=object)
            return np.full(n_rows, default, dtype=float)
        arr = np.asarray(raw)
        if arr.ndim == 0:
            return np.full(n_rows, arr.item(), dtype=object if dtype is object else float)
        flat = arr.reshape(-1) if arr.ndim == 1 else arr
        if flat.ndim == 1:
            out = np.empty(n_rows, dtype=object if dtype is object else float)
            n = min(flat.size, n_rows)
            if dtype is object:
                out[:n] = [str(x) if x is not None else "" for x in flat[:n]]
                out[n:] = ""
            else:
                out[:n] = np.asarray(flat[:n], dtype=float)
                out[n:] = float(flat[n - 1]) if n else default
            return out
        return flat

    logwt = _col("logwt")
    logvol = _col("logvol")
    logz = _col("logz")
    logzerr = _col("logzerr")
    samples_it = _col("samples_it")
    samples_id = _col("samples_id")
    information = _col("information")
    ncall = _col("ncall")
    samples_nlive = _samples_nlive_array(results, n_rows, fallback_nlive)

    uuid_raw = _results_get(results, "samples_uid", None)
    if uuid_raw is None:
        uuids = [""] * n_rows
    else:
        uuids = [str(x) if x is not None else "" for x in np.asarray(uuid_raw).reshape(-1)]
        if len(uuids) < n_rows:
            uuids = uuids + [""] * (n_rows - len(uuids))
        elif len(uuids) > n_rows:
            uuids = uuids[:n_rows]

    samples = np.asarray(_results_get(results, "samples", np.zeros((n_rows, 0))), dtype=float)
    if samples.ndim == 1:
        samples = samples.reshape(n_rows, -1) if samples.size else np.zeros((n_rows, 0))
    if samples.shape[0] != n_rows and samples.size:
        # Best-effort reshape when trailing live points differ
        samples = samples.reshape(-1, samples.shape[-1])[:n_rows]
    samples_u = np.asarray(_results_get(results, "samples_u", np.zeros((n_rows, 0))), dtype=float)
    if samples_u.ndim == 1:
        samples_u = samples_u.reshape(n_rows, -1) if samples_u.size else np.zeros((n_rows, 0))
    if samples_u.shape[0] != n_rows and samples_u.size:
        samples_u = samples_u.reshape(-1, samples_u.shape[-1])[:n_rows]

    ndim_v = int(samples.shape[1]) if samples.ndim == 2 else 0
    ndim_u = int(samples_u.shape[1]) if samples_u.ndim == 2 else 0

    names: list[str] = []
    if param_names:
        for ii, raw_name in enumerate(param_names):
            if ii >= ndim_v:
                break
            text = str(raw_name or "").strip()
            if text and text not in {
                "uuid",
                "log_weight",
                "log_Like",
                "log_PriorVolume",
                "log_Evidence",
                "log_Evidence_err",
                "samples_nlive",
                "ncall",
                "samples_it",
                "samples_id",
                "information",
            }:
                names.append(text)
            else:
                names.append("")
        while len(names) < ndim_v:
            names.append("")

    fieldnames = [
        "uuid",
        "log_weight",
        "log_Like",
        "log_PriorVolume",
        "log_Evidence",
        "log_Evidence_err",
        "samples_nlive",
        "ncall",
        "samples_it",
        "samples_id",
        "information",
    ]
    # Named physical coords first (clean read), then V1 samples_v/u aliases.
    for name in names:
        if name and name not in fieldnames:
            fieldnames.append(name)
    for ii in range(ndim_v):
        fieldnames.append(f"samples_v[{ii}]")
    for ii in range(ndim_u):
        fieldnames.append(f"samples_u[{ii}]")

    rows: list[dict[str, Any]] = []
    for i in range(n_rows):
        row: dict[str, Any] = {
            "uuid": uuids[i],
            "log_weight": float(logwt[i]),
            "log_Like": float(logl[i]),
            "log_PriorVolume": float(logvol[i]),
            "log_Evidence": float(logz[i]),
            "log_Evidence_err": float(logzerr[i]),
            "samples_nlive": float(samples_nlive[i]),
            "ncall": float(ncall[i]),
            "samples_it": float(samples_it[i]),
            "samples_id": float(samples_id[i]),
            "information": float(information[i]),
        }
        for ii in range(ndim_v):
            val = float(samples[i, ii])
            row[f"samples_v[{ii}]"] = val
            if ii < len(names) and names[ii]:
                row[names[ii]] = val
        for ii in range(ndim_u):
            row[f"samples_u[{ii}]"] = float(samples_u[i, ii])
        rows.append(row)

    output = os.path.abspath(str(output_csv))
    parent = os.path.dirname(output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output


class DynestySampler(FeedbackSampler):
    """YAML ``Sampling.Method: Dynesty``.

    Always uses **DynamicNestedSampler**. Static nested sampling is
    ``Sampling.Method: MultiNest`` only — there is no ``Bounds.dynamic`` switch.
    Constructor + run_nested kwargs follow the official dynesty 3.x surface.
    """

    method = "Dynesty"
    # MultiNest subclasses force static NestedSampler via default_dynamic=False.
    default_dynamic: bool = True

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.dynesty")
        self.vars: list = []
        self._dim = 0
        self._nlive = 100
        self._dlogz = 0.5
        self._rseed = 0
        self._run_nested_kwargs: dict[str, Any] = {}
        self._constructor_kwargs: dict[str, Any] = {}
        self._selectionexp: str | None = None
        self._sampler = None  # dynesty Nested/Dynamic sampler
        self._finished = False
        self._summary: dict[str, Any] | None = None
        self._use_dynamic = bool(self.default_dynamic)
        self._dynesty_csv_path: str | None = None

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._dim = len(self.vars)
        if self._dim < 1:
            raise ValueError(f"{self.method} requires at least one Sampling.Variable")
        bounds = sampling.get("Bounds") or {}
        if not isinstance(bounds, Mapping):
            bounds = {}

        # --- meta / HEP-owned ---
        self._nlive = max(2, int(bounds.get("nlive", bounds.get("n_live", 100)) or 100))
        if bounds.get("dlogz") is not None:
            self._dlogz = float(bounds.get("dlogz"))
        elif bounds.get("dlogz_init") is not None:
            self._dlogz = float(bounds.get("dlogz_init"))
        else:
            self._dlogz = 0.5
        self._rseed = int(
            bounds.get("rseed", sampling.get("Seed", sampling.get("seed", 0))) or 0
        )
        self._init_seed_sequence(self._rseed)
        self._selectionexp = sampling.get("selection")
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))

        # Engine fixed by Method (Dynesty=dynamic, MultiNest=static). Ignore any
        # leftover Bounds.dynamic / Dynamic keys (validation rejects them).
        self._use_dynamic = bool(self.default_dynamic)
        if "dynamic" in bounds or "Dynamic" in bounds:
            self._logger.warning(
                "%s ignores Bounds.dynamic/Dynamic (engine is fixed by Method: "
                "Dynesty=DynamicNestedSampler, MultiNest=NestedSampler static); "
                "remove the key from the task YAML",
                self.method,
            )

        # --- official constructor kwargs ---
        self._constructor_kwargs = extract_nested_constructor_kwargs(
            bounds, logger=self._logger
        )
        # nlive: Bounds.nlive is canonical; sampler.nlive may override if set.
        if "nlive" not in self._constructor_kwargs:
            self._constructor_kwargs["nlive"] = self._nlive
        else:
            self._nlive = max(2, int(self._constructor_kwargs["nlive"] or self._nlive))
            self._constructor_kwargs["nlive"] = self._nlive

        # --- official run_nested kwargs (static/dynamic-aware dlogz mapping) ---
        self._run_nested_kwargs = extract_run_nested_kwargs(
            bounds,
            dynamic=self._use_dynamic,
            dlogz_default=self._dlogz,
            logger=self._logger,
        )
        # Default progress to Jarvis logger path; user may set print_progress: false.
        self._run_nested_kwargs.setdefault("print_progress", True)

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
        """Disabled control-side toy logL (must never run in a real scan).

        Historical bug: Dynesty's internal ``sample()`` called ``loglikelihood(v)``
        on the control process, and this method returned a unit-cube Gaussian
        instead of Worker physics — inflating ``ncall`` without SAMPLE archives.
        Production wiring always uses :meth:`_make_redis_loglike`.
        """
        raise RuntimeError(
            f"{self.method}: control-side toy logL is disabled. "
            "All likelihood evaluations must go through RedisEvaluationPool "
            "(pool.evaluate_logl / Worker). If you see this, the Dynesty "
            "loglikelihood bridge was not installed."
        )

    def _make_redis_loglike(self, pool: RedisEvaluationPool) -> Any:
        """loglikelihood callable that always evaluates through Redis Workers.

        Dynesty calls this both via ``pool.map(LogLikelihood, batch)`` and
        *directly* inside ``pool.map(internal_sampler.sample, …)`` (which the
        Redis pool runs on the control process). Direct calls must hit
        Workers — never a control-side toy Gaussian.
        """

        def loglikelihood(params: Any) -> float:
            return float(pool.evaluate_logl(params))

        loglikelihood.__name__ = "loglikelihood"
        return loglikelihood

    def _build_constructor_kwargs(
        self,
        *,
        pool: RedisEvaluationPool,
        rstate: np.random.Generator,
    ) -> dict[str, Any]:
        """Merge user constructor kwargs with HEP-injected runtime handles."""
        kwargs = dict(self._constructor_kwargs)
        # Always Redis — never _local_loglike (toy).
        kwargs["loglikelihood"] = self._make_redis_loglike(pool)
        kwargs["prior_transform"] = _jarvis_prior_transform
        kwargs["ndim"] = self._dim
        kwargs["nlive"] = self._nlive
        kwargs["pool"] = pool
        kwargs["rstate"] = rstate
        # queue_size: user override wins; else batch_size / pool size.
        kwargs.setdefault("queue_size", max(1, self._batch_size))
        return kwargs

    def run_adaptive(self, *, timeout: float = 3600.0) -> int:
        self._require_redis(f"{self.method}.run_adaptive")
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

        # Progress/warnings use this sampler's Jarvis logger. MultiNest reuses
        # the vendored engine but must stamp ·•· MultiNest.Inner, not Dynesty.Inner.
        inner_logger = bind_inner(self._logger, method=self.method)
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
        ctor = self._build_constructor_kwargs(pool=pool, rstate=rstate)
        # Final filter against the live constructor signature (NestedSampler uses __new__).
        factory = getattr(sampler_cls, "__new__", sampler_cls)
        if factory is object.__new__:
            factory = sampler_cls
        # NestedSampler.__new__ has the full parameter list; DynamicNestedSampler.__init__ too.
        if sampler_cls is NestedSampler:
            ctor = filter_kwargs_to_callable(
                ctor, NestedSampler.__new__, context="constructor", logger=self._logger
            )
        else:
            ctor = filter_kwargs_to_callable(
                ctor,
                DynamicNestedSampler.__init__,
                context="constructor",
                logger=self._logger,
            )

        self._logger.info(
            "Starting %s nlive=%d ndim=%d dynamic=%s bound=%s sample=%s",
            sampler_cls.__name__,
            self._nlive,
            self._dim,
            self._use_dynamic,
            ctor.get("bound", "multi"),
            ctor.get("sample", "auto"),
        )
        self._sampler = sampler_cls(**ctor)

        # Attach logger onto nested sampler instance (used by run_nested progress).
        try:
            self._sampler.logger = inner_logger
            for attr in ("sampler", "base_sampler"):
                inner = getattr(self._sampler, attr, None)
                if inner is not None and hasattr(inner, "logger"):
                    inner.logger = inner_logger
                elif inner is not None:
                    setattr(inner, "logger", inner_logger)
        except Exception:
            pass

        run_kwargs = filter_kwargs_to_callable(
            self._run_nested_kwargs,
            self._sampler.run_nested,
            context="run_nested",
            logger=self._logger,
        )
        self._logger.info(
            "run_nested kwargs → %s",
            {k: v for k, v in run_kwargs.items() if k != "print_func"},
        )
        self._sampler.run_nested(**run_kwargs)
        self._finished = True
        self._summary = self._build_summary()
        # V1 finalize path: clean nested CSV for Jarvis-PLOT dynesty_runplot.
        try:
            csv_path = self.save_dynesty_results_to_csv()
            if self._summary is not None and csv_path:
                path_text = str(csv_path)
                if (
                    self.method.lower() == "multinest"
                    or "multinest" in path_text.lower()
                ):
                    self._summary["multinest_result_csv"] = path_text
                    self._summary.pop("dynesty_result_csv", None)
                else:
                    self._summary["dynesty_result_csv"] = path_text
        except Exception as exc:
            self._logger.warning("failed to write nested result CSV: %s", exc)
        # D13.6: sampler_summary.json (logZ, ncall, …) under DATABASE/.
        try:
            from jarvishep2.Sampling.diagnostics_export import export_nested_diagnostics

            written = export_nested_diagnostics(
                self._task_result_dir(), summary=self._summary
            )
            if written and self._summary is not None:
                self._summary.setdefault("diagnostics", {}).update(written)
            for key, path in written.items():
                self._logger.info("%s diagnostics %s → %s", self.method, key, path)
        except Exception as exc:
            self._logger.warning("failed to write nested sampler_summary: %s", exc)
        self.checkpoint_at_barrier(reason="dynesty_finished")
        # Official nested summary block (nlive/niter/ncall/eff/logz) → sampler.log
        try:
            self._log_nested_run_summary()
        except Exception as exc:
            self._logger.warning("failed to log nested run summary: %s", exc)
        self._logger.info(
            "%s finished logZ=%.4f ± %.4f niter=%s ncall=%s",
            self.method,
            self._summary.get("logz", float("nan")),
            self._summary.get("logzerr", float("nan")),
            self._summary.get("niter"),
            self._summary.get("ncall"),
        )
        return int(self._summary.get("ncall") or 0)

    def _log_nested_run_summary(self) -> None:
        """Emit dynesty/MultiNest ``Results.summary()`` block to the sampler logger."""
        if self._sampler is None:
            return
        try:
            results = self._sampler.results
        except Exception:
            return
        summary_fn = getattr(results, "summary", None)
        text = None
        if callable(summary_fn):
            try:
                text = summary_fn()
            except Exception:
                text = None
        if not text:
            # Fallback from our _summary dict when Results.summary is unavailable.
            s = self._summary or {}
            nlive = s.get("nlive", self._nlive)
            niter = s.get("niter", "")
            ncall = s.get("ncall", "")
            logz = s.get("logz", float("nan"))
            logzerr = s.get("logzerr", float("nan"))
            try:
                eff = 100.0 * float(niter) / float(ncall) if ncall else float("nan")
            except Exception:
                eff = float("nan")
            text = (
                "Summary\n"
                "=======\n"
                f"nlive: {int(nlive)}\n"
                f"niter: {int(niter) if niter != '' else niter}\n"
                f"ncall: {int(ncall) if ncall != '' else ncall}\n"
                f"eff(%): {eff:6.3f}\n"
                f"logz: {float(logz):6.3f} +/- {float(logzerr):6.3f}"
            )
        # One multi-line INFO record → sampler.log (module Sampler:dynesty/multinest).
        self._logger.info("%s run summary\n%s", self.method, str(text).rstrip())

    def _task_result_dir(self) -> str:
        return str(
            self.config.get("task_result_dir")
            or (self.config.get("Runtime") or {}).get("task_result_dir")
            or os.getcwd()
        )

    def dynesty_result_csv_path(self, task_result_dir: str | None = None) -> str:
        """Canonical path: ``<task_result_dir>/DATABASE/dynesty_result.csv``."""
        root = os.path.abspath(str(task_result_dir or self._task_result_dir()))
        return os.path.join(root, "DATABASE", "dynesty_result.csv")

    def save_dynesty_results_to_csv(
        self,
        output_csv: str | None = None,
    ) -> str | None:
        """Export clean nested dead-point CSV (not the full logL-call archive).

        Writes ``DATABASE/dynesty_result.csv`` for Jarvis-PLOT and post-run
        analysis. Row count is ``niter`` (nested samples), not Worker ncall.

        Returns the written path, or ``None`` if the sampler has no results yet.
        """
        if self._sampler is None:
            self._logger.warning("save_dynesty_results_to_csv: no sampler results")
            return None
        try:
            results = self._sampler.results
        except Exception as exc:
            self._logger.warning("save_dynesty_results_to_csv: cannot read results: %s", exc)
            return None
        path = output_csv or self.dynesty_result_csv_path()
        param_names = [
            str(getattr(var, "name", "") or "").strip()
            for var in (self.vars or [])
        ]
        written = export_dynesty_results_csv(
            results,
            path,
            fallback_nlive=self._nlive,
            param_names=param_names or None,
        )
        self._dynesty_csv_path = written
        self._logger.info(
            "%s clean nested CSV saved → %s (niter dead points; not all logL calls)",
            self.method,
            written,
        )
        return written

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
                "use_dynamic": self._use_dynamic,
                "constructor_kwargs": dict(self._constructor_kwargs),
                "run_nested_kwargs": dict(self._run_nested_kwargs),
                "dynesty_csv_path": self._dynesty_csv_path,
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
        if "use_dynamic" in state:
            self._use_dynamic = bool(state.get("use_dynamic"))
        ctor = state.get("constructor_kwargs")
        if isinstance(ctor, Mapping):
            self._constructor_kwargs = dict(ctor)
        run = state.get("run_nested_kwargs")
        if isinstance(run, Mapping):
            self._run_nested_kwargs = dict(run)
        path = state.get("dynesty_csv_path")
        self._dynesty_csv_path = str(path) if path else None


def create_dynesty() -> DynestySampler:
    return DynestySampler()


__all__ = [
    "DYNAMIC_RUN_NESTED_KEYS",
    "DynestySampler",
    "NESTED_CONSTRUCTOR_USER_KEYS",
    "STATIC_RUN_NESTED_KEYS",
    "create_dynesty",
    "export_dynesty_results_csv",
    "extract_nested_constructor_kwargs",
    "extract_run_nested_kwargs",
    "filter_kwargs_to_callable",
    "_jarvis_prior_transform",
]
