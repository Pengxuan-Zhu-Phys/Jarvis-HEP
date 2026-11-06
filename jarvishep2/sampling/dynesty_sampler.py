#!/usr/bin/env python3
"""Dynesty nested sampling on the V2 Redis runtime (D13.5, route A).

Uses the **vendored** dynesty package under
``jarvishep2.sampling.Source.Dynesty`` (stock 3.1.0 + Jarvis UUID channel),
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
import pickle
from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any
from uuid import uuid4

import numpy as np

from jarvishep2.sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.sampling.redis_evaluation_pool import RedisEvaluationPool
from jarvishep2.redis_queue import RedisQueue
from jarvishep2.sampling.sampling_utils import evaluate_selection, physical_from_u
from jarvishep2.sampling.variables import load_variables
from jarvishep2.contracts.nested import NESTED_CONSTRUCTOR_USER_KEYS
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_checkpoint_config, get_runtime_block
from jarvishep2.sample import Sample
from jarvishep2.sampling.sample_checkpoint import FORMAT, VERSION, config_identity

# Side-file basename next to ``runtime.pickle``: a copy of the engine snapshot
# (live points, saved_run, rstate, bounds, internal_state, …) that
# ``dynesty.utils.restore_sampler`` can read.  ``--resume`` never reads it.
NESTED_ENGINE_FILENAME = "nested_engine.pkl"
# ``native_sampler_blob`` holds the bytes dynesty's ``save_sampler`` writes.
NATIVE_SAMPLER_FORMAT = "dynesty_save_v1"
_DYNESTY_SAVE_FORMAT_VERSION = 1

# NestedSampler.run_nested (static) — user-facing only.
# Checkpoint cadence/path/resume are HEP-owned (EnvReqs.V2 + runtime).
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
    }
)

# DynamicNestedSampler.run_nested — user-facing only (same HEP ownership rule).
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
    }
)

# Injected by HEP at run time — not authored under Sampling.Bounds.run_nested.
HEP_OWNED_RUN_NESTED_KEYS: frozenset[str] = frozenset(
    {
        "checkpoint_file",
        "checkpoint_every",
        "resume",
    }
)

# Bounds keys owned by HEP / meta (never forwarded as constructor kwargs).
_BOUNDS_META_KEYS: frozenset[str] = frozenset(
    {
        "nlive",
        "seed",
        "dlogz",
        "dlogz_init",
        "run_nested",
        "sampler",
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


class NestedRedisLogL:
    """Pickle-friendly loglikelihood bridge for dynesty / MultiNest.

    Dynesty's native ``save`` / ``restore`` (and our runtime checkpoint) pickle
    the full engine — including this callable. Nested closures over
    :class:`RedisEvaluationPool` are **not** pickleable (they close over Redis
    sockets). This class strips the pool on pickle and re-attaches it on resume
    (V1 ``NestedLikelihoodBridge`` pattern).

    Do **not** use dill here: the design contract is stdlib pickle for the
    engine blob plus explicit runtime reattachment (see DESIGN_RESUME §16.1 and
    DESIGN_SAMPLERS §3.4).
    """

    __slots__ = ("_pool",)
    # Dynesty / RedisEvaluationPool detect logL callables by name / label.
    __name__ = "loglikelihood"

    def __init__(self, pool: RedisEvaluationPool | None = None) -> None:
        self._pool = pool

    def attach_pool(self, pool: RedisEvaluationPool | None) -> None:
        self._pool = pool

    def __call__(self, params: Any) -> float:
        if self._pool is None:
            raise RuntimeError(
                "NestedRedisLogL has no attached RedisEvaluationPool; "
                "call attach_pool() after restore before run_nested(resume=True)"
            )
        return float(self._pool.evaluate_logl(params))

    def __getstate__(self) -> dict[str, Any]:
        # Never pickle Redis / pool handles.
        return {}

    def __setstate__(self, state: Mapping[str, Any] | None) -> None:
        self._pool = None

    def __repr__(self) -> str:
        attached = self._pool is not None
        return f"NestedRedisLogL(pool_attached={attached})"


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
    2. Nested ``Bounds.sampler`` block (preferred)

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
    for block_key in ("sampler",):
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
    # Checkpoint cadence lives only under EnvReqs.V2.checkpoint.heartbeat.
    for banned in HEP_OWNED_RUN_NESTED_KEYS:
        if banned in run:
            run.pop(banned, None)
            if logger is not None:
                logger.warning(
                    "nested sampler: ignoring Sampling.Bounds.run_nested.%s "
                    "(HEP-owned; set EnvReqs.V2.checkpoint.heartbeat for "
                    "checkpoint interval in seconds, resume path is automatic)",
                    banned,
                )

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
        from jarvishep2.sampling.Source.Dynesty.py.dynesty.utils import (
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


class DynestySampler(CheckpointedSampler):
    """YAML ``Sampling.Method: Dynesty``.

    Always uses **DynamicNestedSampler**. Static nested sampling is
    ``Sampling.Method: MultiNest`` only — there is no ``Bounds.dynamic`` switch.
    Constructor + run_nested kwargs follow the official dynesty 3.x surface.
    """

    method = "Dynesty"
    # MultiNest subclasses force static NestedSampler via default_dynamic=False.
    default_dynamic: bool = True
    # The vendored engine snapshot in runtime.pickle is the resume source.
    # Workers still publish logL; HDF5 is not part of this checkpoint.
    uses_sample_checkpoint = True

    def __init__(self) -> None:
        super().__init__()
        self._seed = 0
        self._seed_seq: np.random.SeedSequence | None = None
        self._batch_size = 16
        self._pending_uuids: set[str] = set()
        self._generation = 0
        self._on_failure = "reject"
        self._logger = get_jarvis_logger("sampler.dynesty")
        self.vars: list = []
        self._mapper_pipeline = None
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
        # Native engine resume (dynesty save format inside runtime.pickle).
        self._native_sampler_blob: bytes | None = None
        self._engine_checkpoint_path: str | None = None
        self._pool_call_index = 0
        self._pool: RedisEvaluationPool | None = None
        self._loglike_bridge: NestedRedisLogL | None = None
        self._checkpoint_every_sec = 60.0
        self._parallel_workers = 1
        # Only True after import_runtime_state (Core --resume).
        self._resume_engine_requested = False
        self._wrapped_dynesty_save = False
        self._pickle_writer: Any | None = None
        self._owns_checkpoint = False
        self._config_id = ""

    @property
    def _accepted_index(self) -> int:
        """Next ``sample_index``: the number of likelihood tasks pushed so far."""
        if self._pool is not None:
            return int(self._pool._call_index)
        return int(self._pool_call_index or 0)

    def set_resume_mode(self, enabled: bool) -> None:
        super().set_resume_mode(enabled)
        if not enabled and self._resume_engine_requested:
            # A discarded checkpoint must not leave its engine or cursor behind.
            self._resume_engine_requested = False
            self._sampler = None
            self._native_sampler_blob = None
            self._pool_call_index = 0
            self._finished = False
            self._summary = None
            self._last_safe_state = None
            self._checkpoint_snapshots.clear()

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        from jarvishep2.mapper import MapperPipeline

        self._mapper_pipeline = MapperPipeline.from_config(self.config)
        self._dim = len(self.vars)
        if self._dim < 1:
            raise ValueError(f"{self.method} requires at least one Sampling.Variable")
        bounds = sampling.get("Bounds") or {}
        if not isinstance(bounds, Mapping):
            bounds = {}

        # --- meta / HEP-owned ---
        self._nlive = max(2, int(bounds.get("nlive", 100) or 100))
        if bounds.get("dlogz") is not None:
            self._dlogz = float(bounds.get("dlogz"))
        elif bounds.get("dlogz_init") is not None:
            self._dlogz = float(bounds.get("dlogz_init"))
        else:
            self._dlogz = 0.5
        # seed lives under Bounds only.
        self._rseed = int(bounds.get("seed", 0) or 0)
        self._init_seed_sequence(self._rseed)
        self._selectionexp = sampling.get("selection")
        workers = int(runtime.get("workers", 1) or 1)
        self._parallel_workers = max(1, workers)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        # Single source of truth: EnvReqs.V2.checkpoint.heartbeat (seconds).
        ckpt = get_checkpoint_config(self.config)
        self._checkpoint_every_sec = float(ckpt["heartbeat"])
        self._checkpoint_heartbeat_sec = self._checkpoint_every_sec

        # Engine is fixed by Method (Dynesty=dynamic, MultiNest=static).
        self._use_dynamic = bool(self.default_dynamic)

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
        # Strip any HEP-owned keys that slipped past extract (defense in depth).
        for banned in HEP_OWNED_RUN_NESTED_KEYS:
            self._run_nested_kwargs.pop(banned, None)
        self._engine_checkpoint_path = self._default_engine_checkpoint_path()
        self._config_id = config_identity(self.config)

    def _init_seed_sequence(self, seed: int) -> None:
        self._seed = int(seed)
        self._seed_seq = np.random.SeedSequence(self._seed)

    def _ensure_seed_sequence(self) -> np.random.SeedSequence:
        if self._seed_seq is None:
            self._seed_seq = np.random.SeedSequence(self._seed)
        return self._seed_seq

    def _require_redis(self, what: str) -> RedisQueue:
        if self.redis is None:
            raise RuntimeError(f"{what} requires redis")
        return self.redis

    def at_safe_barrier(self) -> bool:
        """Nested sampling checkpoints via dynesty save, not generation barriers."""
        return True

    def checkpoint_at_barrier(
        self,
        *,
        reason: str = "dynesty_finished",
        timeout: float = 3600.0,
    ) -> bool:
        """Finished run. One last engine snapshot; HDF5 is not consulted."""
        del reason, timeout
        return self.save_checkpoint(force=True)

    def _feedback_export_state(self) -> dict[str, Any]:
        """Keep historical checkpoint keys so --resume payloads stay compatible."""
        return {
            "generation": int(self._generation or 0),
            "pending_uuids": sorted(self._pending_uuids),
            "seed": self._seed,
            "submitted_uuids": list(self._submitted_uuids),
            "control_state": self._checkpoint_control_state(),
            "on_failure": self._on_failure,
        }

    def _feedback_import_state(self, state: Mapping[str, Any]) -> None:
        self._generation = int(state.get("generation", 0) or 0)
        self._pending_uuids = {str(item) for item in (state.get("pending_uuids") or [])}
        self._seed = int(state.get("seed", self._seed) or self._seed)
        self._seed_seq = np.random.SeedSequence(self._seed)
        self._on_failure = str(state.get("on_failure") or self._on_failure or "reject")
        self._import_checkpoint_control_state(state)
        pending = self._pending_uuids
        if not pending:
            self._submitted_uuids = []
        else:
            self._submitted_uuids = [
                str(uuid) for uuid in self._submitted_uuids if str(uuid) in pending
            ]

    def _build_sample_for_pool(self, payload: np.ndarray, uuid: str) -> Sample:
        """Build a Sample from unit-cube coords (uuid already assigned)."""
        u = np.asarray(payload, dtype=np.float64).reshape(-1)
        if u.size != self._dim:
            raise ValueError(f"{self.method}: expected {self._dim} coordinates, got {u.size}")
        sample = self._build_sample(u)
        sample.uuid = uuid
        return sample

    def _accepts_pool_payload(self, u: np.ndarray) -> bool:
        if u.size != self._dim:
            raise ValueError(f"{self.method}: expected {self._dim} coordinates, got {u.size}")
        if not self._selectionexp:
            return True
        return evaluate_selection(
            self._selectionexp,
            physical_from_u(u, self.vars, self._mapper_pipeline),
            context=self._expression_context,
        )

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

    def _make_redis_loglike(self, pool: RedisEvaluationPool) -> NestedRedisLogL:
        """Return a pickle-friendly logL bridge bound to *pool*.

        Dynesty calls this both via ``pool.map(LogLikelihood, batch)`` and
        *directly* inside ``pool.map(internal_sampler.sample, …)`` (which the
        Redis pool runs on the control process). Direct calls must hit
        Workers — never a control-side toy Gaussian.
        """
        bridge = NestedRedisLogL(pool)
        self._loglike_bridge = bridge
        return bridge

    def _build_constructor_kwargs(
        self,
        *,
        pool: RedisEvaluationPool,
        rstate: np.random.Generator,
    ) -> dict[str, Any]:
        """Merge user constructor kwargs with HEP-injected runtime handles."""
        kwargs = dict(self._constructor_kwargs)
        # Always Redis — never _local_loglike (toy). Pickle-safe bridge.
        kwargs["loglikelihood"] = self._make_redis_loglike(pool)
        kwargs["prior_transform"] = _jarvis_prior_transform
        kwargs["ndim"] = self._dim
        kwargs["nlive"] = self._nlive
        kwargs["pool"] = pool
        kwargs["rstate"] = rstate
        # Transport chunk size must not limit the number of concurrent walks.
        kwargs.setdefault("queue_size", self._parallel_workers)
        return kwargs

    # ------------------------------------------------------------------ resume
    def _default_engine_checkpoint_path(self) -> str:
        """``checkpoints/<scan>/<Method>/nested_engine.pkl`` next to runtime.pickle."""
        task_root = str(
            self.config.get("task_root")
            or (self.config.get("Runtime") or {}).get("task_root")
            or os.getcwd()
        )
        scan_block = self.config.get("Scan")
        scan_mapping = scan_block if isinstance(scan_block, Mapping) else {}
        scan_name = str(
            self.config.get("scan_name") or scan_mapping.get("name") or "scan"
        )
        method_dir = str(self.method or "Dynesty")
        return os.path.join(
            os.path.abspath(task_root),
            "checkpoints",
            scan_name,
            method_dir,
            NESTED_ENGINE_FILENAME,
        )

    def _engine_path(self) -> str:
        if self._engine_checkpoint_path:
            return str(self._engine_checkpoint_path)
        path = self._default_engine_checkpoint_path()
        self._engine_checkpoint_path = path
        return path

    def _iter_nested_sampler_layers(self, sampler: Any) -> list[Any]:
        """Outer DynamicNestedSampler + inner static / batch samplers."""
        layers = [sampler]
        for attr in ("sampler", "batch_sampler", "base_sampler"):
            inner = getattr(sampler, attr, None)
            if inner is not None and inner not in layers:
                layers.append(inner)
        return layers

    def _reattach_native_runtime(
        self,
        sampler: Any,
        *,
        pool: RedisEvaluationPool,
        bridge: NestedRedisLogL,
        logger: Any | None = None,
    ) -> None:
        """Re-bind pool / logL / prior_transform / logger after pickle restore.

        Dynesty's ``__getstate__`` strips ``pool`` and ``mapper``; the logL
        bridge comes back with ``_pool is None``. Without this step,
        ``run_nested(resume=True)`` would crash or silently use a dead handle.
        """
        bridge.attach_pool(pool)
        self._loglike_bridge = bridge
        queue_size = int(self._constructor_kwargs.get("queue_size") or self._parallel_workers)
        use_pool = dict(self._constructor_kwargs.get("use_pool") or {})
        for layer in self._iter_nested_sampler_layers(sampler):
            layer.use_pool = dict(use_pool)
            for attr, key in (
                ("use_pool_ptform", "prior_transform"),
                ("use_pool_logl", "loglikelihood"),
                ("use_pool_evolve", "propose_point"),
                ("use_pool_update", "update_bound"),
                ("use_pool_stopfn", "stop_function"),
            ):
                setattr(layer, attr, bool(use_pool.get(key, True)))
            layer.queue_size = queue_size if layer.use_pool_evolve else 1
            try:
                layer.pool = pool
                layer.mapper = pool.map
            except Exception:
                pass
            try:
                layer.prior_transform = _jarvis_prior_transform
            except Exception:
                pass
            # loglikelihood may be LogLikelihood wrapper or the bridge itself.
            logl_obj = getattr(layer, "loglikelihood", None)
            if logl_obj is None:
                continue
            if isinstance(logl_obj, NestedRedisLogL):
                logl_obj.attach_pool(pool)
            elif hasattr(logl_obj, "loglikelihood"):
                inner = getattr(logl_obj, "loglikelihood", None)
                if isinstance(inner, NestedRedisLogL):
                    inner.attach_pool(pool)
                else:
                    # Replaced with a fresh pickle-safe bridge.
                    try:
                        logl_obj.loglikelihood = bridge
                    except Exception:
                        pass
            if logger is not None:
                try:
                    layer.logger = logger
                except Exception:
                    pass
        if logger is not None:
            try:
                sampler.logger = logger
            except Exception:
                pass

    def _strip_dynesty_save_hook(self, engine: Any | None = None) -> None:
        """Remove non-pickleable instance ``save`` override before engine pickle."""
        target = engine if engine is not None else self._sampler
        if target is None:
            return
        d = getattr(target, "__dict__", None)
        if isinstance(d, dict):
            d.pop("save", None)
            d.pop("_jarvis_raw_save", None)
        self._wrapped_dynesty_save = False

    def _pickle_native_sampler(self, sampler: Any | None = None) -> bytes | None:
        """Serialize the engine in dynesty's ``save_sampler`` format (stdlib pickle)."""
        from jarvishep2.sampling.Source.Dynesty.py.dynesty.utils import DYNESTY_VERSION

        target = sampler if sampler is not None else self._sampler
        if target is None:
            return None
        reinstall = target is self._sampler
        self._strip_dynesty_save_hook(target)
        try:
            return pickle.dumps(
                {
                    "sampler": target,
                    "version": DYNESTY_VERSION,
                    "format_version": _DYNESTY_SAVE_FORMAT_VERSION,
                },
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        except Exception as exc:
            self._logger.warning(
                "%s: failed to pickle native dynesty engine for checkpoint: %s",
                self.method,
                exc,
            )
            return None
        finally:
            if reinstall:
                self._install_dynesty_save_hook()

    def _install_dynesty_save_hook(self) -> None:
        """Turn dynesty's timed ``save`` into a Jarvis engine snapshot.

        Dynesty calls ``save(checkpoint_file)`` on its **sampling thread**
        between iterations, every ``checkpoint_every`` seconds and at the end
        of the run.  No likelihood task is in flight then, so the pool cursor
        is exact.  The hook is stripped before every pickle so the local
        function is never written into the engine.
        """
        engine = self._sampler
        if engine is None or self._wrapped_dynesty_save:
            return
        if not callable(getattr(type(engine), "save", None)):
            return
        sampler_self = self

        def _save_with_jarvis_state(fname: str, *args: Any, **kwargs: Any) -> None:
            # Dynesty's timer is the heartbeat. The hook returns immediately
            # when the checkpoint process is busy, and never waits for HDF5.
            _ = fname, args, kwargs
            try:
                sampler_self.save_checkpoint()
            except Exception as exc:
                sampler_self._logger.warning(
                    "%s: engine checkpoint was not saved -> %s", sampler_self.method, exc
                )

        engine.save = _save_with_jarvis_state  # type: ignore[method-assign]
        self._wrapped_dynesty_save = True

    def checkpoint_state(self) -> dict[str, Any] | None:
        """Engine bytes at the current pool cursor.

        Called only when the checkpoint process can accept a snapshot, or when
        the run forces one. Pickle stays on this thread so the bytes match the
        engine between iterations; the child does the fsync.
        """
        cursor = self._accepted_index
        blob = self._native_sampler_blob
        if self._sampler is not None:
            blob = self._pickle_native_sampler()
            if blob is None:
                return None
        self._pool_call_index = cursor
        state = self._engine_state(blob=blob, cursor=cursor)
        state["next_sample_index"] = int(cursor)
        return state

    def checkpoint_runtime_state(self, *, safe: bool) -> dict[str, Any]:
        """Live engine snapshot. ``safe`` does not consult DATABASE."""
        del safe
        state = self.checkpoint_state()
        if state is None:
            return self._initial_engine_state()
        return state

    def _initial_engine_state(self) -> dict[str, Any]:
        """Recovery point before the first likelihood task: no engine, cursor 0."""
        state = self._engine_state(blob=None, cursor=0)
        state["finished"] = False
        state["summary"] = None
        state["next_sample_index"] = 0
        return state

    def _checkpoint_envelope(self, state: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "format": FORMAT,
            "version": VERSION,
            "method": self.method.lower(),
            "config_id": self._config_id,
            "integrity": {
                "unit_execution": (self.config.get("_unit_execution") or {}).get(
                    "execution_identity"
                ),
            },
            "sampler_state": dict(state),
            # The child writes this sibling from the blob already in the payload.
            "sidecars": [
                {
                    "name": NESTED_ENGINE_FILENAME,
                    "field": "sampler_state.native_sampler_blob",
                }
            ],
        }

    def prepare_fresh_points(self, checkpoint_file: str) -> None:
        """Remove an earlier run's engine copy. Resume never reads that file."""
        path = os.path.join(
            os.path.dirname(os.path.abspath(checkpoint_file)), NESTED_ENGINE_FILENAME
        )
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass

    def restore_checkpoint(self, state: Mapping[str, Any]) -> None:
        if not isinstance(state, Mapping):
            raise ValueError(f"{self.method} checkpoint has no state")
        self.import_runtime_state(state)

    def configure_checkpoint(self, *, checkpoint_file: str, writer: Any | None = None) -> None:
        """Attach the checkpoint process. Dynesty's own timer picks the safe point."""
        from jarvishep2.runtime.checkpoint_process import CheckpointProcess

        self._checkpoint_file = str(checkpoint_file or "")
        self._engine_checkpoint_path = os.path.join(
            os.path.dirname(os.path.abspath(self._checkpoint_file or "runtime.pickle")),
            NESTED_ENGINE_FILENAME,
        )
        policy = get_checkpoint_config(self.config)
        if not policy["enabled"] or self.config.get("_check_modules_sample_layout"):
            if not policy["enabled"]:
                self._logger.info(
                    "%s checkpoint disabled (EnvReqs.V2.checkpoint.enabled=false)",
                    self.method,
                )
            return
        if not checkpoint_file:
            raise ValueError("sampler checkpoint path is required")
        if not self._resume_engine_requested:
            self.prepare_fresh_points(checkpoint_file)
        if not self._config_id:
            self._config_id = config_identity(self.config)
        self._owns_checkpoint = writer is None
        if writer is None:
            try:
                writer = CheckpointProcess(
                    checkpoint_file,
                    self._logger,
                    interval=policy["heartbeat"],
                    scan_name=str(self.config.get("scan_name") or ""),
                )
            except Exception as exc:
                self._logger.warning("Checkpoint process unavailable -> %s", exc)
                return
            writer.start()
        self._pickle_writer = writer
        self._logger.info(
            "%s checkpoint: engine snapshot every %.0fs (dynesty); "
            "the checkpoint process writes runtime.pickle",
            self.method,
            float(self._checkpoint_every_sec),
        )
        self.save_checkpoint(force=True)

    def save_checkpoint(self, *, force: bool = False) -> bool:
        """Pickle the engine only when the checkpoint process can take it.

        A busy process returns before the pickle, so the sampling thread does
        not stall on disk. ``force`` waits, then submits the finished engine.
        """
        writer = self._pickle_writer
        if writer is None:
            return False
        if force and not writer.wait(timeout=600.0):
            self._logger.warning(
                "%s: checkpoint process was still busy; the snapshot was not submitted",
                self.method,
            )
            return False
        if not writer.available():
            return False
        state = self.checkpoint_state()
        if state is None:
            return False
        return writer.submit(self._checkpoint_envelope(state))

    def shutdown_checkpointing(self) -> None:
        writer = self._pickle_writer
        owns = self._owns_checkpoint
        self._pickle_writer = None
        self._owns_checkpoint = False
        super().shutdown_checkpointing()
        if writer is not None and owns:
            writer.finish()

    def _load_native_from_blob(self, blob: bytes) -> Any:
        """Return the engine in a ``dynesty_save_v1`` blob; ValueError if unreadable."""
        try:
            saved = pickle.loads(blob)
        except Exception as exc:
            raise ValueError(f"{self.method} engine snapshot cannot be read: {exc}") from exc
        if (
            not isinstance(saved, Mapping)
            or saved.get("format_version") != _DYNESTY_SAVE_FORMAT_VERSION
            or saved.get("sampler") is None
        ):
            raise ValueError(f"{self.method} engine snapshot has an unknown layout")
        from jarvishep2.sampling.Source.Dynesty.py.dynesty.utils import DYNESTY_VERSION

        if saved.get("version") != DYNESTY_VERSION:
            self._logger.warning(
                "%s: engine snapshot comes from dynesty %s; this release has %s",
                self.method,
                saved.get("version"),
                DYNESTY_VERSION,
            )
        return saved["sampler"]

    @staticmethod
    def _is_dynamic_engine(sampler: Any) -> bool:
        """True for DynamicNestedSampler (has batch/base state machine)."""
        if sampler is None:
            return False
        # DynamicNestedSampler exposes internal_state + nested .sampler.
        if type(sampler).__name__ == "DynamicNestedSampler":
            return True
        return hasattr(sampler, "internal_state") and hasattr(sampler, "sampler")

    def _engine_matches_method(self, sampler: Any) -> bool:
        """Dynesty requires Dynamic engine; MultiNest requires static NestedSampler."""
        is_dynamic = self._is_dynamic_engine(sampler)
        want_dynamic = bool(self.default_dynamic)
        return is_dynamic is want_dynamic

    @staticmethod
    def _engine_looks_finished(sampler: Any) -> bool:
        """True when dynesty would no-op ``run_nested(resume=True)``.

        * Dynamic: ``internal_state == RUN_DONE``
        * Static: ``added_live`` after a completed run with ``add_live=True``
        """
        if sampler is None:
            return False
        state = getattr(sampler, "internal_state", None)
        if state is not None:
            name = getattr(state, "name", None) or str(state)
            if str(name).endswith("RUN_DONE") or str(name) == "RUN_DONE":
                return True
        # Static NestedSampler (class name often ``Sampler``): resume is a no-op
        # once live points were merged via add_live.
        if bool(getattr(sampler, "added_live", False)) and not DynestySampler._is_dynamic_engine(
            sampler
        ):
            return True
        return False

    def _resolve_resume_engine(
        self,
        *,
        pool: RedisEvaluationPool,
        logger: Any | None = None,
    ) -> tuple[Any | None, bool]:
        """Reattach the engine imported from runtime.pickle for ``run_nested(resume=True)``.

        Only runs after ``import_runtime_state`` (Core ``--resume``).  The
        engine is the snapshot in that pickle; ``nested_engine.pkl`` is never
        a resume source.
        """
        if not self._resume_engine_requested or self._finished:
            return None, False
        sampler = self._sampler
        if sampler is None and self._native_sampler_blob:
            sampler = self._load_native_from_blob(self._native_sampler_blob)
        if sampler is None:
            return None, False
        bridge = NestedRedisLogL(pool)
        self._reattach_native_runtime(sampler, pool=pool, bridge=bridge, logger=logger)
        if self._engine_looks_finished(sampler):
            self._logger.info(
                "%s resume: the engine is already complete (RUN_DONE / added_live); "
                "not re-entering run_nested",
                self.method,
            )
            # run_adaptive skips run_nested and writes the result files again.
            self._finished = True
            return sampler, False
        self._logger.info(
            "%s resume: restored native dynesty engine (it=%s ncall=%s cursor=%d)",
            self.method,
            getattr(sampler, "it", "?"),
            getattr(sampler, "ncall", "?"),
            int(pool._call_index),
        )
        return sampler, True

    def run_adaptive(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Run nested sampling with a per-generation timeout setting.

        ``--resume`` continues the engine snapshot from runtime.pickle
        with ``run_nested(resume=True)``.
        """
        if timeout is not None:
            generation_timeout = timeout
        timeout = generation_timeout
        self._require_redis(f"{self.method}.run_adaptive")
        self._ensure_seed_sequence()
        if self._finished:
            self._resume_engine_requested = False
            return int((self._summary or {}).get("ncall") or 0)

        from jarvishep2.sampling.Source.Dynesty.py.dynesty import (
            DynamicNestedSampler,
            NestedSampler,
        )
        from jarvishep2.sampling.Source.Dynesty.py.dynesty.jarvis_logging import (
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
            accepts_payload=self._accepts_pool_payload,
            njobs=self._parallel_workers,
            on_submit=self.note_archive_submissions,
        )
        self._logger.info(
            "%s parallelism: workers=%d queue_size=%s transport_batch=%d propose_point=%s; "
            "each walk is sequential and each map waits for its slowest walk",
            self.method, pool.size,
            self._constructor_kwargs.get("queue_size", self._parallel_workers),
            self._batch_size,
            (self._constructor_kwargs.get("use_pool") or {}).get("propose_point", True),
        )
        # The pool cursor continues from the snapshot, so new tasks take the
        # sample_index values after the committed cursor.
        pool._call_index = max(0, int(self._pool_call_index or 0))
        self._pool = pool

        resume = False
        wanted_resume = bool(self._resume_engine_requested)
        restored, resume = self._resolve_resume_engine(pool=pool, logger=inner_logger)
        # Consume the one-shot resume request after the resolve attempt so a
        # later accidental re-entry cannot silently re-load a finished engine.
        self._resume_engine_requested = False
        if restored is not None:
            self._sampler = restored
            self._wrapped_dynesty_save = False  # restored object has stock .save
            self._logger.warning(
                "%s engine resumed -> iteration -> %s | ncall -> %s | state -> %s "
                "| run_nested resume -> %s",
                self.method,
                getattr(restored, "it", "?"),
                getattr(restored, "ncall", "?"),
                getattr(restored, "internal_state", "?"),
                resume,
            )
        else:
            if wanted_resume:
                self._logger.info(
                    "%s resume: the checkpoint holds no engine yet; the nested run "
                    "starts again from its seed and replaces the earlier rows",
                    self.method,
                )
            rstate = np.random.default_rng(self._rseed)
            sampler_cls = DynamicNestedSampler if self._use_dynamic else NestedSampler
            ctor = self._build_constructor_kwargs(pool=pool, rstate=rstate)
            # Final filter against the live constructor signature.
            if sampler_cls is NestedSampler:
                ctor = filter_kwargs_to_callable(
                    ctor,
                    NestedSampler.__new__,
                    context="constructor",
                    logger=self._logger,
                )
            else:
                ctor = filter_kwargs_to_callable(
                    ctor,
                    DynamicNestedSampler.__init__,
                    context="constructor",
                    logger=self._logger,
                )

            self._logger.info(
                "%s engine started fresh -> %s | nlive -> %d | ndim -> %d "
                "| bound -> %s | sample -> %s",
                self.method,
                sampler_cls.__name__,
                self._nlive,
                self._dim,
                ctor.get("bound", "multi"),
                ctor.get("sample", "auto"),
            )
            self._sampler = sampler_cls(**ctor)
            resume = False
            self._wrapped_dynesty_save = False

        # Attach logger onto nested sampler instance (used by run_nested progress).
        try:
            self._sampler.logger = inner_logger
            for attr in ("sampler", "base_sampler", "batch_sampler"):
                inner = getattr(self._sampler, attr, None)
                if inner is not None:
                    try:
                        inner.logger = inner_logger
                    except Exception:
                        setattr(inner, "logger", inner_logger)
        except Exception:
            pass

        # dynesty's timed save is the safe point between iterations.
        self._install_dynesty_save_hook()

        if not self._finished:
            self._run_nested_engine(pool=pool, resume=resume)

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
            self._logger.warning(
                "%s result CSV was not written -> %s", self.method, exc, exc_info=True
            )
        # D13.6: sampler_summary.json (logZ, ncall, …) under DATABASE/.
        try:
            from jarvishep2.sampling.diagnostics_export import export_nested_diagnostics

            written = export_nested_diagnostics(
                self._task_result_dir(), summary=self._summary
            )
            if written and self._summary is not None:
                self._summary.setdefault("diagnostics", {}).update(written)
            for key, path in written.items():
                self._logger.warning("%s diagnostics (%s) saved to %s", self.method, key, path)
        except Exception as exc:
            self._logger.warning(
                "%s sampler summary was not written -> %s", self.method, exc, exc_info=True
            )
        if get_checkpoint_config(self.config)["enabled"]:
            self.save_checkpoint(force=True)
        # Official nested summary block (nlive/niter/ncall/eff/logz) → sampler.log
        try:
            self._log_nested_run_summary()
        except Exception as exc:
            self._logger.warning(
                "%s run summary could not be logged -> %s", self.method, exc, exc_info=True
            )
        return int(self._summary.get("ncall") or 0)

    def _run_nested_engine(self, *, pool: RedisEvaluationPool, resume: bool) -> None:
        run_kwargs = filter_kwargs_to_callable(
            self._run_nested_kwargs,
            self._sampler.run_nested,
            context="run_nested",
            logger=self._logger,
        )
        # HEP injects the checkpoint interval from EnvReqs.V2 — never from
        # Sampling YAML.  The save hook ignores the file name.
        if get_checkpoint_config(self.config)["enabled"]:
            run_kwargs["checkpoint_file"] = self._engine_path()
            run_kwargs["checkpoint_every"] = float(self._checkpoint_every_sec)
        if resume:
            run_kwargs["resume"] = True
        else:
            # Fresh run: never pass resume=True (would no-op if state is weird).
            run_kwargs.pop("resume", None)

        self._logger.info(
            "run_nested kwargs → %s",
            {
                k: v
                for k, v in run_kwargs.items()
                if k not in {"print_func"}
            },
        )
        # An interrupt keeps the last landed runtime.pickle. The live engine
        # may be ahead of that file; resume continues from the file.
        try:
            self._sampler.run_nested(**run_kwargs)
            self._finished = True
        finally:
            self._pool_call_index = int(pool._call_index)

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

    # ------------------------------------------------------------------ logging
    def log_settings_rows(self) -> list[tuple[str, Any]]:
        rows: list[tuple[str, Any]] = [
            ("nlive", self._nlive),
            ("dlogz", self._dlogz),
            ("dynamic", "yes" if self._use_dynamic else "no"),
        ]
        for key in ("bound", "sample"):
            if key in self._constructor_kwargs:
                rows.append((key, self._constructor_kwargs[key]))
        if self._selectionexp:
            rows.append(("selection", self._selectionexp))
        rows.append(
            (
                "engine checkpoint",
                self._engine_path() if get_checkpoint_config(self.config)["enabled"] else "off",
            )
        )
        return rows

    def log_summary_rows(self) -> list[tuple[str, Any]]:
        summary = self._summary or {}
        if not summary.get("finished") or "niter" not in summary:
            return [("iterations", getattr(self._sampler, "it", "n/a"))]
        niter = int(summary.get("niter") or 0)
        ncall = int(summary.get("ncall") or 0)
        rows: list[tuple[str, Any]] = [
            ("nlive", summary.get("nlive", self._nlive)),
            ("iterations", niter),
            ("ncall", ncall),
            ("efficiency", self._efficiency_text(summary, niter, ncall)),
            (
                "logZ",
                f"{float(summary.get('logz', float('nan'))):.4f} ± "
                f"{float(summary.get('logzerr', float('nan'))):.4f}",
            ),
        ]
        for key in ("dynesty_result_csv", "multinest_result_csv"):
            if summary.get(key):
                rows.append(("result file", summary[key]))
        for key, path in dict(summary.get("diagnostics") or {}).items():
            rows.append((f"diagnostics {key}", path))
        return rows

    @staticmethod
    def _efficiency_text(summary: Mapping[str, Any], niter: int, ncall: int) -> str:
        # Same number as the engine's own summary block ("eff(%)").
        eff = summary.get("eff")
        if eff is not None and np.isfinite(float(eff)):
            return f"{float(eff):.2f} %"
        return f"{100.0 * niter / ncall:.2f} %" if ncall else "n/a"

    def log_stop_reason(self) -> str | None:
        if self._finished:
            # run_nested returns normally on convergence, budgets and custom
            # stopping rules. This vendored API does not report which fired.
            budgets = [key for key in ("maxiter", "maxcall", "maxbatch")
                       if self._run_nested_kwargs.get(key) is not None]
            detail = f"; configured budgets: {', '.join(budgets)}" if budgets else ""
            return f"engine returned (termination criterion not reported{detail})"
        return None

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
        already_reported = written is not None and written == self._dynesty_csv_path
        self._dynesty_csv_path = written
        # The runtime rewrites the same file after the run; report it once.
        (self._logger.debug if already_reported else self._logger.warning)(
            "%s results (one row per dead point) saved to %s",
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
            try:
                eff = float(res["eff"]) if "eff" in res.keys() else None
            except Exception:
                eff = None
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
                "eff": eff,
                "nlive": self._nlive,
                "samples_uid": samples_uid,
            }
        except Exception as exc:
            self._logger.warning(
                "%s summary could not be built -> %s", self.method, exc, exc_info=True
            )
            return {"method": self.method, "finished": self._finished, "error": str(exc)}

    def summary(self) -> dict[str, Any]:
        if self._summary is None:
            self._summary = self._build_summary()
        return dict(self._summary)

    def _engine_state(self, *, blob: bytes | None, cursor: int) -> dict[str, Any]:
        """Checkpoint fields for an engine snapshot taken at pool cursor *cursor*.

        ``native_sampler_blob`` holds dynesty's ``save_sampler`` bytes.  Redis,
        logger and pool are never pickled; they are reattached on resume.
        """
        state = self._feedback_export_state()
        # Strip non-pickleable keys from constructor kwargs (callables).
        safe_ctor = {
            k: v
            for k, v in dict(self._constructor_kwargs).items()
            if k
            not in {
                "loglikelihood",
                "prior_transform",
                "pool",
                "rstate",
            }
        }
        safe_run = {
            k: v
            for k, v in dict(self._run_nested_kwargs).items()
            if k not in {"print_func"}
        }
        state.update(
            {
                "method": self.method,
                "finished": self._finished,
                "summary": deepcopy(self._summary),
                "nlive": self._nlive,
                "dlogz": self._dlogz,
                "rseed": self._rseed,
                "use_dynamic": self._use_dynamic,
                "constructor_kwargs": safe_ctor,
                "run_nested_kwargs": safe_run,
                "dynesty_csv_path": self._dynesty_csv_path,
                "checkpoint_index": int(cursor),
                "pool_call_index": int(cursor),
                "engine_checkpoint_path": self._engine_path(),
                "native_sampler_blob": blob,
                "native_sampler_format": NATIVE_SAMPLER_FORMAT,
                "nested_engine_filename": NESTED_ENGINE_FILENAME,
            }
        )
        return state

    def export_runtime_state(self) -> dict[str, Any]:
        """Live engine state at the current pool cursor.

        This does not write a file. Resume reads the snapshot the checkpoint
        process has already landed.
        """
        blob = self._native_sampler_blob
        if self._sampler is not None:
            blob = self._pickle_native_sampler()
        return self._engine_state(blob=blob, cursor=self._accepted_index)

    def _validated_resume_engine(self, state: Mapping[str, Any]) -> Any | None:
        """Return the snapshot's engine, or raise ValueError if it cannot continue."""
        blob = state.get("native_sampler_blob")
        if blob is None:
            if self._state_index(state) > 0 and not state.get("finished"):
                raise ValueError(
                    f"{self.method} checkpoint at cursor {self._state_index(state)} has no engine"
                )
            return None
        if (
            state.get("native_sampler_format") != NATIVE_SAMPLER_FORMAT
            or not isinstance(blob, (bytes, bytearray))
        ):
            raise ValueError(
                f"{self.method} checkpoint engine format "
                f"{state.get('native_sampler_format')!r} comes from an earlier release"
            )
        if state.get("finished"):
            return None
        engine = self._load_native_from_blob(bytes(blob))
        if not self._engine_matches_method(engine):
            expected = "DynamicNestedSampler" if self.default_dynamic else "static NestedSampler"
            raise ValueError(
                f"{self.method} checkpoint holds a {type(engine).__name__} engine; "
                f"{self.method} needs {expected}"
            )
        return engine

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        """Restore the committed engine snapshot for ``run_adaptive``.

        Raises ValueError for a snapshot this release cannot continue (earlier
        format, unreadable engine, or the other method's engine).  Core then
        discards the checkpoint and starts fresh.
        """
        engine = self._validated_resume_engine(state)
        self._feedback_import_state(state)
        # Arm engine restore for the next run_adaptive (Core --resume only).
        self._resume_engine_requested = True
        self._finished = bool(state.get("finished", False))
        self._summary = state.get("summary")
        self._nlive = int(state.get("nlive", self._nlive) or self._nlive)
        self._dlogz = float(state.get("dlogz", self._dlogz) or self._dlogz)
        self._rseed = int(state.get("rseed", self._rseed) or self._rseed)
        if "use_dynamic" in state:
            self._use_dynamic = bool(state.get("use_dynamic"))
        runtime_ctor = {
            k: self._constructor_kwargs[k]
            for k in ("queue_size", "use_pool") if k in self._constructor_kwargs
        }
        ctor = state.get("constructor_kwargs")
        if isinstance(ctor, Mapping):
            # Drop any legacy non-serializable leftovers.
            self._constructor_kwargs = {
                str(k): v
                for k, v in ctor.items()
                if str(k)
                not in {
                    "loglikelihood",
                    "prior_transform",
                    "pool",
                    "rstate",
                }
            }
        for key in ("queue_size", "use_pool"):
            self._constructor_kwargs.pop(key, None)
        self._constructor_kwargs.update(runtime_ctor)
        run = state.get("run_nested_kwargs")
        if isinstance(run, Mapping):
            self._run_nested_kwargs = {
                str(k): v for k, v in run.items() if str(k) != "print_func"
            }
        path = state.get("dynesty_csv_path")
        self._dynesty_csv_path = str(path) if path else None
        self._pool_call_index = self._state_index(state)
        blob = state.get("native_sampler_blob")
        self._native_sampler_blob = bytes(blob) if blob is not None else None
        # Pool reattachment happens in run_adaptive.
        self._sampler = engine


def create_dynesty() -> DynestySampler:
    return DynestySampler()


__all__ = [
    "DYNAMIC_RUN_NESTED_KEYS",
    "DynestySampler",
    "HEP_OWNED_RUN_NESTED_KEYS",
    "NESTED_CONSTRUCTOR_USER_KEYS",
    "NESTED_ENGINE_FILENAME",
    "NestedRedisLogL",
    "STATIC_RUN_NESTED_KEYS",
    "create_dynesty",
    "export_dynesty_results_csv",
    "extract_nested_constructor_kwargs",
    "extract_run_nested_kwargs",
    "filter_kwargs_to_callable",
    "_jarvis_prior_transform",
]
