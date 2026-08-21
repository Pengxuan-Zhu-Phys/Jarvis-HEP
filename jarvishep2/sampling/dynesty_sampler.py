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
import tempfile
import time
from collections.abc import Mapping, Sequence
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
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample

# Side-file basename next to Jarvis ``state.pkl`` for dynesty's native engine
# pickle (live points, saved_run, rstate, bounds, internal_state, …).
NESTED_ENGINE_FILENAME = "nested_engine.pkl"

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
    # Checkpoint cadence lives only under EnvReqs.V2.checkpoint_heartbeat_sec.
    for banned in HEP_OWNED_RUN_NESTED_KEYS:
        if banned in run:
            run.pop(banned, None)
            if logger is not None:
                logger.warning(
                    "nested sampler: ignoring Sampling.Bounds.run_nested.%s "
                    "(HEP-owned; set EnvReqs.V2.checkpoint_heartbeat_sec for "
                    "checkpoint interval, resume path is automatic)",
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
        # Native engine resume (dynesty save/restore + runtime checkpoint).
        self._native_sampler_blob: bytes | None = None
        self._engine_checkpoint_path: str | None = None
        self._pool_call_index = 0
        self._loglike_bridge: NestedRedisLogL | None = None
        self._checkpoint_every_sec = 60.0
        # Only True after import_runtime_state / arm_engine_resume (Core --resume).
        # Prevents a leftover nested_engine.pkl from poisoning a fresh run.
        self._resume_engine_requested = False
        self._wrapped_dynesty_save = False

    def arm_engine_resume(self) -> None:
        """Core ``--resume``: allow loading ``nested_engine.pkl`` even without blob."""
        self._resume_engine_requested = True

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
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        # Single source of truth: EnvReqs.V2.checkpoint_heartbeat_sec (via runtime).
        # Used for dynesty nested_engine.pkl cadence.
        heartbeat = float(runtime.get("checkpoint_heartbeat_sec", 30.0) or 30.0)
        self._checkpoint_every_sec = max(1.0, heartbeat)
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

    def checkpoint_at_barrier(self, *, reason: str = "dynesty_finished") -> bool:
        """Persist nested engine + state.pkl after a completed (or interrupted) run."""
        redis = self.redis
        scan_block = self.config.get("Scan")
        scan_mapping = scan_block if isinstance(scan_block, Mapping) else {}
        scan_name = str(
            self.config.get("scan_name") or scan_mapping.get("name") or "scan"
        )
        if redis is not None and self._submitted_uuids:
            deadline = time.monotonic() + 5.0
            while time.monotonic() < deadline:
                persisted = redis.get_archived_uuids(scan_name)
                self.set_persisted_uuids(persisted)
                if set(self._submitted_uuids) <= persisted:
                    break
                time.sleep(0.05)
        return self.persist_runtime_checkpoint(force=True, reason=reason)

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

    def _build_sample_for_pool(self, payload: np.ndarray, uuid: str) -> Sample:
        """Build a Sample from unit-cube coords (uuid already assigned)."""
        u = np.asarray(payload, dtype=np.float64).reshape(-1)
        if u.size != self._dim:
            # If physical coords slipped through, clip to dim
            u = u[: self._dim]
        sample = self._build_sample(u)
        sample.uuid = uuid
        if self._selectionexp:
            physical = physical_from_u(u, self.vars, self._mapper_pipeline)
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
        # queue_size: user override wins; else batch_size / pool size.
        kwargs.setdefault("queue_size", max(1, self._batch_size))
        return kwargs

    # ------------------------------------------------------------------ resume
    def _default_engine_checkpoint_path(self) -> str:
        """``checkpoints/<scan>/<Method>/nested_engine.pkl`` next to state.pkl."""
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
        for layer in self._iter_nested_sampler_layers(sampler):
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
        """Serialize the dynesty engine with stdlib pickle (no dill)."""
        target = sampler if sampler is not None else self._sampler
        if target is None:
            return None
        reinstall = target is self._sampler
        self._strip_dynesty_save_hook(target)
        try:
            return pickle.dumps(target, protocol=pickle.HIGHEST_PROTOCOL)
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

    def _save_engine_side_file(self, sampler: Any | None = None) -> str | None:
        """Atomic dynesty-native save → ``nested_engine.pkl`` (sampling-thread safe)."""
        target = sampler if sampler is not None else self._sampler
        if target is None:
            return None
        path = self._engine_path()
        parent = os.path.dirname(os.path.abspath(path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        reinstall = target is self._sampler
        self._strip_dynesty_save_hook(target)
        try:
            # Stock class method — pickle-safe.
            type(target).save(target, path)
            return path
        except Exception as exc:
            self._logger.warning(
                "%s: dynesty engine save to %s failed: %s", self.method, path, exc
            )
            return None
        finally:
            if reinstall:
                self._install_dynesty_save_hook()

    def _install_dynesty_save_hook(self) -> None:
        """Wrap dynesty ``.save`` so each engine pickle also updates Jarvis state.pkl.

        Dynesty calls ``save(checkpoint_file)`` on its **sampling thread** every
        ``checkpoint_every`` seconds (and at run end).  The hook is stripped
        before every pickle so the local function is never written into the
        engine checkpoint.
        """
        engine = self._sampler
        if engine is None or self._wrapped_dynesty_save:
            return
        if not callable(getattr(type(engine), "save", None)):
            return
        sampler_self = self

        def _save_with_jarvis_state(fname: str, *args: Any, **kwargs: Any) -> None:
            # Strip before dynesty pickles ``self`` inside save_sampler.
            sampler_self._strip_dynesty_save_hook(engine)
            try:
                type(engine).save(engine, fname, *args, **kwargs)
            finally:
                sampler_self._install_dynesty_save_hook()
            try:
                blob = sampler_self._pickle_native_sampler(engine)
                if blob is not None:
                    sampler_self._native_sampler_blob = blob
            except Exception as exc:
                sampler_self._logger.warning(
                    "%s: pickle after dynesty save failed: %s", sampler_self.method, exc
                )
            cb = getattr(sampler_self, "_save_checkpoint_callback", None)
            if cb is not None:
                try:
                    sampler_self._last_safe_state = None
                    cb(reason="dynesty_engine_checkpoint")
                except Exception as exc:
                    sampler_self._logger.warning(
                        "%s: state.pkl write after dynesty save failed: %s",
                        sampler_self.method,
                        exc,
                    )

        engine.save = _save_with_jarvis_state  # type: ignore[method-assign]
        self._wrapped_dynesty_save = True

    def checkpoint_runtime_state(self, *, safe: bool) -> dict[str, Any]:
        """Always snapshot the live dynesty engine (not a stale barrier copy).

        Nested sampling does not use FeedbackSampler generation barriers.  The
        default ``CheckpointedSampler.checkpoint_runtime_state`` may return the
        post-configure snapshot when ``safe=False``, freezing pre-run state and
        making interrupt checkpoints useless.  Always re-export.
        """
        state = self.export_runtime_state()
        # Avoid deepcopy of multi-MB blobs twice; export already returns a new dict.
        self._last_safe_state = state
        return dict(state)

    def configure_checkpoint(
        self,
        *,
        checkpoint_file: str,
        save_callback,
    ) -> None:
        """Wire Jarvis state.pkl writes; mid-run cadence = dynesty checkpoint_every."""
        from jarvishep2.sampling.runtime_checkpoint import (
            CHECKPOINT_HEARTBEAT_SEC,
            CheckpointHeartbeat,
        )

        self._checkpoint_file = str(checkpoint_file)
        self._save_checkpoint_callback = save_callback
        # Capture current engine (restored or empty) as the first safe snapshot.
        try:
            self._last_safe_state = self.export_runtime_state()
        except Exception:
            self._last_safe_state = {"method": self.method, "finished": self._finished}
        # Heartbeat only requests a flag for enum samplers.  Dynesty never polls
        # it; keep a long interval no-op so CheckpointedSampler wiring stays happy.
        if self._checkpoint_heartbeat is not None:
            self._checkpoint_heartbeat.stop()
        interval = float(
            getattr(self, "_checkpoint_heartbeat_sec", None)
            or getattr(self, "_checkpoint_every_sec", None)
            or CHECKPOINT_HEARTBEAT_SEC
        )
        self._checkpoint_heartbeat = CheckpointHeartbeat(
            interval_sec=max(5.0, interval),
            save_callback=lambda reason="checkpoint_heartbeat": None,
        )
        self._checkpoint_heartbeat.start()
        self._logger.info(
            "%s checkpoint: nested_engine.pkl every %.0fs (dynesty); "
            "state.pkl on each engine save + interrupt/finish",
            self.method,
            float(self._checkpoint_every_sec),
        )

    def _load_native_from_blob(self, blob: bytes) -> Any | None:
        try:
            return pickle.loads(blob)
        except Exception as exc:
            self._logger.warning(
                "%s: failed to unpickle native engine blob: %s", self.method, exc
            )
            return None

    def _load_native_from_engine_file(self, path: str, pool: Any = None) -> Any | None:
        if not path or not os.path.isfile(path):
            return None
        try:
            from jarvishep2.sampling.Source.Dynesty.py.dynesty.utils import (
                restore_sampler,
            )

            return restore_sampler(path, pool=pool)
        except Exception as exc:
            self._logger.warning(
                "%s: failed to restore engine from %s: %s", self.method, path, exc
            )
            return None

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
        """Load a previously interrupted dynesty/MultiNest engine if available.

        Only runs when :attr:`_resume_engine_requested` is set by
        ``import_runtime_state`` / ``arm_engine_resume`` (Core ``--resume``).
        A leftover ``nested_engine.pkl`` alone must **not** hijack a fresh start.

        Preference order (newest / most durable first):

        1. Side-file ``nested_engine.pkl`` — written by dynesty's own timer
           between iterations (survives kill -9 when the last timed save completed)
        2. Live ``self._sampler`` already unpickled by ``import_runtime_state``
        3. Embedded pickle blob in the runtime ``state.pkl``
        """
        if not self._resume_engine_requested:
            return None, False
        if self._finished:
            return None, False

        bridge = NestedRedisLogL(pool)
        self._loglike_bridge = bridge

        candidates: list[tuple[str, Any | None]] = []
        engine_path = self._engine_path()
        if engine_path and os.path.isfile(engine_path):
            candidates.append(("engine_file", None))  # load lazily
        if self._sampler is not None:
            candidates.append(("live", self._sampler))
        if self._native_sampler_blob:
            candidates.append(("blob", None))

        for source, obj in candidates:
            sampler = obj
            if source == "engine_file":
                sampler = self._load_native_from_engine_file(engine_path, pool=pool)
            elif source == "blob":
                sampler = self._load_native_from_blob(self._native_sampler_blob or b"")
            if sampler is None:
                continue
            if not self._engine_matches_method(sampler):
                expected = (
                    "DynamicNestedSampler"
                    if self.default_dynamic
                    else "static NestedSampler"
                )
                self._logger.warning(
                    "%s resume: skipping engine from %s — type %s is not %s "
                    "(Method selects the engine; do not mix Dynesty/MultiNest "
                    "checkpoints under the same scan name)",
                    self.method,
                    source,
                    type(sampler).__name__,
                    expected,
                )
                continue
            if self._engine_looks_finished(sampler):
                self._logger.info(
                    "%s resume: engine from %s is already complete "
                    "(RUN_DONE / added_live); not re-entering run_nested",
                    self.method,
                    source,
                )
                self._reattach_native_runtime(
                    sampler, pool=pool, bridge=bridge, logger=logger
                )
                self._sampler = sampler
                self._finished = True
                if self._summary is None:
                    try:
                        self._summary = self._build_summary()
                    except Exception:
                        pass
                return sampler, False
            self._reattach_native_runtime(
                sampler, pool=pool, bridge=bridge, logger=logger
            )
            self._logger.info(
                "%s resume: restored native dynesty engine from %s (it=%s ncall=%s)",
                self.method,
                source,
                getattr(sampler, "it", "?"),
                getattr(sampler, "ncall", "?"),
            )
            return sampler, True
        return None, False

    def run_adaptive(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Run nested sampling with a per-generation timeout setting.

        Supports mid-run resume via dynesty native pickle
        (``nested_engine.pkl`` + runtime ``state.pkl`` blob) and
        ``run_nested(resume=True)``.
        """
        if timeout is not None:
            generation_timeout = timeout
        timeout = generation_timeout
        self._require_redis(f"{self.method}.run_adaptive")
        self._ensure_seed_sequence()
        if self._finished and self._sampler is not None and not self._resume_engine_requested:
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
        )
        # Resume call_index so any non-uuid fallback UUIDs stay deterministic.
        pool._call_index = max(0, int(self._pool_call_index or 0))

        resume = False
        wanted_resume = bool(self._resume_engine_requested)
        restored, resume = self._resolve_resume_engine(pool=pool, logger=inner_logger)
        # Consume the one-shot resume request after the resolve attempt so a
        # later accidental re-entry cannot silently re-load a finished engine.
        self._resume_engine_requested = False
        if self._finished and restored is not None:
            return int((self._summary or {}).get("ncall") or 0)
        if restored is not None:
            self._sampler = restored
            self._wrapped_dynesty_save = False  # restored object has stock .save
            self._logger.warning(
                "%s RESUME: continuing dynesty engine it=%s ncall=%s state=%s "
                "(run_nested resume=%s)",
                self.method,
                getattr(restored, "it", "?"),
                getattr(restored, "ncall", "?"),
                getattr(restored, "internal_state", "?"),
                resume,
            )
        else:
            if wanted_resume:
                self._logger.warning(
                    "%s --resume was requested but no native engine was restored "
                    "(engine_file=%s exists=%s, blob=%s). Starting a FRESH nested "
                    "run — progress will restart from iter≈0. Check "
                    "checkpoints/<scan>/%s/nested_engine.pkl and state.pkl.",
                    self.method,
                    self._engine_path(),
                    os.path.isfile(self._engine_path()),
                    "yes" if self._native_sampler_blob else "no",
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

            self._logger.warning(
                "%s START FRESH %s nlive=%d ndim=%d dynamic=%s bound=%s sample=%s",
                self.method,
                sampler_cls.__name__,
                self._nlive,
                self._dim,
                self._use_dynamic,
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

        # Mid-run: dynesty.save → nested_engine.pkl + Jarvis state.pkl together.
        self._install_dynesty_save_hook()

        run_kwargs = filter_kwargs_to_callable(
            self._run_nested_kwargs,
            self._sampler.run_nested,
            context="run_nested",
            logger=self._logger,
        )
        # HEP injects checkpoint path/interval from EnvReqs.V2 — never from Sampling YAML.
        run_kwargs["checkpoint_file"] = self._engine_path()
        run_kwargs["checkpoint_every"] = float(self._checkpoint_every_sec)
        # Dynesty's save_sampler writes ``path.tmp`` then renames — parent dir
        # must exist (stock dynesty does not create it).
        ckpt = run_kwargs.get("checkpoint_file")
        if ckpt:
            parent = os.path.dirname(os.path.abspath(str(ckpt)))
            if parent:
                os.makedirs(parent, exist_ok=True)
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
        try:
            self._sampler.run_nested(**run_kwargs)
            self._finished = True
        except BaseException:
            # Interrupt / crash path: persist engine so --resume can continue.
            try:
                self._save_engine_side_file()
                blob = self._pickle_native_sampler()
                if blob is not None:
                    self._native_sampler_blob = blob
            except Exception as save_exc:
                self._logger.warning(
                    "%s: emergency engine checkpoint failed: %s",
                    self.method,
                    save_exc,
                )
            raise
        finally:
            try:
                self._pool_call_index = int(getattr(pool, "_call_index", 0) or 0)
            except Exception:
                pass

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
            from jarvishep2.sampling.diagnostics_export import export_nested_diagnostics

            written = export_nested_diagnostics(
                self._task_result_dir(), summary=self._summary
            )
            if written and self._summary is not None:
                self._summary.setdefault("diagnostics", {}).update(written)
            for key, path in written.items():
                self._logger.info("%s diagnostics %s → %s", self.method, key, path)
        except Exception as exc:
            self._logger.warning("failed to write nested sampler_summary: %s", exc)
        # Final engine snapshot (completed) + Jarvis state.pkl.
        try:
            self._save_engine_side_file()
            blob = self._pickle_native_sampler()
            if blob is not None:
                self._native_sampler_blob = blob
        except Exception as exc:
            self._logger.warning("%s: final engine snapshot failed: %s", self.method, exc)
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
        """Export Jarvis feedback fields + pickled dynesty engine (live points).

        The native engine is stored two ways:

        * ``native_sampler_blob`` — stdlib pickle bytes inside ``state.pkl``
        * ``nested_engine.pkl`` side file — dynesty's own ``save()`` format
          (also written periodically by ``run_nested(checkpoint_file=…)``)

        Redis / logger / pool are never pickled; they are reattached on resume.
        """
        # Prefer a fresh snapshot of the live engine when present (interrupt path).
        if self._sampler is not None:
            blob = self._pickle_native_sampler()
            if blob is not None:
                self._native_sampler_blob = blob
            try:
                self._save_engine_side_file()
            except Exception as exc:
                self._logger.warning(
                    "%s: engine side-file write during export failed: %s",
                    self.method,
                    exc,
                )

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
        engine_path = self._engine_path()
        state.update(
            {
                "method": self.method,
                "finished": self._finished,
                "summary": self._summary,
                "nlive": self._nlive,
                "dlogz": self._dlogz,
                "rseed": self._rseed,
                "use_dynamic": self._use_dynamic,
                "constructor_kwargs": safe_ctor,
                "run_nested_kwargs": safe_run,
                "dynesty_csv_path": self._dynesty_csv_path,
                "call_index": int(self._pool_call_index or 0),
                "pool_call_index": int(self._pool_call_index or 0),
                "engine_checkpoint_path": engine_path,
                "native_sampler_blob": self._native_sampler_blob,
                "native_sampler_format": "pickle",
                "nested_engine_filename": NESTED_ENGINE_FILENAME,
            }
        )
        return state

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        """Restore feedback metadata + queue native engine for ``run_adaptive``."""
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
        run = state.get("run_nested_kwargs")
        if isinstance(run, Mapping):
            self._run_nested_kwargs = {
                str(k): v for k, v in run.items() if str(k) != "print_func"
            }
        path = state.get("dynesty_csv_path")
        self._dynesty_csv_path = str(path) if path else None
        self._pool_call_index = int(
            state.get("pool_call_index", state.get("call_index", 0)) or 0
        )
        engine_path = state.get("engine_checkpoint_path")
        if engine_path:
            self._engine_checkpoint_path = str(engine_path)
        elif not self._engine_checkpoint_path:
            self._engine_checkpoint_path = self._default_engine_checkpoint_path()

        blob = state.get("native_sampler_blob")
        if isinstance(blob, (bytes, bytearray)):
            self._native_sampler_blob = bytes(blob)
            # Eagerly unpickle so export/interrupt have a live object if needed;
            # pool reattachment happens in run_adaptive.
            if not self._finished:
                native = self._load_native_from_blob(self._native_sampler_blob)
                if native is not None:
                    self._sampler = native
        elif blob is not None:
            # Legacy: native object embedded directly (already unpickled by outer load).
            try:
                if hasattr(blob, "run_nested") or hasattr(blob, "results"):
                    self._sampler = blob
                    self._native_sampler_blob = self._pickle_native_sampler(blob)
            except Exception:
                self._logger.warning(
                    "%s: ignoring unreadable native_sampler in checkpoint",
                    self.method,
                )


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
