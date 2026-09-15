"""Bounded, read-only sampler telemetry for the Monitor Core board."""

from __future__ import annotations

import json
import math
import time
from collections.abc import Mapping, Sequence
from typing import Any


SAMPLER_STATUS_SCHEMA = 1
SAMPLER_STATUS_MAX_BYTES = 4096


def _integer(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError, OverflowError):
        return None


def _number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return result if math.isfinite(result) else None


def _length(value: Any) -> int | None:
    try:
        return len(value)
    except (TypeError, ValueError):
        return None


def _state(sampler: Any) -> str:
    if bool(getattr(sampler, "_converged", False)):
        return "converged"
    if bool(getattr(sampler, "_finished", False)):
        return "complete"
    if sampler is None:
        return "initializing"
    return "running"


def _progress(kind: str, current: Any, target: Any) -> dict[str, Any]:
    return {
        "kind": kind,
        "current": _integer(current),
        "target": _integer(target),
    }


def _mcmc_status(sampler: Any, method: str) -> tuple[dict[str, Any], dict[str, Any], str]:
    registry = getattr(sampler, "_registry", None)
    iterations: list[int] = []
    if registry is not None:
        try:
            iterations = [
                max(0, int(chain.engine.iterations))
                for chain in registry.all()
            ]
        except (AttributeError, TypeError, ValueError):
            iterations = []
    state = _state(sampler) if iterations else "initializing"
    kind = "generation" if method in {"EnsembleMCMC", "DEMCMC", "PTEnsemble"} else "iteration"
    progress = _progress(
        kind,
        min(iterations) if iterations else None,
        getattr(sampler, "_niters", None),
    )
    proposed = _integer(getattr(sampler, "_total_proposed", 0)) or 0
    accepted = _integer(getattr(sampler, "_total_accepted", 0)) or 0
    metrics: dict[str, Any] = {
        "chains": _integer(getattr(sampler, "_nchains", None)),
        "accept_rate": accepted / proposed if proposed > 0 else None,
    }
    if method == "ToyMCMC":
        scale = getattr(sampler, "_proposal_scales", None)
        if isinstance(scale, Sequence) and not isinstance(scale, (str, bytes)):
            scale = scale[0] if scale else None
        metrics["proposal_scale"] = _number(scale)
    elif method == "AMMCMC":
        if not bool(getattr(sampler, "_adapt_enabled", False)):
            metrics["adapt_state"] = "disabled"
        else:
            start = _integer(getattr(sampler, "_adapt_start_iter", 0)) or 0
            floor = min(iterations) if iterations else 0
            metrics["adapt_state"] = "active" if floor >= start else "warmup"
    elif method == "DRAM":
        base_steps = sum(iterations)
        retries = max(0, proposed - base_steps)
        metrics["retry_rate"] = retries / proposed if proposed > 0 else None
    elif method == "EnsembleMCMC":
        metrics["walkers"] = metrics.pop("chains")
        metrics["stretch_a"] = _number(getattr(sampler, "_stretch_a", None))
    elif method == "DEMCMC":
        metrics["walkers"] = metrics.pop("chains")
        metrics["de_gamma"] = _number(getattr(sampler, "_de_gamma", None))
    elif method in {"PTMCMC", "PTEnsemble"}:
        metrics["swap_accepts"] = _integer(getattr(sampler, "_swap_accepts", 0))
        metrics["swap_attempts"] = _integer(getattr(sampler, "_swap_attempts", 0))
    return progress, metrics, state


def _nested_status(sampler: Any) -> tuple[dict[str, Any], dict[str, Any], str]:
    native = getattr(sampler, "_sampler", None)
    niter = _integer(getattr(native, "it", None)) if native is not None else None
    # Only read native scalar counters. Accessing ``results.ncall`` here can
    # walk a run-sized array on every heartbeat, so unavailable stays honest.
    ncall = _integer(getattr(native, "ncall", None)) if native is not None else None
    efficiency = niter / ncall if niter is not None and ncall and ncall > 0 else None
    metrics = {
        "nlive": _integer(getattr(sampler, "_nlive", None)),
        "niter": niter,
        "ncall": ncall,
        "efficiency": efficiency,
        "dlogz_target": _number(getattr(sampler, "_dlogz", None)),
    }
    state = _state(sampler) if native is not None else "initializing"
    return {"kind": "evidence", "current": None, "target": None}, metrics, state


def sampler_status_snapshot(sampler: Any) -> dict[str, Any]:
    """Read cheap scalar attributes without advancing or summarizing a sampler."""
    method = str(getattr(sampler, "method", None) or type(sampler).__name__)
    state = _state(sampler)
    metrics: dict[str, Any] = {}
    if method == "Random":
        progress = _progress("finite", getattr(sampler, "_index", None), getattr(sampler, "_maxp", None))
        metrics.update(
            accepted=_integer(getattr(sampler, "_accepted_index", None)),
            seed=_integer(getattr(sampler, "_seed", None)),
        )
    elif method in {"Grid", "Bridson"}:
        points = getattr(sampler, "_P", None)
        target = _length(points)
        if target is None:
            info = getattr(sampler, "info", None)
            target = _integer(info.get("NSamples")) if isinstance(info, Mapping) else None
        progress = _progress("finite", getattr(sampler, "_index", None), target)
        if method == "Bridson":
            metrics.update(
                radius=_number(getattr(sampler, "_radius", None)),
                max_attempt=_integer(getattr(sampler, "_k", None)),
            )
    elif method == "CSV":
        progress = _progress(
            "finite",
            getattr(sampler, "_runtime_csv_cursor", None),
            getattr(sampler, "_source_row_total", None),
        )
        metrics["accepted"] = _integer(getattr(sampler, "_accepted_index", None))
    elif method == "AdaptiveBridson":
        progress = _progress(
            "generation",
            getattr(sampler, "_generation", None),
            getattr(sampler, "_max_generations", None),
        )
        metrics.update(
            radius=_number(getattr(sampler, "_radius", None)),
            core=_length(getattr(sampler, "_live_core_indices", None)),
            open=_integer(getattr(sampler, "_open_brackets", None)),
        )
        if getattr(sampler, "_stop_reason", None):
            state = "converged" if bool(getattr(sampler, "_converged", False)) else "partial"
    elif method in {
        "MCMC",
        "ToyMCMC",
        "AMMCMC",
        "DRAM",
        "EnsembleMCMC",
        "DEMCMC",
        "PTMCMC",
        "PTEnsemble",
    }:
        progress, metrics, state = _mcmc_status(sampler, method)
    elif method in {"Dynesty", "MultiNest"}:
        progress, metrics, state = _nested_status(sampler)
    else:
        progress = {"kind": "none", "current": None, "target": None}

    return {
        "schema": SAMPLER_STATUS_SCHEMA,
        "ts": time.time(),
        "method": method[:64],
        "state": state[:32],
        "progress": progress,
        "metrics": {key: value for key, value in metrics.items() if value is not None},
    }


def sampler_status_json(sampler: Any) -> str:
    """Serialize one status field, enforcing the Core-board 4 KiB contract."""
    try:
        payload = sampler_status_snapshot(sampler)
    except Exception:
        # Observability must never prevent the existing Core lease board from
        # refreshing. A minimal payload is still enough to remain honest.
        try:
            method = str(getattr(sampler, "method", "Unknown"))[:64]
        except Exception:
            method = "Unknown"
        payload = {
            "schema": SAMPLER_STATUS_SCHEMA,
            "ts": time.time(),
            "method": method,
            "state": "unavailable",
            "progress": {"kind": "none", "current": None, "target": None},
            "metrics": {},
        }
    encoded = json.dumps(payload, separators=(",", ":"), sort_keys=True)
    if len(encoded.encode("utf-8")) <= SAMPLER_STATUS_MAX_BYTES:
        return encoded
    fallback = {
        "schema": SAMPLER_STATUS_SCHEMA,
        "ts": payload["ts"],
        "method": payload["method"],
        "state": payload["state"],
        "progress": {"kind": "none", "current": None, "target": None},
        "metrics": {},
    }
    return json.dumps(fallback, separators=(",", ":"), sort_keys=True)


__all__ = [
    "SAMPLER_STATUS_MAX_BYTES",
    "SAMPLER_STATUS_SCHEMA",
    "sampler_status_json",
    "sampler_status_snapshot",
]
