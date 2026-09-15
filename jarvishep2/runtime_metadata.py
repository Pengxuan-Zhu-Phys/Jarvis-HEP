"""Runtime scan metadata published through the scan's own Redis instance."""

from __future__ import annotations

import json
import os
import tempfile
import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.sampler_catalog import family_of


RUNTIME_METADATA_KEY = "hep:runtime:metadata_path"


def _short_scalar(value: Any) -> Any | None:
    """Project one safe sampler knob into runtime metadata."""
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, str):
        return value[:128]
    if isinstance(value, (list, tuple)) and len(value) <= 16:
        projected = [_short_scalar(item) for item in value]
        if all(item is not None for item in projected):
            return projected
    return None


def sampler_metadata_from_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Return a non-secret, allowlisted description of ``Sampling``."""
    raw_sampling = config.get("Sampling")
    sampling = dict(raw_sampling) if isinstance(raw_sampling, Mapping) else {}
    method = str(sampling.get("Method") or "Unknown").strip() or "Unknown"
    raw_bounds = sampling.get("Bounds")
    bounds = dict(raw_bounds) if isinstance(raw_bounds, Mapping) else {}
    variables = sampling.get("Variables")
    dimensions = len(variables) if isinstance(variables, list) else 0

    allowlist: dict[str, tuple[str, ...]] = {
        "Random": ("point_number", "seed"),
        "Bridson": ("radius", "max_attempt", "seed"),
        "AdaptiveBridson": (
            "initial_radius",
            "min_radius",
            "max_generations",
            "seed",
        ),
        "MCMC": ("num_chains", "num_iters", "proposal_scale", "seed"),
        "ToyMCMC": ("num_chains", "num_iters", "proposal_scale", "seed"),
        "AMMCMC": (
            "num_chains",
            "num_iters",
            "adapt_enabled",
            "adapt_start_iter",
            "adapt_window",
            "seed",
        ),
        "DRAM": (
            "num_chains",
            "num_iters",
            "dr_steps",
            "dr_scale_factors",
            "seed",
        ),
        "EnsembleMCMC": ("num_chains", "num_iters", "stretch_a", "seed"),
        "DEMCMC": ("num_chains", "num_iters", "de_gamma", "seed"),
        "PTMCMC": (
            "num_chains",
            "num_iters",
            "temperature_ladder",
            "exchange_interval",
            "seed",
        ),
        "PTEnsemble": (
            "num_chains",
            "num_iters",
            "temperature_ladder",
            "exchange_interval",
            "stretch_a",
            "seed",
        ),
        "Dynesty": ("nlive", "dlogz", "dlogz_init", "seed"),
        "MultiNest": ("nlive", "dlogz", "dlogz_init", "seed"),
    }
    knobs: dict[str, Any] = {}
    for key in allowlist.get(method, ()):
        if key not in bounds:
            continue
        projected = _short_scalar(bounds[key])
        if projected is not None:
            knobs[key] = projected

    if method == "CSV":
        path = str(bounds.get("path") or "").strip()
        if path:
            knobs["basename"] = os.path.basename(path)[:128]
        uuid_column = _short_scalar(bounds.get("uuid_column"))
        if uuid_column is not None:
            knobs["uuid_column"] = uuid_column
    elif method == "Grid":
        shape: list[int] = []
        for variable in variables if isinstance(variables, list) else []:
            if not isinstance(variable, Mapping):
                shape = []
                break
            distribution = variable.get("distribution")
            if not isinstance(distribution, Mapping):
                shape = []
                break
            parameters = distribution.get("parameters")
            if not isinstance(parameters, Mapping) or "num" not in parameters:
                shape = []
                break
            try:
                shape.append(max(1, int(parameters["num"])))
            except (TypeError, ValueError):
                shape = []
                break
        if shape:
            knobs["shape"] = shape
            total = 1
            for count in shape:
                total *= count
            knobs["total"] = total

    return {
        "method": method,
        "family": family_of(method),
        "dimensions": dimensions,
        "config": knobs,
    }


def _calculator_pool_catalog(
    config: Mapping[str, Any],
) -> tuple[dict[str, int], dict[str, dict[str, Any]]]:
    """Return monitor-safe physical pool metadata from a task card.

    Calculator occupancy cannot safely be reconstructed from the stale
    ``hep:calculator:status`` cache while the monitor is not watching.  Keep
    the static pool catalogue with the scan metadata instead.  This helper is
    intentionally best-effort: failing to describe a pool must never prevent a
    scan from publishing its identity record.
    """
    try:
        from jarvishep2.calculator_modes import (
            expand_calculator_modes,
            shared_mode_groups,
        )
        from jarvishep2.calculator_pools import resolve_calculator_pools

        calculators = config.get("Calculators")
        block = dict(calculators) if isinstance(calculators, Mapping) else {}
        raw_modules = block.get("Modules") or []
        modules = (
            expand_calculator_modes(raw_modules)
            if isinstance(raw_modules, list)
            else []
        )
        worker_config: dict[str, Any] = {"calculator_modules": modules}
        pools = block.get("Pools") or block.get("pools")
        if isinstance(pools, Mapping):
            worker_config["calculator_pools"] = dict(pools)
        if "make_parallel" in block:
            worker_config["calculator_make_parallel"] = block["make_parallel"]

        physical = resolve_calculator_pools(worker_config)
        groups = shared_mode_groups(modules)
    except (ImportError, TypeError, ValueError):
        return {}, {}

    normal = {
        name: int(slots)
        for name, slots in physical.items()
        if name not in groups and int(slots) > 0
    }
    shared = {
        parent: {"modes": list(modes), "n": int(physical[parent])}
        for parent, modes in groups.items()
        if parent in physical and int(physical[parent]) > 0
    }
    return normal, shared


def calculator_pools_from_metadata(
    metadata: Mapping[str, Any],
) -> tuple[list[str] | None, dict[str, int] | None]:
    """Extract the physical calculator catalogue from a runtime record.

    Older scans do not carry these optional fields, so ``(None, None)`` keeps
    the reader's existing discovery fallback intact.
    """
    slots: dict[str, int] = {}
    pools = metadata.get("calculator_pools")
    if isinstance(pools, Mapping):
        for name, value in pools.items():
            text = str(name or "").strip()
            try:
                count = int(value)
            except (TypeError, ValueError):
                continue
            if text and count > 0:
                slots[text] = count
    shared = metadata.get("calculator_shared")
    if isinstance(shared, Mapping):
        for parent, spec in shared.items():
            text = str(parent or "").strip()
            if not text or not isinstance(spec, Mapping):
                continue
            try:
                count = int(spec.get("n"))
            except (TypeError, ValueError):
                continue
            if count > 0:
                slots.setdefault(text, count)
    names = list(slots)
    return (names, slots) if names else (None, None)


def write_scan_metadata(*, config: Mapping[str, Any], info: Mapping[str, Any], redis: Mapping[str, Any]) -> str:
    """Write the non-secret scan identity record and return its absolute path."""
    task_result_dir = os.path.abspath(str(info.get("task_result_dir") or config.get("task_result_dir") or os.getcwd()))
    directory = os.path.join(task_result_dir, ".jarvis2")
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, "runtime.json")
    calculator_pools, calculator_shared = _calculator_pool_catalog(config)
    payload: dict[str, Any] = {
        "schema": 1,
        "written_at": time.time(),
        "scan_name": str(info.get("scan_name") or config.get("scan_name") or "scan"),
        "task_yaml": str(config.get("task_yaml") or ""),
        "project_root": str(config.get("project_root") or config.get("task_root") or ""),
        "task_result_dir": task_result_dir,
        "redis": {"host": str(redis["host"]), "port": int(redis["port"]), "db": int(redis["db"])},
        "control_pid": os.getpid(),
        "calculator_pools": calculator_pools,
        "calculator_shared": calculator_shared,
        "sampler": sampler_metadata_from_config(config),
    }
    fd, temporary = tempfile.mkstemp(prefix=".runtime.", suffix=".json", dir=directory)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise
    return path


def read_scan_metadata(path: str, *, redis: Mapping[str, Any], expected_scan: str) -> dict[str, Any] | None:
    """Read and validate the record advertised by Redis before using it."""
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, ValueError, TypeError):
        return None
    if not isinstance(payload, dict) or payload.get("schema") != 1:
        return None
    if str(payload.get("scan_name") or "") != str(expected_scan):
        return None
    recorded = payload.get("redis")
    if not isinstance(recorded, Mapping):
        return None
    if (str(recorded.get("host")), int(recorded.get("port", -1)), int(recorded.get("db", -1))) != (
        str(redis.get("host")), int(redis.get("port", -2)), int(redis.get("db", -2))
    ):
        return None
    return payload


__all__ = [
    "RUNTIME_METADATA_KEY",
    "calculator_pools_from_metadata",
    "read_scan_metadata",
    "sampler_metadata_from_config",
    "write_scan_metadata",
]
