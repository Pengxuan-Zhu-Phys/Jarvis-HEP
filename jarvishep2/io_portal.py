#!/usr/bin/env python3
"""Thin bridge from Jarvis-HEP V2 calculator I/O to Jarvis-Portal.

Format adapters live in Portal. This module owns expression evaluation, registry
caching, context construction, and HEP-facing error messages.
"""

from __future__ import annotations

import asyncio
import math
from collections.abc import Mapping
from typing import Any

import numpy as np
import sympy as sp

_REGISTRY = None
_PORTAL_IMPORT_ERROR: Exception | None = None


class UnsupportedIOTypeError(ValueError):
    """Raised when a calculator IO ``type`` has no registered Portal adapter."""


def _import_portal():
    global _PORTAL_IMPORT_ERROR
    try:
        from jarvis_portal import (  # type: ignore import-not-found
            IOContext,
            MissingAdapterError,
            create_entry_point_registry,
        )
    except ImportError as exc:
        _PORTAL_IMPORT_ERROR = exc
        raise ImportError(
            "Calculator I/O requires Jarvis-HEP-Portal. "
            "Install it with `pip install Jarvis-HEP-Portal` "
            "(or `pip install -e ../Jarvis-Portal` for local development)."
        ) from exc
    return IOContext, MissingAdapterError, create_entry_point_registry


def get_io_registry():
    """Return the process-local Portal registry (builtins + entry points)."""
    global _REGISTRY
    if _REGISTRY is not None:
        return _REGISTRY
    _, _, create_entry_point_registry = _import_portal()
    _REGISTRY = create_entry_point_registry()
    return _REGISTRY


def reset_io_registry_for_tests() -> None:
    """Drop the cached registry (test isolation only)."""
    global _REGISTRY
    _REGISTRY = None


def available_io_formats(direction: str | None = None) -> list[str]:
    """List format names registered for the given direction (or all)."""
    return get_io_registry().available_formats(direction)


def evaluate_io_expression(expression: str, values: Mapping[str, Any]) -> Any:
    """Evaluate a Dump expression with HEP-owned sympy semantics (``Pi`` / ``pi``)."""
    text = str(expression).strip()
    if not text:
        raise ValueError("IO expression must not be empty.")
    expr = sp.sympify(text, locals={"Pi": sp.pi, "pi": sp.pi})
    substitutions: dict[Any, float] = {}
    missing: list[str] = []
    for symbol in expr.free_symbols:
        name = str(symbol)
        if name in {"Pi", "pi"}:
            substitutions[symbol] = math.pi
            continue
        numeric = _coerce_numeric_param(values.get(name))
        if numeric is None:
            missing.append(name)
            continue
        substitutions[symbol] = numeric
    if missing:
        raise ValueError(f"Dump expression '{text}' misses parameters: {sorted(missing)}")
    return float(expr.evalf(subs=substitutions))


def _coerce_numeric_param(value: Any) -> float | None:
    if isinstance(value, bool):
        return float(int(value))
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(np.asarray(value).item() if isinstance(value, np.generic) else value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return None
    return None


def build_io_context(
    *,
    sample_info: Mapping[str, Any] | None = None,
    pack_id: str | None = None,
    module: str | None = None,
    project_root: str | None = None,
    evaluate_expression=None,
    runtime_values: Mapping[str, Any] | None = None,
    logger: Any = None,
):
    """Construct a Portal ``IOContext`` from Worker/sample fields."""
    IOContext, _, _ = _import_portal()
    info = dict(sample_info or {})
    sample_uuid = info.get("uuid")
    sample_save_dir = info.get("save_dir")
    if sample_save_dir is not None:
        sample_save_dir = str(sample_save_dir)
    if logger is None and "logger" in info:
        logger = info.get("logger")
    evaluator = evaluate_expression if evaluate_expression is not None else evaluate_io_expression
    return IOContext(
        logger=logger,
        sample_uuid=str(sample_uuid) if sample_uuid is not None else None,
        pack_id=str(pack_id) if pack_id is not None else None,
        sample_save_dir=sample_save_dir,
        project_root=project_root,
        module=module,
        evaluate_expression=evaluator,
        runtime_values=dict(runtime_values or {}),
    )


def _normalize_type(spec: Mapping[str, Any]) -> str:
    raw = str(spec.get("type", "")).strip()
    if not raw:
        raise ValueError("Calculator IO spec requires a non-empty 'type' field.")
    return raw


def _adapter_spec(spec: Mapping[str, Any], *, path: str) -> dict[str, Any]:
    payload = dict(spec)
    payload["path"] = str(path)
    return payload


async def write_io_input(
    spec: Mapping[str, Any],
    data: Mapping[str, Any],
    *,
    context,
    path: str | None = None,
) -> dict[str, Any]:
    """Write one input file via Portal and return adapter observables."""
    format_name = _normalize_type(spec)
    resolved_path = str(path if path is not None else spec.get("path", ""))
    if not resolved_path:
        raise ValueError(f"IO input spec '{spec.get('name', '')}' requires a path.")
    registry = get_io_registry()
    _, MissingAdapterError, _ = _import_portal()
    try:
        adapter = registry.get(format_name, "input")
    except MissingAdapterError as exc:
        available = ", ".join(registry.available_formats("input")) or "none"
        raise UnsupportedIOTypeError(
            f"unknown IO type {format_name!r} for input; registered: {available}"
        ) from exc
    return await adapter.write_input(context, _adapter_spec(spec, path=resolved_path), dict(data))


async def read_io_output(
    spec: Mapping[str, Any],
    *,
    context,
    path: str | None = None,
) -> dict[str, Any]:
    """Read one output file via Portal and return observables."""
    format_name = _normalize_type(spec)
    resolved_path = str(path if path is not None else spec.get("path", ""))
    if not resolved_path:
        raise ValueError(f"IO output spec '{spec.get('name', '')}' requires a path.")
    registry = get_io_registry()
    _, MissingAdapterError, _ = _import_portal()
    try:
        adapter = registry.get(format_name, "output")
    except MissingAdapterError as exc:
        available = ", ".join(registry.available_formats("output")) or "none"
        raise UnsupportedIOTypeError(
            f"unknown IO type {format_name!r} for output; registered: {available}"
        ) from exc
    return await adapter.read_output(context, _adapter_spec(spec, path=resolved_path))


def _run_coro_sync(coro: Any) -> Any:
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)
    raise RuntimeError(
        "Synchronous Portal I/O cannot run inside an active event loop; "
        "use the async write_io_input/read_io_output helpers instead."
    )


def write_io_input_sync(
    spec: Mapping[str, Any],
    data: Mapping[str, Any],
    *,
    context,
    path: str | None = None,
) -> dict[str, Any]:
    return _run_coro_sync(write_io_input(spec, data, context=context, path=path))


def read_io_output_sync(
    spec: Mapping[str, Any],
    *,
    context,
    path: str | None = None,
) -> dict[str, Any]:
    return _run_coro_sync(read_io_output(spec, context=context, path=path))


__all__ = [
    "UnsupportedIOTypeError",
    "available_io_formats",
    "build_io_context",
    "evaluate_io_expression",
    "get_io_registry",
    "read_io_output",
    "read_io_output_sync",
    "reset_io_registry_for_tests",
    "write_io_input",
    "write_io_input_sync",
]
