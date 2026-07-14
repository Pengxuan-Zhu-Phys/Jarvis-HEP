#!/usr/bin/env python3
"""Thin bridge from Jarvis-HEP V2 calculator I/O to Jarvis-Portal.

Format adapters live in Portal. This module owns expression evaluation, registry
caching, context construction, and HEP-facing error messages.
"""

from __future__ import annotations

import asyncio
from collections.abc import Mapping
from typing import Any

import numpy as np

from jarvishep2.expression import (
    ExpressionContext,
    MissingExpressionVariablesError,
)

_REGISTRY = None
_PORTAL_IMPORT_ERROR: Exception | None = None
_IO_EXPRESSIONS: ExpressionContext | None = None
_IO_OPERAS_LOADED = False


class UnsupportedIOTypeError(ValueError):
    """Raised when a calculator IO ``type`` has no registered Portal adapter."""


def _import_portal():
    """Import the **V2** Portal product surface (not V1 top-level factories)."""
    global _PORTAL_IMPORT_ERROR
    try:
        # Prefer explicit V2 interface package so V1/V2 surfaces never share factories.
        from jarvis_portal.context import IOContext  # type: ignore import-not-found
        from jarvis_portal.core.registry import (  # type: ignore import-not-found
            MissingAdapterError,
        )
        from jarvis_portal.v2 import (  # type: ignore import-not-found
            create_entry_point_registry,
        )
    except ImportError as exc:
        # Older Portal without v2 package: fall back to top-level (V1 surface).
        try:
            from jarvis_portal import (  # type: ignore import-not-found
                IOContext,
                MissingAdapterError,
                create_entry_point_registry,
            )
        except ImportError as inner:
            _PORTAL_IMPORT_ERROR = inner
            raise ImportError(
                "Calculator I/O requires Jarvis-HEP-Portal>=1.4.1 (V2 surface; "
                "pyslha/xslha are core deps). "
                "Install it with `pip install -U Jarvis-HEP-Portal` "
                "(or `pip install -e ../Jarvis-Portal` for local development)."
            ) from inner
        _PORTAL_IMPORT_ERROR = exc
        return IOContext, MissingAdapterError, create_entry_point_registry
    return IOContext, MissingAdapterError, create_entry_point_registry


def get_io_registry():
    """Return the process-local **V2** Portal registry (builtins + entry points)."""
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


def _format_expression_input_pairs(values: Mapping[str, Any]) -> str:
    """Render ``[x : 1.0, y : 2.0]`` like V1 sample-log Dump lines."""
    parts = [f"{key} : {value}" for key, value in values.items()]
    return ", ".join(parts)


def evaluate_io_expression(
    expression: str,
    values: Mapping[str, Any],
    *,
    context: ExpressionContext | None = None,
    logger: Any = None,
) -> Any:
    """Evaluate a cached Dump expression with the shared HEP expression semantics."""
    text = str(expression).strip()
    if not text:
        raise ValueError("IO expression must not be empty.")
    global _IO_EXPRESSIONS, _IO_OPERAS_LOADED
    if context is None:
        if _IO_EXPRESSIONS is None:
            _IO_EXPRESSIONS = ExpressionContext()
        if not _IO_OPERAS_LOADED:
            from jarvishep2.operas_functions import (
                build_operas_expression_context,
                expression_uses_operas_function,
            )

            if expression_uses_operas_function(text):
                _IO_EXPRESSIONS = build_operas_expression_context(required=True)
                _IO_OPERAS_LOADED = True
        context = _IO_EXPRESSIONS
    compiled = context.compile(text)
    numeric_values: dict[str, float] = {}
    missing: list[str] = []
    for name in compiled.variable_names:
        numeric = _coerce_numeric_param(values.get(name))
        if numeric is None:
            missing.append(name)
            continue
        numeric_values[name] = numeric
    if missing:
        raise ValueError(f"Dump expression '{text}' misses parameters: {sorted(missing)}")
    try:
        value = float(compiled.evaluate(numeric_values))
    except MissingExpressionVariablesError as exc:
        raise ValueError(
            f"Dump expression '{text}' misses parameters: {list(exc.missing)}"
        ) from exc
    if logger is not None:
        # V1 sample-log shape for Portal Dump evaluations.
        logger.info(
            "Evaluating: expression \n\t-> {} \n    with input \t -> [{}] "
            "\n    Output \t\t-> {}".format(
                text,
                _format_expression_input_pairs(numeric_values),
                value,
            )
        )
    return value


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
    if evaluate_expression is not None:
        evaluator = evaluate_expression
    else:
        # Close over sample logger so Portal Dump lines land in Sample_running.log.
        def evaluator(expression: str, values: Mapping[str, Any], _logger=logger) -> Any:
            return evaluate_io_expression(expression, values, logger=_logger)

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
