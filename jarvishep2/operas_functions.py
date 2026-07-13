#!/usr/bin/env python3
"""Startup discovery bridge from Jarvis-Operas into HEP expressions.

External functions keep the V1 YAML surface: they are registered in
Jarvis-Operas and referenced by their full ``namespace.function(...)`` name.
No HEP-specific function-registration block is added to task YAML.
"""

from __future__ import annotations

import re
import threading
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

from jarvishep2.expression import ExpressionContext


@dataclass(frozen=True)
class OperasExpressionSnapshot:
    """Startup snapshot used to construct an expression context."""

    parse_locals: Mapping[str, Any] = field(default_factory=dict)
    numeric_functions: Mapping[str, Callable[..., Any]] = field(default_factory=dict)
    registered_names: tuple[str, ...] = ()
    entrypoints_loaded: tuple[str, ...] = ()

    def create_context(self, *, maxsize: int | None = 128) -> ExpressionContext:
        return ExpressionContext(
            functions=self.numeric_functions,
            parse_locals=self.parse_locals,
            maxsize=maxsize,
        )


_DISCOVERY_LOCK = threading.RLock()
_DISCOVERED_REGISTRIES: set[int] = set()
_QUALIFIED_FUNCTION_CALL = re.compile(r"\b[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)+\s*\(")


def expression_uses_operas_function(value: Any) -> bool:
    """Return whether a config payload contains a qualified function call."""
    if isinstance(value, str):
        return bool(_QUALIFIED_FUNCTION_CALL.search(value))
    if isinstance(value, Mapping):
        return any(expression_uses_operas_function(item) for item in value.values())
    if isinstance(value, (list, tuple)):
        return any(expression_uses_operas_function(item) for item in value)
    return False


def operas_expression_functions_required(expressions: Any) -> bool:
    """Avoid bootstrapping Operas for expressions using only the internal core."""
    return expression_uses_operas_function(expressions)


def _discover_entrypoints_once(
    registry: Any,
    discover_entrypoints: Callable[..., Any],
) -> tuple[str, ...]:
    registry_key = id(registry)
    with _DISCOVERY_LOCK:
        if registry_key in _DISCOVERED_REGISTRIES:
            return ()
        loaded = tuple(sorted(str(name) for name in discover_entrypoints(registry)))
        _DISCOVERED_REGISTRIES.add(registry_key)
        return loaded


def discover_operas_expression_functions(
    *,
    registry: Any | None = None,
    discover_plugins: bool = True,
    required: bool = False,
) -> OperasExpressionSnapshot:
    """Discover registered functions and freeze their numerical implementations."""
    try:
        from jarvis_operas import (
            build_register_dicts,
            build_sympy_dicts,
            discover_entrypoints,
            get_global_operas_registry,
        )
    except ImportError as exc:
        if required:
            raise ImportError(
                "A qualified Jarvis-Operas expression requires Jarvis-Operas; "
                "install the Operas extra."
            ) from exc
        return OperasExpressionSnapshot()

    target_registry = registry or get_global_operas_registry()
    loaded: tuple[str, ...] = ()
    if discover_plugins:
        loaded = _discover_entrypoints_once(target_registry, discover_entrypoints)
    register_dict = build_register_dicts(target_registry)
    parse_locals, numeric_functions = build_sympy_dicts(register_dict, include_all=True)
    return OperasExpressionSnapshot(
        parse_locals=dict(parse_locals),
        numeric_functions=dict(numeric_functions),
        registered_names=tuple(sorted(register_dict)),
        entrypoints_loaded=loaded,
    )


def build_operas_expression_context(
    *,
    registry: Any | None = None,
    discover_plugins: bool = True,
    maxsize: int | None = 128,
    required: bool = False,
) -> ExpressionContext:
    return discover_operas_expression_functions(
        registry=registry,
        discover_plugins=discover_plugins,
        required=required,
    ).create_context(maxsize=maxsize)


def ensure_operas_registry_discovered(registry: Any | None = None) -> Any:
    """Return a registry after persisted and entry-point functions are visible."""
    try:
        from jarvis_operas import discover_entrypoints, get_global_operas_registry
    except ImportError as exc:
        raise ImportError(
            "Jarvis-Operas is required for a registered Operas operator; "
            "install the Operas extra."
        ) from exc
    target_registry = registry or get_global_operas_registry()
    _discover_entrypoints_once(target_registry, discover_entrypoints)
    return target_registry


def reset_operas_discovery_for_tests() -> None:
    with _DISCOVERY_LOCK:
        _DISCOVERED_REGISTRIES.clear()


__all__ = [
    "OperasExpressionSnapshot",
    "build_operas_expression_context",
    "discover_operas_expression_functions",
    "ensure_operas_registry_discovered",
    "expression_uses_operas_function",
    "operas_expression_functions_required",
]
