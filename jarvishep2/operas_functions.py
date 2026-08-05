#!/usr/bin/env python3
"""Startup discovery bridge from Jarvis-Operas into HEP expressions.

External functions keep the V1 YAML surface: they are registered in
Jarvis-Operas and referenced by their full ``namespace.name`` form.
Constants (D23) use the bare form ``pdg.mZ``; callables keep
``namespace.function(...)``.  No HEP-specific function-registration block
is added to task YAML.
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
    # D23: diagnostic/log only — not consumed by ExpressionContext.
    constants: Mapping[str, float] = field(default_factory=dict)

    def create_context(self, *, maxsize: int | None = 128) -> ExpressionContext:
        return ExpressionContext(
            functions=self.numeric_functions,
            parse_locals=self.parse_locals,
            maxsize=maxsize,
        )


_DISCOVERY_LOCK = threading.RLock()
_DISCOVERED_REGISTRIES: set[int] = set()

# Qualified Operas reference: namespace.short (function call *or* bare constant).
# Guards reject path-like dotted fragments that would force a ~1 s cold discovery
# (design §8.1): ``/usr/local/bin/foo.sh``, ``./run.sh``, ``results/out.dat``,
# ``C:\\tools\\a.exe``.  Residual bare ``run.sh`` is unreachable for config walks
# because only expression/selection/target_expression keys are scanned.
_QUALIFIED_REFERENCE = re.compile(
    r"(?<![/\\.\w])[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)+(?![\w.]*[/\\])"
)
# Backward-compatible name (pre-D23 required a trailing ``(``).
_QUALIFIED_FUNCTION_CALL = _QUALIFIED_REFERENCE

# Only these mapping keys hold YAML expression text. Command/path strings such as
# calculator ``cmd`` / ``installation`` must never trigger Operas discovery.
# ``expression`` also covers Sampling.Mapper list items ({name, expression}).
_EXPRESSION_TEXT_KEYS = frozenset(
    {
        "expression",
        "selection",
        "target_expression",
    }
)
# Container keys whose *values* are expression-bearing structures even when the
# container key itself is not an expression field (D23.9: Sampling.Mapper list).
_EXPRESSION_CONTAINER_KEYS = frozenset(
    {
        "Mapper",
        "mapper",
    }
)


def _string_uses_operas_reference(text: str) -> bool:
    return bool(_QUALIFIED_REFERENCE.search(text))


# Pre-D23 alias.
_string_uses_operas_function = _string_uses_operas_reference


def _mapper_block_uses_operas_reference(value: Any) -> bool:
    """Scan Sampling.Mapper list or legacy name→expr map for qualified refs."""
    if isinstance(value, str):
        return _string_uses_operas_reference(value)
    if isinstance(value, Mapping):
        # Canonical: {name, expression}; also tolerate legacy free-form map values.
        if "expression" in value and isinstance(value.get("expression"), str):
            if _string_uses_operas_reference(str(value.get("expression") or "")):
                return True
        for item in value.values():
            if _mapper_block_uses_operas_reference(item):
                return True
        return False
    if isinstance(value, (list, tuple)):
        return any(_mapper_block_uses_operas_reference(item) for item in value)
    return False


def expression_uses_operas_reference(
    value: Any,
    *,
    _as_expression_text: bool = True,
) -> bool:
    """Return whether *expression text* contains a qualified Operas reference.

    Matches both call form (``math.add(x, 2)``) and bare constants
    (``pdg.mZ``).  A bare string is treated as expression text only when the
    caller already holds a single formula (``_as_expression_text=True``, the
    default) — e.g. ``evaluate_selection`` / Dump. Nested config walks set the
    flag to ``False`` except under expression-bearing keys
    (``expression`` / ``selection`` / ``target_expression``), so calculator
    ``cmd`` / ``path`` / install strings never force Operas discovery.
    """
    if isinstance(value, str):
        return bool(_as_expression_text and _string_uses_operas_reference(value))
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_name = str(key)
            if key_name in _EXPRESSION_TEXT_KEYS:
                if expression_uses_operas_reference(item, _as_expression_text=True):
                    return True
            elif key_name in _EXPRESSION_CONTAINER_KEYS:
                # Mapper list / legacy name→expr map (D23.9).
                if _mapper_block_uses_operas_reference(item):
                    return True
            elif isinstance(item, (Mapping, list, tuple)):
                if expression_uses_operas_reference(item, _as_expression_text=False):
                    return True
        return False
    if isinstance(value, (list, tuple)):
        return any(
            expression_uses_operas_reference(
                item, _as_expression_text=_as_expression_text
            )
            for item in value
        )
    return False


# Pre-D23 public name kept as alias (exported in ``__all__``).
expression_uses_operas_function = expression_uses_operas_reference


def operas_expression_references_required(expressions: Any) -> bool:
    """Avoid bootstrapping Operas for payloads that only use the internal core."""
    return expression_uses_operas_reference(expressions)


# Pre-D23 public name kept as alias.
operas_expression_functions_required = operas_expression_references_required


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
    """Discover registered functions/constants and freeze numerical implementations."""
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
                "A qualified Jarvis-Operas expression requires Jarvis-Operas, "
                "which is a core dependency of jarvishep2."
            ) from exc
        return OperasExpressionSnapshot()

    # D23: optional on older Operas installs that lack the constant table API.
    try:
        from jarvis_operas import build_constant_dicts
    except ImportError:  # pragma: no cover - pre-D23 Operas
        build_constant_dicts = None  # type: ignore[assignment]

    target_registry = registry or get_global_operas_registry()
    loaded: tuple[str, ...] = ()
    if discover_plugins:
        loaded = _discover_entrypoints_once(target_registry, discover_entrypoints)
    register_dict = build_register_dicts(target_registry)
    if build_constant_dicts is not None:
        constant_dict = build_constant_dicts(target_registry)
        parse_locals, numeric_functions = build_sympy_dicts(
            register_dict,
            constants=constant_dict,
            include_all=True,
        )
    else:  # pragma: no cover - pre-D23 Operas
        constant_dict = {}
        parse_locals, numeric_functions = build_sympy_dicts(
            register_dict, include_all=True
        )
    return OperasExpressionSnapshot(
        parse_locals=dict(parse_locals),
        numeric_functions=dict(numeric_functions),
        registered_names=tuple(sorted(register_dict)),
        entrypoints_loaded=loaded,
        constants=dict(constant_dict),
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
            "Jarvis-Operas is required for a registered Operas operator "
            "(core dependency of jarvishep2)."
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
    "expression_uses_operas_reference",
    "expression_uses_operas_function",  # alias
    "operas_expression_references_required",
    "operas_expression_functions_required",  # alias
]
