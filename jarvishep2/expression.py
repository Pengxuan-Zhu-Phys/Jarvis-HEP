#!/usr/bin/env python3
"""Shared compile-once expression runtime for Jarvis-HEP V2."""

from __future__ import annotations

import threading
from collections import namedtuple
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import sympy as sp
from sympy.utilities.lambdify import lambdify

from jarvishep2.inner_func import EXPRESSION_CONSTANTS, NUMERIC_MODULES


ExpressionCacheInfo = namedtuple(
    "ExpressionCacheInfo",
    ("hits", "misses", "maxsize", "currsize"),
)


class MissingExpressionVariablesError(KeyError):
    """Raised when a compiled expression cannot bind all of its variables."""

    def __init__(self, expression: str, missing: Iterable[str]) -> None:
        self.expression = str(expression)
        self.missing = tuple(sorted(str(name) for name in missing))
        super().__init__(
            f"Expression '{self.expression}' misses variables: {list(self.missing)}"
        )


@dataclass(frozen=True)
class CompiledExpression:
    """An immutable SymPy expression and its reusable numerical callable."""

    expression: str
    variable_names: tuple[str, ...]
    _callable: Callable[..., Any] = field(repr=False, compare=False)

    def evaluate(self, values: Mapping[str, Any]) -> Any:
        """Bind values by name and evaluate without reparsing or recompiling."""
        missing = [name for name in self.variable_names if name not in values]
        if missing:
            raise MissingExpressionVariablesError(self.expression, missing)
        return self._callable(*(values[name] for name in self.variable_names))

    __call__ = evaluate


class ExpressionContext:
    """Own an expression language namespace and a process-local compile cache.

    Contexts are intentionally process-local: Workers construct their runtime after
    ``spawn``, while samplers and Portal I/O keep a control-process context. Cache keys
    include the expression text and explicitly declared symbol names.
    """

    def __init__(
        self,
        *,
        functions: Mapping[str, Any] | None = None,
        constants: Mapping[str, Any] | None = None,
        parse_locals: Mapping[str, Any] | None = None,
        maxsize: int | None = 128,
    ) -> None:
        self._functions = dict(NUMERIC_MODULES)
        if functions:
            self._functions.update({str(name): value for name, value in functions.items()})
        self._constants: dict[str, Any] = dict(EXPRESSION_CONSTANTS)
        if constants:
            self._constants.update({str(name): value for name, value in constants.items()})
        self._parse_locals = {
            str(name): value for name, value in (parse_locals or {}).items()
        }
        self._maxsize = None if maxsize is None else max(1, int(maxsize))
        self._cache: dict[tuple[str, tuple[str, ...]], CompiledExpression] = {}
        self._cache_order: list[tuple[str, tuple[str, ...]]] = []
        self._hits = 0
        self._misses = 0
        self._lock = threading.RLock()

    def compile(
        self,
        expression: str,
        *,
        symbols: Iterable[str] = (),
    ) -> CompiledExpression:
        """Compile once for this context and return the cached object thereafter."""
        text = str(expression).strip()
        if not text:
            raise ValueError("Expression must not be empty.")
        symbol_names = tuple(sorted({str(name) for name in symbols}))
        key = (text, symbol_names)
        with self._lock:
            cached = self._cache.get(key)
            if cached is not None:
                self._hits += 1
                return cached

            parse_locals = {name: sp.Symbol(name) for name in symbol_names}
            # Namespace objects supplied by Jarvis-Operas must win over a
            # same-named observable so ``namespace.function(x)`` remains
            # callable. Constants remain reserved above both.
            parse_locals.update(self._parse_locals)
            # V1 constants are reserved language names and win over an
            # accidentally identical observable name.
            parse_locals.update(self._constants)
            parsed = sp.sympify(text, locals=parse_locals)
            variable_names = tuple(sorted(str(symbol) for symbol in parsed.free_symbols))
            numerical = lambdify(
                list(variable_names),
                parsed,
                modules=[dict(self._functions), "numpy"],
            )
            compiled = CompiledExpression(text, variable_names, numerical)
            self._misses += 1
            self._cache[key] = compiled
            self._cache_order.append(key)
            if self._maxsize is not None and len(self._cache_order) > self._maxsize:
                oldest = self._cache_order.pop(0)
                self._cache.pop(oldest, None)
            return compiled

    def evaluate(
        self,
        expression: str,
        values: Mapping[str, Any],
        *,
        symbols: Iterable[str] | None = None,
    ) -> Any:
        """Compile/cache an expression and evaluate it against a value mapping."""
        normalized_values = {str(name): value for name, value in values.items()}
        declared = normalized_values.keys() if symbols is None else symbols
        return self.compile(expression, symbols=declared).evaluate(normalized_values)

    def update_functions(self, functions: Mapping[str, Any]) -> None:
        """Extend the numerical namespace and invalidate previously compiled callables."""
        with self._lock:
            self._functions.update({str(name): value for name, value in functions.items()})
            self.clear_cache()

    def update_namespace(
        self,
        *,
        functions: Mapping[str, Any] | None = None,
        parse_locals: Mapping[str, Any] | None = None,
    ) -> None:
        """Extend parsing and numerical namespaces as one atomic cache update."""
        with self._lock:
            if functions:
                self._functions.update(
                    {str(name): value for name, value in functions.items()}
                )
            if parse_locals:
                self._parse_locals.update(
                    {str(name): value for name, value in parse_locals.items()}
                )
            self.clear_cache()

    def clear_cache(self) -> None:
        with self._lock:
            self._cache.clear()
            self._cache_order.clear()
            self._hits = 0
            self._misses = 0

    def cache_info(self) -> ExpressionCacheInfo:
        with self._lock:
            return ExpressionCacheInfo(
                self._hits,
                self._misses,
                self._maxsize,
                len(self._cache),
            )


__all__ = [
    "CompiledExpression",
    "ExpressionCacheInfo",
    "ExpressionContext",
    "MissingExpressionVariablesError",
]
