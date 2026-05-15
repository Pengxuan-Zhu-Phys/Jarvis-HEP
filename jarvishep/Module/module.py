#!/usr/bin/env python3
from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify

from jarvishep.base import Base
from jarvishep.inner_func import build_expression_context, update_const, update_funcs

class Module(Base): 
    def __init__(self, name, inputs=None, outputs=None, selection: str | None = None):
        super().__init__()
        self.logger = None
        self.name = name
        self.input = inputs if inputs else []
        self.output = outputs if outputs else []
        self.inputs = {}
        self.outputs = {}
        self.required_modules = []
        self._funcs = {}
        self.selection = None
        self._selection_deps = ()
        self._selection_fn = None
        if selection is not None:
            self.set_selection(selection)

    def execute(self):
        raise NotImplementedError("This method should be implemented by subclasses.")

    def set_funcs(self, funcs):
        self._funcs = dict(funcs or {})
        if self.selection:
            self._compile_selection()

    def set_selection(self, expression: str | None):
        if expression is None or str(expression).strip() == "":
            self.selection = None
            self._selection_deps = ()
            self._selection_fn = None
            return
        if not isinstance(expression, str):
            raise ValueError(
                f"Invalid module selection for '{self.name}': expected string or null, got {expression!r}"
            )
        self.selection = str(expression)
        self._compile_selection()

    def selection_checker(self, available_keys: Iterable[str]) -> tuple[bool, set[str]]:
        if not self.selection:
            return True, set()
        avail = {str(key) for key in available_keys}
        missing = set(self._selection_deps) - avail
        return len(missing) == 0, missing

    def evaluate_selection(self, values: Mapping[str, Any]) -> bool:
        if not self.selection:
            return True
        if self._selection_fn is None:
            self._compile_selection()
        ok, missing = self.selection_checker(values.keys())
        if not ok:
            raise KeyError(
                f"Module selection for '{self.name}' is missing required observables: {sorted(missing)}"
            )
        args = []
        for dep in self._selection_deps:
            args.append(self._ensure_scalar(values[dep], dep, "input"))
        result = self._selection_fn(*args)
        return bool(self._ensure_scalar(result, self.selection, "output"))

    def _compile_selection(self) -> None:
        funcs = update_funcs(dict(self._funcs or {}))
        parse_locals, numeric_modules = build_expression_context(
            funcs=funcs,
            consts=update_const({}),
        )
        sym_expr = sp.sympify(self.selection, locals=parse_locals)
        self._selection_deps = tuple(sorted(str(symbol) for symbol in sym_expr.free_symbols))
        sym_locals = {key: sp.Symbol(key) for key in self._selection_deps}
        sym_locals.update(parse_locals)
        expr = sp.sympify(self.selection, locals=sym_locals)
        symbols = [sym_locals[key] for key in self._selection_deps]
        self._selection_fn = lambdify(symbols, expr, modules=[numeric_modules, "numpy", math])

    @staticmethod
    def _ensure_scalar(value: Any, field: str, kind: str) -> Any:
        arr = np.asarray(value)
        if arr.ndim != 0:
            raise ValueError(
                f"Module selection expects scalar {kind} for '{field}', got {type(value).__name__} with shape {arr.shape}."
            )
        return arr.item()
