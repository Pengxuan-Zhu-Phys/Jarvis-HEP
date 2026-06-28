#!/usr/bin/env python3
"""Log-likelihood evaluation inside Workers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify


class LogLikelihoodEvaluator:
    """Compile and evaluate configured LogLikelihood expressions."""

    def __init__(self, expressions: Sequence[Mapping[str, Any]] | None) -> None:
        self._compiled: list[tuple[str, list[str], Any]] = []
        parse_locals = {
            name: sp.Symbol(name)
            for name in ("x", "y", "z", "shift", "LogL", "LogL_Z")
        }
        numeric_modules = {"sin": np.sin, "cos": np.cos, "exp": np.exp, "log": np.log}
        for item in expressions or []:
            if not isinstance(item, Mapping):
                continue
            name = str(item.get("name", "LogL"))
            expression = str(item.get("expression", "")).strip()
            if not expression:
                continue
            expr = sp.sympify(expression, locals=parse_locals)
            var_names = [str(sym) for sym in expr.free_symbols]
            num_expr = lambdify(var_names, expr, modules=[numeric_modules, "numpy"])
            self._compiled.append((name, var_names, num_expr))

    def evaluate(self, observables: Mapping[str, Any]) -> dict[str, float]:
        """Return likelihood terms computed from observables."""
        values: dict[str, float] = {}
        payload = dict(observables)
        eval_values = dict(payload)
        explicit_total = any(name == "LogL" for name, _, _ in self._compiled)
        total_loglikelihood = 0.0
        for name, var_names, num_expr in self._compiled:
            symbol_values = {
                key: eval_values[key] for key in var_names if key in eval_values
            }
            missing = [key for key in var_names if key not in symbol_values]
            if missing:
                raise KeyError(
                    f"LogLikelihood expression '{name}' misses observables: {missing}"
                )
            result = num_expr(**symbol_values)
            if isinstance(result, np.generic):
                result = result.item()
            likelihood = float(result)
            values[name] = likelihood
            eval_values[name] = likelihood
            if name == "LogL":
                total_loglikelihood = likelihood
            elif not explicit_total:
                total_loglikelihood += likelihood
        values["LogL"] = float(total_loglikelihood)
        return values

    def calculate(self, sample_info: dict[str, Any]) -> float:
        """Worker-facing evaluation that writes into sample_info."""
        observables = sample_info.get("observables", {})
        if not isinstance(observables, dict):
            raise TypeError("sample_info['observables'] must be a dict")
        values = self.evaluate(observables)
        observables.update(values)
        sample_info["observables"] = observables
        likelihood = values.get("LogL")
        if likelihood is None and values:
            likelihood = next(iter(values.values()))
        sample_info["likelihood"] = float(likelihood)
        return float(likelihood)


__all__ = ["LogLikelihoodEvaluator"]