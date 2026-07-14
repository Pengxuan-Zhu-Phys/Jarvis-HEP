#!/usr/bin/env python3
"""Log-likelihood evaluation inside Workers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from jarvishep2.expression import (
    CompiledExpression,
    ExpressionContext,
    MissingExpressionVariablesError,
)


class LogLikelihoodEvaluator:
    """Compile and evaluate configured LogLikelihood expressions."""

    def __init__(
        self,
        expressions: Sequence[Mapping[str, Any]] | None,
        *,
        expression_context: ExpressionContext | None = None,
    ) -> None:
        if expression_context is None:
            from jarvishep2.operas_functions import (
                build_operas_expression_context,
                expression_uses_operas_function,
            )

            expression_context = (
                build_operas_expression_context(required=True)
                if expression_uses_operas_function(expressions)
                else ExpressionContext()
            )
        self._context = expression_context
        self._compiled: list[tuple[str, CompiledExpression]] = []
        known_symbols = ("x", "y", "z", "shift", "calc_z", "LogL", "LogL_Z")
        for item in expressions or []:
            if not isinstance(item, Mapping):
                continue
            name = str(item.get("name", "LogL"))
            expression = str(item.get("expression", "")).strip()
            if not expression:
                continue
            compiled = self._context.compile(expression, symbols=known_symbols)
            self._compiled.append((name, compiled))

    def evaluate(self, observables: Mapping[str, Any]) -> dict[str, float]:
        """Return likelihood terms computed from observables."""
        values: dict[str, float] = {}
        payload = dict(observables)
        eval_values = dict(payload)
        explicit_total = any(name == "LogL" for name, _ in self._compiled)
        total_loglikelihood = 0.0
        for name, compiled in self._compiled:
            try:
                result = compiled.evaluate(eval_values)
            except MissingExpressionVariablesError as exc:
                raise KeyError(
                    f"LogLikelihood expression '{name}' misses observables: {list(exc.missing)}"
                ) from exc
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

    def _sample_logger(self, sample_info: Mapping[str, Any]) -> Any | None:
        parent = sample_info.get("logger") if isinstance(sample_info, Mapping) else None
        if parent is None:
            return None
        base_name = ""
        if isinstance(sample_info, Mapping):
            base_name = str(sample_info.get("logger_name") or "").strip()
        if not base_name:
            uuid = str(sample_info.get("uuid") or "UNKNOWN") if isinstance(sample_info, Mapping) else "UNKNOWN"
            base_name = f"Sample@{uuid}"
        binder = getattr(parent, "bind", None)
        if callable(binder):
            return binder(module=f"{base_name} (Likelihood)")
        return parent

    def calculate(self, sample_info: dict[str, Any]) -> float:
        """Worker-facing evaluation that writes into sample_info."""
        observables = sample_info.get("observables", {})
        if not isinstance(observables, dict):
            raise TypeError("sample_info['observables'] must be a dict")
        slogger = self._sample_logger(sample_info)
        # Evaluate term-by-term so each expression can be sample-logged like V1.
        payload = dict(observables)
        eval_values = dict(payload)
        explicit_total = any(name == "LogL" for name, _ in self._compiled)
        total_loglikelihood = 0.0
        values: dict[str, float] = {}
        for name, compiled in self._compiled:
            try:
                result = compiled.evaluate(eval_values)
            except MissingExpressionVariablesError as exc:
                raise KeyError(
                    f"LogLikelihood expression '{name}' misses observables: {list(exc.missing)}"
                ) from exc
            if isinstance(result, np.generic):
                result = result.item()
            likelihood = float(result)
            values[name] = likelihood
            eval_values[name] = likelihood
            if name == "LogL":
                total_loglikelihood = likelihood
            elif not explicit_total:
                total_loglikelihood += likelihood
            if slogger is not None:
                used = {
                    key: eval_values[key]
                    for key in compiled.variable_names
                    if key in eval_values
                }
                input_text = ", ".join(f"{key} : {val}" for key, val in used.items())
                slogger.info(
                    f"Evaluating   {name}: \n"
                    f"   expression \t-> {compiled.expression}\n"
                    f"   with input \t-> [{input_text}] \n"
                    f"   Output \t\t-> {likelihood}"
                )
        values["LogL"] = float(total_loglikelihood)
        observables.update(values)
        sample_info["observables"] = observables
        likelihood = values.get("LogL")
        if likelihood is None and values:
            likelihood = next(iter(values.values()))
        sample_info["likelihood"] = float(likelihood)
        return float(likelihood)


__all__ = ["LogLikelihoodEvaluator"]
