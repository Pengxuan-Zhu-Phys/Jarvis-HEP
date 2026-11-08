#!/usr/bin/env python3
"""Log-likelihood evaluation inside Workers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import inf, isnan
from numbers import Real
from typing import Any

import numpy as np

from jarvishep2.expression import (
    CompiledExpression,
    ExpressionContext,
    MissingExpressionVariablesError,
)


def _sanitize_eval_values(values: Mapping[str, Any]) -> dict[str, Any]:
    """Map JSON/portal nulls to NaN so Gaussian kernels do not TypeError."""
    out: dict[str, Any] = dict(values)
    for key, value in list(out.items()):
        if value is None:
            out[key] = float("nan")
    return out


class LogLikelihoodEvaluator:
    """Compile and evaluate configured LogLikelihood expressions."""

    def __init__(
        self,
        expressions: Sequence[Mapping[str, Any]] | None,
        *,
        expression_context: ExpressionContext | None = None,
        strict: bool = False,
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
        self._strict = strict
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

    def _evaluate_terms(
        self, payload: Mapping[str, Any], values: dict[str, float], logger: Any = None,
    ) -> None:
        if not self._compiled:
            if self._strict:
                raise ValueError("sampling requires a LogLikelihood expression")
            values["LogL"] = -inf
            return

        inputs = _sanitize_eval_values(payload)
        explicit_total = any(name == "LogL" for name, _ in self._compiled)
        total = 0.0
        for name, compiled in self._compiled:
            try:
                result = compiled.evaluate(inputs)
            except MissingExpressionVariablesError as exc:
                raise KeyError(
                    f"LogLikelihood expression '{name}' misses observables: {list(exc.missing)}"
                ) from exc
            if self._strict:
                likelihood = require_log_density(result)
            else:
                try:
                    likelihood = _finite_or_neginf(float(result))
                except (TypeError, ValueError, OverflowError):
                    likelihood = -inf

            if logger is not None:
                input_text = ", ".join(
                    f"{key} : {inputs[key]}" for key in compiled.variable_names if key in inputs
                )
                logger.info(
                    f"Evaluating   {name}: \n"
                    f"   expression \t-> {compiled.expression}\n"
                    f"   with input \t-> [{input_text}] \n"
                    f"   Output \t\t-> {likelihood}"
                )
            values[name] = inputs[name] = likelihood
            if name == "LogL":
                total = likelihood
            elif not explicit_total:
                total += likelihood
        values["LogL"] = require_log_density(total) if self._strict else _finite_or_neginf(total)

    def evaluate(self, observables: Mapping[str, Any]) -> dict[str, float]:
        """Evaluate the ordered terms and return the final LogL."""
        values: dict[str, float] = {}
        self._evaluate_terms(observables, values)
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
        """Worker-facing evaluation that writes into sample_info.

        On failure, preserve evaluated terms for diagnostics and re-raise so
        the Worker records a failed sample.
        """
        observables = sample_info.get("observables", {})
        if not isinstance(observables, dict):
            raise TypeError("sample_info['observables'] must be a dict")
        values: dict[str, float] = {}
        try:
            self._evaluate_terms(observables, values, self._sample_logger(sample_info))
        except Exception:
            values["LogL"] = -inf
            raise
        finally:
            if "LogL" in values:
                observables.update(values)
                sample_info["observables"] = observables
                sample_info["likelihood"] = values["LogL"]
        return values["LogL"]


def require_log_density(value: Any) -> float:
    """Return a real score, preserving negative infinity as zero density."""
    if isinstance(value, Real) and not isinstance(value, (bool, np.bool_)):
        try:
            score = float(value)
        except OverflowError:
            pass
        else:
            if not isnan(score) and score != inf and (score != -inf or value == -inf):
                return score
    raise ValueError(f"invalid log density: {value!r}; expected a real number or -inf")


def _finite_or_neginf(value: float) -> float:
    return -inf if isnan(value) or value == inf else value


__all__ = ["LogLikelihoodEvaluator", "require_log_density"]
