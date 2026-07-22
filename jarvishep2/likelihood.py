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
        """Return likelihood terms computed from observables.

        Unusable totals are written as ``-inf`` (no expressions, non-finite
        results, null inputs). Callers that need fail-loud on missing symbols
        still raise.
        """
        values: dict[str, float] = {}
        payload = dict(observables)
        eval_values = _sanitize_eval_values(payload)
        if not self._compiled:
            values["LogL"] = float("-inf")
            return values
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
            if result is None:
                likelihood = float("-inf")
            else:
                try:
                    likelihood = _finite_or_neginf(float(result))
                except (TypeError, ValueError):
                    likelihood = float("-inf")
            values[name] = likelihood
            eval_values[name] = likelihood
            if name == "LogL":
                total_loglikelihood = likelihood
            elif not explicit_total:
                total_loglikelihood += likelihood
        values["LogL"] = _finite_or_neginf(float(total_loglikelihood))
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

        On missing symbols / evaluation failure, writes ``LogL = -inf`` and
        re-raises so the Worker can mark the sample Failed while feedback still
        has a defined logL (D13.8).
        """
        observables = sample_info.get("observables", {})
        if not isinstance(observables, dict):
            raise TypeError("sample_info['observables'] must be a dict")
        slogger = self._sample_logger(sample_info)
        # Evaluate term-by-term so each expression can be sample-logged like V1.
        payload = dict(observables)
        eval_values = _sanitize_eval_values(payload)
        if not self._compiled:
            observables["LogL"] = float("-inf")
            sample_info["observables"] = observables
            sample_info["likelihood"] = float("-inf")
            return float("-inf")
        explicit_total = any(name == "LogL" for name, _ in self._compiled)
        total_loglikelihood = 0.0
        values: dict[str, float] = {}
        try:
            for name, compiled in self._compiled:
                try:
                    result = compiled.evaluate(eval_values)
                except MissingExpressionVariablesError as exc:
                    raise KeyError(
                        f"LogLikelihood expression '{name}' misses observables: "
                        f"{list(exc.missing)}"
                    ) from exc
                if isinstance(result, np.generic):
                    result = result.item()
                if result is None:
                    likelihood = float("-inf")
                else:
                    try:
                        likelihood = _finite_or_neginf(float(result))
                    except (TypeError, ValueError):
                        # Null/non-numeric residual from a calculator portal → unusable.
                        likelihood = float("-inf")
                values[name] = likelihood
                eval_values[name] = likelihood
                if name == "LogL":
                    total_loglikelihood = likelihood
                elif not explicit_total:
                    total_loglikelihood += likelihood
                if slogger is not None:
                    used = {
                        key: payload.get(key, eval_values.get(key))
                        for key in compiled.variable_names
                        if key in eval_values or key in payload
                    }
                    input_text = ", ".join(
                        f"{key} : {val}" for key, val in used.items()
                    )
                    slogger.info(
                        f"Evaluating   {name}: \n"
                        f"   expression \t-> {compiled.expression}\n"
                        f"   with input \t-> [{input_text}] \n"
                        f"   Output \t\t-> {likelihood}"
                    )
            values["LogL"] = _finite_or_neginf(float(total_loglikelihood))
            observables.update(values)
            sample_info["observables"] = observables
            likelihood = float(values["LogL"])
            sample_info["likelihood"] = likelihood
            return likelihood
        except Exception:
            # Guarantee a defined total for feedback/archive before re-raise.
            # Missing symbols / hard eval bugs still fail the sample; null physics
            # values are sanitized above and should not reach this path.
            observables["LogL"] = float("-inf")
            observables.update({k: v for k, v in values.items()})
            sample_info["observables"] = observables
            sample_info["likelihood"] = float("-inf")
            raise


def _finite_or_neginf(value: float) -> float:
    """Map non-finite floats to ``-inf`` (D13.8 likelihood contract)."""
    x = float(value)
    if x != x or x == float("inf") or x == float("-inf"):  # nan or ±inf
        # Preserve -inf; collapse +inf/nan to -inf (unusable likelihood).
        if x == float("-inf"):
            return float("-inf")
        return float("-inf")
    return x


__all__ = ["LogLikelihoodEvaluator"]
