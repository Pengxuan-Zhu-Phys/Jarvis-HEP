#!/usr/bin/env python3
"""Nuisance expression registry and pass-condition (D13.4).

Ports V1 ``Module/nuisance_LogLikelihood.py`` / ``nuisance_passCondition.py``
onto the shared V2 :class:`ExpressionContext` (compile-once, no jarvishep).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from jarvishep2.expression import CompiledExpression, ExpressionContext, MissingExpressionVariablesError


@dataclass(frozen=True)
class NuisanceTerm:
    """One compiled nuisance LogL or pass-condition expression."""

    name: str
    expression: str
    compiled: CompiledExpression

    @property
    def deps(self) -> tuple[str, ...]:
        return self.compiled.variable_names

    def can_eval(self, available_keys: Iterable[str]) -> tuple[bool, set[str]]:
        avail = {str(k) for k in available_keys}
        missing = set(self.deps) - avail
        return (not missing), missing

    def eval(self, values: Mapping[str, Any]) -> Any:
        return self.compiled.evaluate(values)

    def eval_bool(self, values: Mapping[str, Any]) -> bool:
        return bool(self.eval(values))

    def eval_float(self, values: Mapping[str, Any]) -> float:
        value = self.eval(values)
        return float(value)


class NuisanceExpressionRegistry:
    """Compile-once registry of named nuisance LogL terms."""

    def __init__(self, context: ExpressionContext | None = None) -> None:
        self._context = context or ExpressionContext()
        self._terms: dict[str, NuisanceTerm] = {}

    @property
    def names(self) -> list[str]:
        return list(self._terms.keys())

    def __contains__(self, name: object) -> bool:
        return str(name) in self._terms

    def get(self, name: str) -> NuisanceTerm:
        return self._terms[str(name)]

    def set_config(self, name: str, expression: str) -> NuisanceTerm:
        compiled = self._context.compile(str(expression))
        term = NuisanceTerm(name=str(name), expression=str(expression), compiled=compiled)
        self._terms[term.name] = term
        return term

    def load_from_config(self, items: Sequence[Mapping[str, Any]] | None) -> None:
        for item in items or []:
            if not isinstance(item, Mapping):
                continue
            name = item.get("name")
            expr = item.get("expression")
            if name is None or expr is None:
                continue
            self.set_config(str(name), str(expr))

    def evaluate_all(self, values: Mapping[str, Any]) -> dict[str, float]:
        out: dict[str, float] = {}
        for name, term in self._terms.items():
            try:
                out[name] = term.eval_float(values)
            except MissingExpressionVariablesError:
                raise
        return out

    def total(self, values: Mapping[str, Any]) -> float:
        terms = self.evaluate_all(values)
        return float(sum(terms.values())) if terms else 0.0


class NuisancePassConditionRegistry:
    """Compile-once registry of named pass-condition predicates."""

    def __init__(self, context: ExpressionContext | None = None) -> None:
        self._context = context or ExpressionContext()
        self._terms: dict[str, NuisanceTerm] = {}

    @property
    def names(self) -> list[str]:
        return list(self._terms.keys())

    def set_config(self, name: str, expression: str) -> NuisanceTerm:
        compiled = self._context.compile(str(expression))
        term = NuisanceTerm(name=str(name), expression=str(expression), compiled=compiled)
        self._terms[term.name] = term
        return term

    def load_from_config(self, items: Sequence[Mapping[str, Any]] | None) -> None:
        for item in items or []:
            if not isinstance(item, Mapping):
                continue
            name = item.get("name")
            expr = item.get("expression")
            if name is None or expr is None:
                continue
            self.set_config(str(name), str(expr))

    def evaluate_all(self, values: Mapping[str, Any]) -> dict[str, bool]:
        return {name: term.eval_bool(values) for name, term in self._terms.items()}

    def all_pass(self, values: Mapping[str, Any]) -> bool:
        if not self._terms:
            return True
        results = self.evaluate_all(values)
        return all(results.values())


def extract_nuisance_config(config: Mapping[str, Any] | None) -> dict[str, Any] | None:
    """Return the Nuisance block from Sampling.Nuisance or top-level Nuisance."""
    if not isinstance(config, Mapping):
        return None
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    block = None
    if isinstance(sampling, Mapping):
        block = sampling.get("Nuisance") or sampling.get("nuisance")
    if block is None:
        block = config.get("Nuisance") or config.get("nuisance")
    if not isinstance(block, Mapping):
        return None
    return dict(block)


def parse_nuisance_variable(block: Mapping[str, Any]) -> tuple[str, float, float]:
    """Return (name, zmin, zmax) for the first nuisance variable (Profile1D)."""
    vars_list = block.get("Variables") or block.get("variables") or []
    if not isinstance(vars_list, list) or not vars_list:
        raise ValueError("Nuisance.Variables must contain at least one variable")
    var = vars_list[0]
    if not isinstance(var, Mapping):
        raise ValueError("Nuisance.Variables[0] must be a mapping")
    name = str(var.get("name") or "nuisance").strip()
    dist = var.get("distribution") if isinstance(var.get("distribution"), Mapping) else {}
    params = dist.get("parameters") if isinstance(dist.get("parameters"), Mapping) else {}
    zmin = float(params.get("min", 0.0))
    zmax = float(params.get("max", 1.0))
    if zmin > zmax:
        zmin, zmax = zmax, zmin
    if zmin == zmax:
        zmax = zmin + 1.0
    return name, zmin, zmax


__all__ = [
    "NuisanceExpressionRegistry",
    "NuisancePassConditionRegistry",
    "NuisanceTerm",
    "extract_nuisance_config",
    "parse_nuisance_variable",
]
