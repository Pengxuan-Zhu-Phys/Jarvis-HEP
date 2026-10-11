"""Compile task expressions at the validation gate without evaluating samples."""

from collections.abc import Mapping
from typing import Any

from jarvishep2.expression import ExpressionContext
from jarvishep2.operas_functions import (
    build_operas_expression_context,
    expression_uses_operas_reference,
)


def _expressions(value: Any, path: str = ""):
    if isinstance(value, Mapping):
        for key, item in value.items():
            child = f"{path}.{key}" if path else str(key)
            # Mapper has its own closed-namespace compiler and diagnostics.
            if key in {"Mapper", "mapper", "kwargs", "sampler", "run_nested", "options"}:
                continue
            if key in {"expression", "selection", "target_expression"} and isinstance(item, str):
                yield child, item
            elif isinstance(item, (Mapping, list)):
                yield from _expressions(item, child)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _expressions(item, f"{path}[{index}]")


def _symbols(value: Any) -> set[str]:
    names: set[str] = set()
    if isinstance(value, Mapping):
        name = value.get("name")
        if isinstance(name, str) and name.isidentifier():
            names.add(name)
        for item in value.values():
            names.update(_symbols(item))
    elif isinstance(value, list):
        for item in value:
            names.update(_symbols(item))
    return names


def validate_expressions(config: Mapping[str, Any]):
    from jarvishep2.task_validation import issue

    # Limit the walk to execution blocks; arbitrary plug-in/CLI options are opaque.
    sampling = config.get("Sampling")
    blocks = {key: config[key] for key in ("Operas", "Calculators") if key in config}
    if isinstance(sampling, Mapping):
        blocks["Sampling"] = {key: sampling[key] for key in ("LogLikelihood", "selection") if key in sampling}
        bounds = sampling.get("Bounds")
        if isinstance(bounds, Mapping) and "target_expression" in bounds:
            blocks["Sampling"]["Bounds"] = {"target_expression": bounds["target_expression"]}
    expressions = list(_expressions(blocks))
    symbols = _symbols(config)
    internal = ExpressionContext()
    external = None
    issues = []
    for path, expression in expressions:
        try:
            context = internal
            if expression_uses_operas_reference(expression):
                if external is None:
                    external = build_operas_expression_context(required=True)
                context = external
            context.compile(expression, symbols=symbols)
        except Exception as exc:
            issues.append(issue(
                "error", "JV2-EXP-001", path,
                f"Cannot compile expression {expression!r}: {type(exc).__name__}: {exc}",
            ))
    return issues
