#!/usr/bin/env python3
from __future__ import annotations

from inspect import signature
from typing import Any, Mapping

import sympy as sp
from sympy.utilities.lambdify import lambdify

from jarvishep.IOs.IOs import InputFile, OutputFile
from jarvishep.inner_func import build_expression_context, update_const


class PortalInputFile(InputFile):
    async def write(self, param_values):
        adapter = _portal_adapter(self, "input")
        context = _portal_context(
            owner=self,
            runtime_values=param_values,
            evaluate_expression=lambda expression, values: _evaluate_expression(
                self,
                expression,
                values,
            ),
        )
        return await adapter.write_input(context, self._portal_spec("actions"), param_values)

    def _portal_spec(self, payload_key):
        return _portal_spec(self, payload_key, self.variables)


class PortalOutputFile(OutputFile):
    async def read(self):
        adapter = _portal_adapter(self, "output")
        context = _portal_context(owner=self, runtime_values={})
        return await adapter.read_output(context, self._portal_spec("variables"))

    def _portal_spec(self, payload_key):
        return _portal_spec(self, payload_key, self.variables)


def _portal_spec(owner, payload_key: str, payload: Any) -> dict[str, Any]:
    spec = dict(getattr(owner, "spec", None) or {})
    spec.update(
        {
            "name": owner.name,
            "path": owner.path,
            "type": owner.file_type,
            "save": owner.save,
            payload_key: payload,
        }
    )
    for key in ("header", "columns", "comment"):
        if hasattr(owner, key):
            spec[key] = getattr(owner, key)
    return spec


def _portal_context(owner, runtime_values, evaluate_expression=None):
    try:
        from jarvis_portal import IOContext
    except ImportError as exc:
        raise ImportError(
            "Portal-backed calculator IO requires Jarvis-HEP-Portal. "
            "Install it with `pip install Jarvis-HEP-Portal`."
        ) from exc

    funcs = dict(owner.funcs or {})
    kwargs = {
        "logger": owner.logger,
        "io_manager": owner.io_manager,
        "sample_uuid": owner.sample_uuid,
        "pack_id": owner.PackID,
        "sample_save_dir": owner.sample_save_dir,
        "module": owner.module,
        "resolve_path": owner.decode_path,
        "evaluate_expression": evaluate_expression,
        "runtime_values": runtime_values,
    }
    params = signature(IOContext).parameters
    if "funcs" in params:
        kwargs["funcs"] = funcs
    accepted_kwargs = {key: value for key, value in kwargs.items() if key in params}
    context = IOContext(**accepted_kwargs)
    for key, value in kwargs.items():
        if key not in params:
            object.__setattr__(context, key, value)
    if "funcs" not in params:
        object.__setattr__(context, "funcs", funcs)
    return context


def _portal_adapter(owner, direction: str):
    adapter = getattr(owner, "portal_adapter", None)
    if adapter is not None:
        return adapter
    try:
        from jarvis_portal import get as portal_get
    except ImportError as exc:
        raise ImportError(
            "Portal-backed calculator IO requires Jarvis-HEP-Portal. "
            "Install it with `pip install Jarvis-HEP-Portal`."
        ) from exc
    return portal_get(owner.file_type, direction)


def _evaluate_expression(owner, expression: str, values: Mapping[str, Any]) -> Any:
    parse_locals, numeric_modules = build_expression_context(
        funcs=dict(owner.funcs or {}),
        consts=update_const({}),
    )
    expr = sp.sympify(expression, locals=parse_locals)
    symbol_names = [str(par) for par in expr.free_symbols]
    symbol_name_set = set(symbol_names)
    num_expr = lambdify(symbol_names, expr, modules=[numeric_modules, "numpy"])
    symbol_values = {
        str(key): value
        for key, value in values.items()
        if str(key) in symbol_name_set
    }
    value = num_expr(**symbol_values)
    if owner.logger is not None:
        owner.logger.info(
            "Evaluating: expression \n\t-> {} \n    with input \t -> [{}] "
            "\n    Output \t\t-> {}".format(
                expr,
                ", ".join("{} : {}".format(key, val) for key, val in symbol_values.items()),
                value,
            )
        )
    return value
