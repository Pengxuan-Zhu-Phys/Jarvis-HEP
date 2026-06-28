#!/usr/bin/env python3
"""In-process Operas backend for Jarvis-HEP V2 Workers."""

from __future__ import annotations

import asyncio
import importlib
import inspect
import queue
import threading
from collections.abc import Callable, Mapping
from typing import Any

import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify


def _resolve_dotted_callable(path: str) -> Callable[..., Any]:
    module_path, _, attr = str(path).rpartition(".")
    if not module_path or not attr:
        raise ValueError(f"invalid operator path: {path}")
    module = importlib.import_module(module_path)
    target = getattr(module, attr)
    if not callable(target):
        raise TypeError(f"operator '{path}' is not callable")
    return target


def _resolve_entry(data: Mapping[str, Any], entry: str) -> Any:
    current: Any = data
    for part in entry.split("."):
        if not isinstance(current, Mapping) or part not in current:
            return None
        current = current[part]
    return current


class OperasModule:
    """Lightweight Operas executor with Worker-side preload."""

    def __init__(self, name: str, config: Mapping[str, Any]) -> None:
        self.name = str(name)
        self.config = dict(config)
        self.operator = str(config["operator"])
        self.input = list(config.get("input", []) or [])
        self.output = list(config.get("output", []) or [])
        self.kwargs = dict(config.get("kwargs", {}) or {})
        self.call_mode = str(config.get("call_mode", "call")).strip().lower()
        self.timeout_sec = self._normalize_timeout(
            config.get("timeout_sec", config.get("timeout"))
        )
        self._func: Callable[..., Any] | None = None
        self._parse_locals, self._numeric_modules = self._build_expression_context()

    @staticmethod
    def _normalize_timeout(value: Any) -> float | None:
        if value is None:
            return None
        timeout = float(value)
        return timeout if timeout > 0 else None

    @staticmethod
    def _build_expression_context() -> tuple[dict[str, Any], dict[str, Any]]:
        parse_locals = {
            name: sp.Symbol(name)
            for name in (
                "x",
                "y",
                "z",
                "shift",
                "LogL",
            )
        }
        numeric_modules = {"sin": np.sin, "cos": np.cos, "exp": np.exp, "log": np.log}
        return parse_locals, numeric_modules

    def preload(self) -> None:
        """Import and cache the operator once per Worker."""
        if self._func is not None:
            return
        self._func = _resolve_dotted_callable(self.operator)

    def _build_input_observables(self, observables: Mapping[str, Any]) -> dict[str, Any]:
        payload = dict(observables)
        if not self.input:
            return payload
        for item in self.input:
            if isinstance(item, str):
                payload[item] = observables.get(item)
                continue
            if not isinstance(item, Mapping) or "name" not in item:
                continue
            target_name = str(item["name"])
            if "expression" in item and isinstance(item["expression"], str):
                expr = sp.sympify(item["expression"], locals=self._parse_locals)
                free_symbols = [str(sym) for sym in expr.free_symbols]
                num_expr = lambdify(free_symbols, expr, modules=[self._numeric_modules, "numpy"])
                symbol_values = {name: payload[name] for name in free_symbols if name in payload}
                missing = [name for name in free_symbols if name not in symbol_values]
                if missing:
                    raise KeyError(
                        f"Operas input expression for '{target_name}' misses observables: {missing}"
                    )
                payload[target_name] = num_expr(**symbol_values)
            else:
                src_key = str(item.get("entry", target_name))
                payload[target_name] = _resolve_entry(payload, src_key)
        return payload

    @staticmethod
    def _run_coro(coro: Any) -> Any:
        try:
            asyncio.get_running_loop()
        except RuntimeError:
            return asyncio.run(coro)
        loop = asyncio.new_event_loop()
        try:
            return loop.run_until_complete(coro)
        finally:
            loop.close()

    def _run_with_timeout(self, callable_obj: Callable[[], Any], label: str) -> Any:
        if self.timeout_sec is None:
            return callable_obj()
        result_queue: queue.Queue[tuple[bool, Any]] = queue.Queue(maxsize=1)

        def _target() -> None:
            try:
                result_queue.put((True, callable_obj()))
            except BaseException as exc:
                result_queue.put((False, exc))

        worker = threading.Thread(
            target=_target,
            name=f"Jarvis2OperasTimeout:{label}",
            daemon=True,
        )
        worker.start()
        try:
            ok, value = result_queue.get(timeout=self.timeout_sec)
        except queue.Empty as exc:
            raise TimeoutError(
                f"Operas call timed out after {self.timeout_sec:g}s -> {label}"
            ) from exc
        if ok:
            return value
        raise value

    def execute(self, observables: Mapping[str, Any], sample_info: Mapping[str, Any]) -> dict[str, Any]:
        """Evaluate the operator and return mapped output observables."""
        if self._func is None:
            self.preload()
        assert self._func is not None

        input_observables = self._build_input_observables(observables)
        call_kwargs = dict(self.kwargs)
        call_kwargs["observables"] = input_observables
        for key, value in input_observables.items():
            call_kwargs.setdefault(str(key), value)

        if self.call_mode == "acall":
            result = self._run_with_timeout(
                lambda: self._run_coro(self._func(**call_kwargs)),
                self.operator,
            )
        else:
            result = self._run_with_timeout(
                lambda: self._func(**call_kwargs),
                self.operator,
            )

        if not isinstance(result, Mapping):
            raise TypeError(
                f"Operas module '{self.name}' requires dict output, got {type(result)}"
            )

        mapped: dict[str, Any] = {}
        output_specs = [
            spec for spec in self.output if isinstance(spec, Mapping) and "name" in spec
        ]
        for spec in output_specs:
            out_name = str(spec["name"])
            src_entry = str(spec.get("entry", out_name))
            mapped[out_name] = _resolve_entry(result, src_entry)
        return mapped


def preload_operas(modules: Mapping[str, Mapping[str, Any]]) -> dict[str, OperasModule]:
    """Build and preload all configured Operas modules."""
    loaded: dict[str, OperasModule] = {}
    for name, config in modules.items():
        module = OperasModule(name, config)
        module.preload()
        loaded[name] = module
    return loaded


__all__ = ["OperasModule", "preload_operas"]