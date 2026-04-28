#!/usr/bin/env python3
from __future__ import annotations

import asyncio
import inspect
import queue
import threading
from collections.abc import Mapping
from typing import Any

from jarvishep.Module.module import Module


class _OperasSampleLoggerBridge:
    """Forward sample logs while silencing low-signal Operas debug chatter."""

    def __init__(self, logger):
        self._logger = logger
        self._options = getattr(logger, "_options", None)

    def bind(self, **extra):
        return type(self)(self._logger.bind(**extra))

    def log(self, level, message, *args, **kwargs):
        if str(level).strip().upper() == "DEBUG":
            return None
        return self._logger.log(level, message, *args, **kwargs)

    def debug(self, *_args, **_kwargs):
        return None

    def info(self, *args, **kwargs):
        return self._logger.info(*args, **kwargs)

    def warning(self, *args, **kwargs):
        return self._logger.warning(*args, **kwargs)

    def error(self, *args, **kwargs):
        return self._logger.error(*args, **kwargs)

    def critical(self, *args, **kwargs):
        return self._logger.critical(*args, **kwargs)


class OperasModule(Module):
    def __init__(self, name, config):
        super().__init__(name)
        self.config = config
        self.type = "Operas"
        self.required_modules = config.get("required_modules", []) or []
        self.operator = config["operator"]
        self.input = config.get("input", []) or []
        self.output = config.get("output", []) or []
        self.kwargs = config.get("kwargs", {}) or {}
        self.call_mode = self._normalize_call_mode(config.get("call_mode", "call"))
        self.timeout_sec = self._normalize_timeout(
            config.get("timeout_sec", config.get("timeout"))
        )
        self._funcs = {}
        self._registry = None
        self._init_io_specs()

    @property
    def funcs(self):
        return self._funcs

    def set_funcs(self, funcs):
        self._funcs = funcs

    def set_logger(self, logger):
        self.logger = logger

    @staticmethod
    def _normalize_call_mode(mode):
        normalized = str(mode).strip().lower()
        if normalized not in {"call", "acall"}:
            raise ValueError(f"Unsupported Operas call_mode '{mode}'. Expected 'call' or 'acall'.")
        return normalized

    @staticmethod
    def _normalize_timeout(value):
        if value is None:
            return None
        try:
            timeout = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"Unsupported Operas timeout_sec '{value}'. Expected a positive number or null.")
        if timeout <= 0:
            return None
        return timeout

    def _init_io_specs(self):
        import sympy as sp
        from jarvishep.inner_func import build_expression_context, update_const, update_funcs

        parse_locals, _ = build_expression_context(
            funcs=update_funcs({}),
            consts=update_const({}),
        )
        for item in self.input:
            if isinstance(item, str):
                self.inputs[item] = None
                continue
            if isinstance(item, dict) and "name" in item:
                if "expression" in item and isinstance(item["expression"], str):
                    expr = sp.sympify(item["expression"], locals=parse_locals)
                    item["_inc"] = [str(sym) for sym in expr.free_symbols]
                    for dep in item["_inc"]:
                        self.inputs[dep] = None
                else:
                    dep_name = item.get("entry", item["name"])
                    self.inputs[str(dep_name)] = None

        for spec in self.output:
            if isinstance(spec, dict) and "name" in spec:
                self.outputs[spec["name"]] = None

    def _get_registry(self):
        if self._registry is None:
            from jarvis_operas import get_global_operas_registry

            self._registry = get_global_operas_registry()
        return self._registry

    @staticmethod
    def _resolve_entry(data: Mapping, entry: str):
        current = data
        for part in entry.split("."):
            if not isinstance(current, Mapping) or part not in current:
                return None
            current = current[part]
        return current

    def _build_input_observables(self, observables, slogger=None):
        import sympy as sp
        from sympy.utilities.lambdify import lambdify
        from jarvishep.inner_func import build_expression_context, update_const

        if not isinstance(observables, Mapping):
            raise TypeError(
                f"Operas module '{self.name}' expects observables mapping, got {type(observables)}"
            )

        payload = dict(observables)
        if not self.input:
            return payload

        parse_locals, numeric_modules = build_expression_context(
            funcs=dict(self.funcs or {}),
            consts=update_const({}),
        )

        for item in self.input:
            if isinstance(item, str):
                payload[item] = observables.get(item)
                continue
            if not isinstance(item, dict) or "name" not in item:
                continue

            target_name = item["name"]
            if "expression" in item and isinstance(item["expression"], str):
                expr = sp.sympify(item["expression"], locals=parse_locals)
                free_symbols = [str(sym) for sym in expr.free_symbols]
                free_symbol_set = set(free_symbols)
                num_expr = lambdify(free_symbols, expr, modules=[numeric_modules, "numpy"])
                symbol_values_strs = {
                    str(k): v for k, v in payload.items() if str(k) in free_symbol_set
                }
                missing = [name for name in free_symbols if name not in symbol_values_strs]
                if missing:
                    raise KeyError(
                        f"Operas input expression for '{target_name}' misses observables: {missing}"
                    )
                value = num_expr(**symbol_values_strs)
                payload[target_name] = value
                target_logger = slogger if slogger is not None else getattr(self, "logger", None)
                if target_logger is not None:
                    target_logger.info(
                        f"Evaluating   {target_name}: \n   expression \t-> {str(item['expression'])} \n   with input \t-> [{', '.join(['{} : {}'.format(kk, vv) for kk, vv in symbol_values_strs.items()])}] \n   Output \t\t-> {value}"
                    )
            else:
                src_key = str(item.get("entry", target_name))
                payload[target_name] = self._resolve_entry(payload, src_key)
        return payload

    @staticmethod
    def _run_coro(coro):
        try:
            asyncio.get_running_loop()
        except RuntimeError:
            return asyncio.run(coro)

        # Fallback for environments with a running loop in current thread.
        loop = asyncio.new_event_loop()
        try:
            return loop.run_until_complete(coro)
        finally:
            loop.close()

    @staticmethod
    def _run_with_timeout(callable_obj, timeout_sec, timeout_label):
        if timeout_sec is None:
            return callable_obj()

        result_queue = queue.Queue(maxsize=1)

        def _target():
            try:
                result_queue.put((True, callable_obj()))
            except BaseException as exc:
                result_queue.put((False, exc))

        worker = threading.Thread(
            target=_target,
            name=f"JarvisOperasTimeoutGuard:{timeout_label}",
            daemon=True,
        )
        worker.start()
        try:
            ok, value = result_queue.get(timeout=timeout_sec)
        except queue.Empty:
            raise TimeoutError(
                f"Operas call timed out after {timeout_sec:g}s -> {timeout_label}"
            )
        if ok:
            return value
        raise value

    def _resolve_operator_signature(self, registry):
        resolved_operator = registry.resolve_name(self.operator)
        declaration = registry.get(resolved_operator)
        fn = declaration
        if not callable(fn):
            candidate = getattr(declaration, "numpy_impl", None)
            if callable(candidate):
                fn = candidate

        try:
            sig = inspect.signature(fn)
        except (TypeError, ValueError):
            # Unknown signature: do not filter kwargs.
            return resolved_operator, None, True

        accepts_var_kwargs = False
        accepted_kwargs = set()
        for param in sig.parameters.values():
            if param.kind == inspect.Parameter.VAR_KEYWORD:
                accepts_var_kwargs = True
            elif param.kind in (
                inspect.Parameter.POSITIONAL_OR_KEYWORD,
                inspect.Parameter.KEYWORD_ONLY,
            ):
                accepted_kwargs.add(param.name)

        return resolved_operator, accepted_kwargs, accepts_var_kwargs

    @staticmethod
    def _filter_call_kwargs(call_kwargs, accepted_kwargs, accepts_var_kwargs):
        if accepted_kwargs is None or accepts_var_kwargs:
            return dict(call_kwargs), []

        filtered = {k: v for k, v in call_kwargs.items() if k in accepted_kwargs}
        dropped = [k for k in call_kwargs.keys() if k not in accepted_kwargs]
        return filtered, dropped

    @staticmethod
    def _normalize_log_value(value: Any):
        try:
            import numpy as np
        except Exception:
            np = None

        if np is not None and isinstance(value, np.generic):
            return value.item()
        if np is not None and isinstance(value, np.ndarray):
            return [OperasModule._normalize_log_value(v) for v in value.tolist()]
        if isinstance(value, Mapping):
            return {
                str(k): OperasModule._normalize_log_value(v)
                for k, v in value.items()
            }
        if isinstance(value, (list, tuple)):
            return [OperasModule._normalize_log_value(v) for v in value]
        return value

    @staticmethod
    def _format_inline_pairs(values: Mapping[str, Any]) -> str:
        rendered = []
        for key, value in values.items():
            normalized = OperasModule._normalize_log_value(value)
            rendered.append(f"{key} : {normalized}")
        return "[" + ", ".join(rendered) + "]"

    @staticmethod
    def _format_dispatch_message(module_name: str, operator: str, call_mode: str) -> str:
        return (
            "Operas input dispatch:\n"
            f"   module \t-> {module_name}\n"
            f"   operator \t-> {operator}\n"
            f"   call_mode \t-> {call_mode}"
        )

    def execute(self, observables, sample_info):
        slogger = sample_info.get("logger", None) if isinstance(sample_info, dict) else None
        registry = self._get_registry()
        input_observables = self._build_input_observables(observables, slogger=slogger)
        call_kwargs = dict(self.kwargs)
        call_kwargs["observables"] = input_observables
        if isinstance(input_observables, Mapping):
            for key, value in input_observables.items():
                call_kwargs.setdefault(str(key), value)

        resolved_operator, accepted_kwargs, accepts_var_kwargs = self._resolve_operator_signature(
            registry
        )
        call_kwargs, dropped_keys = self._filter_call_kwargs(
            call_kwargs,
            accepted_kwargs,
            accepts_var_kwargs,
        )
        if dropped_keys:
            target_logger = slogger if slogger is not None else getattr(self, "logger", None)
            if target_logger is not None:
                target_logger.info(
                    f"Operas module '{self.name}' filtered unsupported kwargs for '{resolved_operator}': {dropped_keys}"
                )

        if slogger is not None:
            slogger.info(
                self._format_dispatch_message(
                    self.name,
                    resolved_operator,
                    self.call_mode,
                )
            )
            normalized_observables = self._normalize_log_value(input_observables)
            slogger.info(
                "Operas input observables:\n"
                f"   with input \t-> {self._format_inline_pairs(normalized_observables)}"
            )
            extra_kwargs = {
                key: value
                for key, value in call_kwargs.items()
                if key not in {"observables", "logger"}
            }
            normalized_kwargs = self._normalize_log_value(extra_kwargs)
            if normalized_kwargs and normalized_kwargs != normalized_observables:
                slogger.info(
                    "Operas input kwargs:\n"
                    f"   with kwargs \t-> {self._format_inline_pairs(normalized_kwargs)}"
                )

        registry_logger = (
            _OperasSampleLoggerBridge(slogger)
            if slogger is not None
            else None
        )
        if registry_logger is not None:
            call_kwargs.setdefault("sample_logger", registry_logger)
            call_kwargs, _ = self._filter_call_kwargs(
                call_kwargs,
                accepted_kwargs,
                accepts_var_kwargs,
            )

        if self.call_mode == "acall":
            result = self._run_with_timeout(
                lambda: self._run_coro(
                    registry.acall(resolved_operator, logger=registry_logger, **call_kwargs)
                ),
                self.timeout_sec,
                resolved_operator,
            )
        else:
            result = self._run_with_timeout(
                lambda: registry.call(resolved_operator, logger=registry_logger, **call_kwargs),
                self.timeout_sec,
                resolved_operator,
            )
        output_specs = [
            spec for spec in self.output if isinstance(spec, dict) and "name" in spec
        ]

        if not isinstance(result, Mapping):
            raise TypeError(
                f"Operas module '{self.name}' requires dict output, got {type(result)}"
            )

        mapped = {}
        for spec in output_specs:
            out_name = spec["name"]
            src_entry = spec.get("entry", out_name)
            mapped[out_name] = self._resolve_entry(result, src_entry)
        return mapped
