#!/usr/bin/env python3
"""In-process Operas backend for Jarvis-HEP V2 Workers.

Operator resolution prefers the Jarvis-Operas registry when the name is
registered there, then falls back to an importlib dotted callable so local
test operators keep working without the Operas package.
"""

from __future__ import annotations

import asyncio
import importlib
import inspect
import time
import queue
import threading
from collections.abc import Callable, Mapping
from typing import Any

from jarvishep2.expression import (
    CompiledExpression,
    ExpressionContext,
    MissingExpressionVariablesError,
)


def _resolve_dotted_callable(path: str) -> Callable[..., Any]:
    module_path, _, attr = str(path).rpartition(".")
    if not module_path or not attr:
        raise ValueError(f"invalid operator path: {path}")
    module = importlib.import_module(module_path)
    target = getattr(module, attr)
    if not callable(target):
        raise TypeError(f"operator '{path}' is not callable")
    return target


def _try_jarvis_operas_registry(operator: str) -> Any | None:
    """Return the global Operas registry if *operator* resolves; else None."""
    try:
        from jarvishep2.operas_functions import ensure_operas_registry_discovered
    except ImportError:
        return None
    registry = ensure_operas_registry_discovered()
    try:
        registry.resolve_name(operator)
    except Exception:
        return None
    return registry


def _snapshot_operas_registry_operator(
    registry: Any,
    operator: str,
    *,
    call_mode: str,
) -> Callable[..., Any]:
    """Resolve the NumPy implementation once, for direct Worker-side reuse."""
    declaration = registry.get(operator)
    func = getattr(declaration, "numpy_impl", None)
    if not callable(func):
        raise TypeError(
            f"Jarvis-Operas operator {operator!r} has no callable NumPy implementation"
        )
    if str(call_mode).strip().lower() != "acall":
        return func

    async def _acall(**kwargs: Any) -> Any:
        result = func(**kwargs)
        if inspect.isawaitable(result):
            return await result
        return result

    return _acall


_VALID_CALL_MODES = frozenset({"call", "acall"})


def normalize_call_mode(value: Any) -> str:
    """Strict ``call`` / ``acall`` validation (D11.4)."""
    mode = str(value if value is not None else "call").strip().lower()
    if mode not in _VALID_CALL_MODES:
        raise ValueError(
            f"Operas call_mode must be one of {sorted(_VALID_CALL_MODES)}, got {value!r}"
        )
    return mode


def filter_operator_kwargs(
    func: Callable[..., Any],
    kwargs: Mapping[str, Any],
    *,
    operator_name: str = "",
) -> dict[str, Any]:
    """Pass only parameters accepted by *func* (unless it has ``**kwargs``).

    Declared parameters missing from *kwargs* raise a friendly error instead of
    a bare TypeError from the call site.
    """
    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):
        return dict(kwargs)

    params = signature.parameters
    accepts_var_keyword = any(
        param.kind == inspect.Parameter.VAR_KEYWORD for param in params.values()
    )
    if accepts_var_keyword:
        return dict(kwargs)

    allowed = {
        name
        for name, param in params.items()
        if param.kind
        in (
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.KEYWORD_ONLY,
        )
    }
    filtered = {key: value for key, value in kwargs.items() if key in allowed}

    missing: list[str] = []
    for name, param in params.items():
        if param.kind in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        ):
            continue
        if param.default is not inspect.Parameter.empty:
            continue
        if name in filtered:
            continue
        missing.append(name)
    if missing:
        label = f" '{operator_name}'" if operator_name else ""
        raise TypeError(
            f"Operas operator{label} missing required parameter(s): {missing}. "
            f"Available inputs: {sorted(kwargs)}"
        )
    return filtered


def resolve_operator(operator: str, *, call_mode: str = "call") -> Callable[..., Any]:
    """Resolve an operator name via importlib, then Jarvis-Operas.

    Resolution order:
    1. importlib dotted path ``package.module.callable`` (fast; no Operas bootstrap)
    2. Jarvis-Operas global registry (if installed and name is registered)

    Importlib is tried first so Worker preload of local test operators does not
    pay the cost of bootstrapping the full Operas catalog.
    """
    name = str(operator).strip()
    if not name:
        raise ValueError("Operas operator name must not be empty")
    mode = normalize_call_mode(call_mode)

    import_error: Exception | None = None
    try:
        func = _resolve_dotted_callable(name)
        if mode == "acall":
            async def _acall(**kwargs: Any) -> Any:
                result = func(**kwargs)
                if inspect.isawaitable(result):
                    return await result
                return result

            return _acall
        return func
    except Exception as exc:
        import_error = exc

    registry = _try_jarvis_operas_registry(name)
    if registry is not None:
        return _snapshot_operas_registry_operator(registry, name, call_mode=mode)

    try:
        import jarvis_operas  # noqa: F401
        operas_hint = (
            "Not found in Jarvis-Operas registry either. "
            "Use a registered name (e.g. helper.eggbox2d) or a Python dotted path."
        )
    except ImportError:
        operas_hint = (
            "Jarvis-Operas is not installed; registry names are unavailable. "
            "Install with `pip install 'jarvishep2[operas]'` or use a Python dotted path."
        )
    raise ValueError(
        f"Cannot resolve Operas operator {name!r}: {import_error}. {operas_hint}"
    ) from import_error


def _resolve_entry(data: Mapping[str, Any], entry: str) -> Any:
    current: Any = data
    for part in entry.split("."):
        if not isinstance(current, Mapping) or part not in current:
            return None
        current = current[part]
    return current


class OperasModule:
    """Lightweight Operas executor with Worker-side preload."""

    def __init__(
        self,
        name: str,
        config: Mapping[str, Any],
        *,
        expression_context: ExpressionContext | None = None,
    ) -> None:
        self.name = str(name)
        self.config = dict(config)
        self.operator = str(config["operator"])
        self.input = list(config.get("input", []) or [])
        self.output = list(config.get("output", []) or [])
        self.kwargs = dict(config.get("kwargs", {}) or {})
        self.call_mode = normalize_call_mode(config.get("call_mode", "call"))
        self.timeout_sec = self._normalize_timeout(
            config.get("timeout_sec", config.get("timeout"))
        )
        self._func: Callable[..., Any] | None = None
        self._bound_signature: inspect.Signature | None = None
        if expression_context is None:
            from jarvishep2.operas_functions import (
                build_operas_expression_context,
                expression_uses_operas_function,
            )

            expression_context = (
                build_operas_expression_context(required=True)
                if expression_uses_operas_function(self.input)
                else ExpressionContext()
            )
        self._expression_context = expression_context
        self._compiled_input_expressions: dict[int, CompiledExpression] = {}
        self._input_expressions_compiled = False

    @staticmethod
    def _normalize_timeout(value: Any) -> float | None:
        if value is None:
            return None
        timeout = float(value)
        return timeout if timeout > 0 else None

    def preload(self) -> None:
        """Compile input expressions and resolve the operator once per Worker."""
        # Re-validate in case config was mutated after construction.
        self.call_mode = normalize_call_mode(self.call_mode)
        self._compile_input_expressions()
        if self._func is None:
            self._func = resolve_operator(self.operator, call_mode=self.call_mode)
            try:
                self._bound_signature = inspect.signature(self._func)
            except (TypeError, ValueError):
                self._bound_signature = None

    def _compile_input_expressions(self) -> None:
        if self._input_expressions_compiled:
            return
        compiled: dict[int, CompiledExpression] = {}
        for index, item in enumerate(self.input):
            if not isinstance(item, Mapping):
                continue
            expression = item.get("expression")
            if not isinstance(expression, str):
                continue
            compiled[index] = self._expression_context.compile(
                expression,
                symbols=("x", "y", "z", "shift", "LogL"),
            )
        self._compiled_input_expressions = compiled
        self._input_expressions_compiled = True

    @staticmethod
    def _normalize_log_value(value: Any) -> Any:
        """Convert NumPy containers/scalars into readable sample-log values."""
        try:
            import numpy as np
        except ImportError:  # pragma: no cover - NumPy is a runtime dependency
            np = None
        if np is not None and isinstance(value, np.generic):
            return value.item()
        if np is not None and isinstance(value, np.ndarray):
            return [OperasModule._normalize_log_value(item) for item in value.tolist()]
        if isinstance(value, Mapping):
            return {str(key): OperasModule._normalize_log_value(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [OperasModule._normalize_log_value(item) for item in value]
        return value

    @classmethod
    def _format_inline_pairs(cls, values: Mapping[str, Any]) -> str:
        return "[" + ", ".join(
            f"{key} : {cls._normalize_log_value(value)}" for key, value in values.items()
        ) + "]"

    def _sample_logger(self, sample_info: Mapping[str, Any]) -> Any | None:
        logger = sample_info.get("logger") if isinstance(sample_info, Mapping) else None
        if logger is None:
            return None
        bind = getattr(logger, "bind", None)
        if callable(bind):
            logger_name = str(sample_info.get("logger_name") or "Sample")
            return bind(module=f"{logger_name} (Operas:{self.name})")
        return logger

    def _build_input_observables(
        self,
        observables: Mapping[str, Any],
        *,
        logger: Any | None = None,
    ) -> dict[str, Any]:
        payload = dict(observables)
        if not self.input:
            return payload
        for index, item in enumerate(self.input):
            if isinstance(item, str):
                payload[item] = observables.get(item)
                continue
            if not isinstance(item, Mapping) or "name" not in item:
                continue
            target_name = str(item["name"])
            if "expression" in item and isinstance(item["expression"], str):
                compiled = self._compiled_input_expressions.get(index)
                if compiled is None:
                    raise RuntimeError(
                        f"Operas input expression for '{target_name}' was not precompiled"
                    )
                try:
                    payload[target_name] = compiled.evaluate(payload)
                except MissingExpressionVariablesError as exc:
                    raise KeyError(
                        f"Operas input expression for '{target_name}' misses observables: "
                        f"{list(exc.missing)}"
                    ) from exc
                if logger is not None:
                    expression_inputs = {
                        name: payload[name]
                        for name in compiled.variable_names
                        if name in payload
                    }
                    logger.info(
                        "Evaluating   {}:\n"
                        "   expression \t-> {}\n"
                        "   with input \t-> {}\n"
                        "   Output \t\t-> {}",
                        target_name,
                        item["expression"],
                        self._format_inline_pairs(expression_inputs),
                        self._normalize_log_value(payload[target_name]),
                    )
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

    def execute(
        self,
        observables: Mapping[str, Any],
        sample_info: Mapping[str, Any],
        *,
        sample_logger: Any | None = None,
    ) -> dict[str, Any]:
        """Evaluate the operator with the Worker-provided Sample logger.

        ``sample_logger`` is an explicit execution dependency, rather than a
        module-global logger.  It is forwarded to Jarvis-Operas NumPy
        implementations under their established ``logger`` keyword.
        """
        if self._func is None:
            self.preload()
        assert self._func is not None

        logger = sample_logger if sample_logger is not None else self._sample_logger(sample_info)
        started_at = time.monotonic()
        try:
            input_observables = self._build_input_observables(observables, logger=logger)
            raw_kwargs: dict[str, Any] = dict(self.kwargs)
            raw_kwargs["observables"] = input_observables
            if logger is not None:
                raw_kwargs.setdefault("logger", logger)
            for key, value in input_observables.items():
                raw_kwargs.setdefault(str(key), value)
            call_kwargs = filter_operator_kwargs(
                self._func,
                raw_kwargs,
                operator_name=self.operator,
            )

            if logger is not None:
                logger.info(
                    "Operas input dispatch:\n"
                    "   module \t-> {}\n"
                    "   operator \t-> {}\n"
                    "   call_mode \t-> {}",
                    self.name,
                    self.operator,
                    self.call_mode,
                )
                logger.info(
                    "Operas input observables:\n"
                    "   with input \t-> {}",
                    self._format_inline_pairs(input_observables),
                )
                if self.kwargs:
                    logger.info(
                        "Operas static kwargs:\n"
                        "   with kwargs \t-> {}",
                        self._format_inline_pairs(self.kwargs),
                    )

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
            if logger is not None:
                elapsed_ms = (time.monotonic() - started_at) * 1000.0
                logger.info(
                    "Operas output:\n"
                    "   result \t-> {}\n"
                    "   mapped \t-> {}\n"
                    "   elapsed \t-> {:.3f} ms",
                    self._normalize_log_value(result),
                    self._format_inline_pairs(mapped),
                    elapsed_ms,
                )
            return mapped
        except Exception as exc:
            if logger is not None:
                logger.error(
                    "Operas call failed:\n"
                    "   module \t-> {}\n"
                    "   operator \t-> {}\n"
                    "   error \t-> {}: {}",
                    self.name,
                    self.operator,
                    type(exc).__name__,
                    exc,
                )
            raise


def preload_operas(
    modules: Mapping[str, Mapping[str, Any]],
    *,
    expression_context: ExpressionContext | None = None,
) -> dict[str, OperasModule]:
    """Build and preload all configured Operas modules."""
    loaded: dict[str, OperasModule] = {}
    for name, config in modules.items():
        module = OperasModule(name, config, expression_context=expression_context)
        module.preload()
        loaded[name] = module
    return loaded


__all__ = [
    "OperasModule",
    "filter_operator_kwargs",
    "normalize_call_mode",
    "preload_operas",
    "resolve_operator",
]
