"""Jarvis-HEP V2 — distributed runtime package (independent from jarvishep V1)."""

from __future__ import annotations

from importlib import import_module
from typing import Any

from jarvishep2._legacy_paths import install_legacy_paths

install_legacy_paths()

from jarvishep2.expression import (
    CompiledExpression,
    ExpressionContext,
    MissingExpressionVariablesError,
)
from jarvishep2.logging import (
    BufferedSampleLogger,
    SampleLogger,
    get_jarvis_logger,
    setup_jarvis_logging,
    shutdown_jarvis_logging,
)
from jarvishep2.sample import (
    ExecutionStep,
    Sample,
    UMapperProtocol,
    VALID_EXECUTION_STEP_TYPES,
    ensure_sample_materialized,
    materialize_failure_artifacts,
)

__version__ = "2.0.6"

_LAZY_ATTRS: dict[str, tuple[str, str]] = {
    "Jarvis2Core": ("jarvishep2.runtime.core", "Jarvis2Core"),
    "RedisQueue": ("jarvishep2.queue.redis_queue", "RedisQueue"),
    "TaskFactory": ("jarvishep2.runtime.factory", "TaskFactory"),
    "TaskValidationError": ("jarvishep2.queue.redis_queue", "TaskValidationError"),
    "Worker": ("jarvishep2.runtime.worker", "Worker"),
    "make_fakeredis_queue": ("jarvishep2.queue.redis_queue", "make_fakeredis_queue"),
}

__all__ = [
    "BufferedSampleLogger",
    "CompiledExpression",
    "ExecutionStep",
    "ExpressionContext",
    "Jarvis2Core",
    "MissingExpressionVariablesError",
    "RedisQueue",
    "Sample",
    "SampleLogger",
    "TaskFactory",
    "TaskValidationError",
    "UMapperProtocol",
    "VALID_EXECUTION_STEP_TYPES",
    "Worker",
    "ensure_sample_materialized",
    "get_jarvis_logger",
    "make_fakeredis_queue",
    "materialize_failure_artifacts",
    "setup_jarvis_logging",
    "shutdown_jarvis_logging",
]


def __getattr__(name: str) -> Any:
    target = _LAZY_ATTRS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr = target
    value = getattr(import_module(module_name), attr)
    globals()[name] = value
    return value
