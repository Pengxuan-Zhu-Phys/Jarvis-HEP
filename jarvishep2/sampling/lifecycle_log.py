"""Lifecycle log lines every sampler must emit (docs/logging-spec.md §2).

The runtime calls these at fixed points of a scan, so built-in and plug-in
samplers produce the same start / settings / ready / checkpoint / result /
summary records under ``Jarvis-HEP.Sampler.<Method>`` in ``sampler.log``.

A sampler customises the content through three optional hooks:

``log_settings_rows() -> list[tuple[str, object]]``
    Method-specific settings appended to the common settings table. Without
    this hook the table lists ``Sampling.Bounds`` as written in the card.
``log_summary_rows() -> list[tuple[str, object]]``
    Method-specific results appended to the common summary table.
``log_stop_reason() -> str | None``
    Why the sampler stopped, e.g. ``"converged (dlogz < 0.5)"``.
"""

from __future__ import annotations

import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.log_kv import format_duration, format_two_column_log
from jarvishep2.logging import get_jarvis_logger

_MAX_BOUNDS_ROWS = 24


def _method(sampler: Any) -> str:
    return str(getattr(sampler, "method", "") or type(sampler).__name__).strip()


def sampler_logger(sampler: Any) -> Any:
    """Logger labelled ``Jarvis-HEP.Sampler.<Method>`` for *sampler*."""
    method = _method(sampler)
    return get_jarvis_logger(f"sampler.{method.lower()}", module=f"Sampler.{method}")


def _variable_names(sampler: Any, config: Mapping[str, Any]) -> list[str]:
    names: list[str] = []
    for var in getattr(sampler, "vars", None) or []:
        name = getattr(var, "name", None)
        if name:
            names.append(str(name))
    if names:
        return names
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    for item in sampling.get("Variables") or []:
        if isinstance(item, Mapping) and item.get("name"):
            names.append(str(item["name"]))
    return names


def _bounds_rows(config: Mapping[str, Any]) -> list[tuple[str, Any]]:
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    bounds = sampling.get("Bounds") if isinstance(sampling.get("Bounds"), Mapping) else {}
    rows: list[tuple[str, Any]] = []
    for key, value in bounds.items():
        if str(key) == "seed":
            continue
        if isinstance(value, Mapping):
            for sub_key, sub_value in value.items():
                if not isinstance(sub_value, (Mapping, list, tuple)):
                    rows.append((f"{key}.{sub_key}", sub_value))
        elif not isinstance(value, (list, tuple)):
            rows.append((str(key), value))
    if len(rows) > _MAX_BOUNDS_ROWS:
        hidden = len(rows) - _MAX_BOUNDS_ROWS
        rows = rows[:_MAX_BOUNDS_ROWS] + [("…", f"{hidden} more in Sampling.Bounds")]
    return rows


def _hook_rows(sampler: Any, name: str, logger: Any) -> list[tuple[str, Any]]:
    hook = getattr(sampler, name, None)
    if not callable(hook):
        return []
    try:
        rows = hook() or []
    except Exception:
        logger.exception("%s Sampler could not report %s", _method(sampler), name)
        return []
    return [(str(key), value) for key, value in rows]


class SamplerLifecycleLog:
    """Emit the required lifecycle records for one sampler and one run."""

    def __init__(self, sampler: Any, *, config: Mapping[str, Any] | None = None) -> None:
        self.sampler = sampler
        self.config = dict(config or getattr(sampler, "config", {}) or {})
        self.method = _method(sampler)
        self.logger = sampler_logger(sampler)
        self._t0 = time.time()
        self._result_logged = False

    # S1 + S2
    def start(self, *, resume: bool, workers: int | None = None) -> None:
        self._t0 = time.time()
        self.logger.warning("Initializing the %s Sampling", self.method)
        names = _variable_names(self.sampler, self.config)
        seed = getattr(self.sampler, "_seed", None)
        rows: list[tuple[str, Any]] = [
            ("method", self.method),
            ("seed", "—" if seed is None else seed),
            ("variables", f"{len(names)} ({', '.join(names)})" if names else "0"),
        ]
        if workers:
            rows.append(("workers", workers))
        rows.append(("resume", "yes" if resume else "no"))
        # A sampler that lists its own settings knows which Bounds matter;
        # otherwise (e.g. a plug-in without the hook) show Bounds as written.
        if callable(getattr(self.sampler, "log_settings_rows", None)):
            rows.extend(_hook_rows(self.sampler, "log_settings_rows", self.logger))
        else:
            rows.extend(_bounds_rows(self.config))
        self.logger.info(format_two_column_log(f"{self.method} Sampler Settings ->", rows))

    # S3
    def ready(self) -> None:
        self.logger.warning("WorkerFactory is ready for %s sampler", self.method)

    # S6
    def checkpoint_saved(self, path: str, reason: str = "") -> None:
        suffix = f" ({reason})" if reason else ""
        self.logger.info("Checkpoint saved to %s%s", path, suffix)

    def checkpoint_loaded(self, path: str) -> None:
        self.logger.info("Checkpoint loaded from %s", path)

    # S8
    def failure(self, doing: str, exc: BaseException) -> None:
        self.logger.error(
            "%s Sampler meets error when %s -> %s",
            self.method,
            doing,
            exc,
            exc_info=(type(exc), exc, exc.__traceback__),
        )

    # S9 + S10
    def result(self, outcome: Any, *, stop_reason: str | None = None) -> None:
        if self._result_logged:
            return
        self._result_logged = True
        elapsed = format_duration(time.time() - self._t0)
        completed = int(getattr(outcome, "completed", 0) or 0)
        failed = int(getattr(outcome, "failed", 0) or 0)
        submitted = int(getattr(outcome, "submitted", 0) or 0)
        self.logger.warning(
            "%s Sampler obtains %d samples in %s", self.method, completed + failed, elapsed
        )
        reason = stop_reason or self._stop_reason(outcome)
        rows: list[tuple[str, Any]] = [
            ("stop reason", reason),
            ("submitted", submitted),
            ("completed", completed),
            ("failed", failed),
            ("elapsed", elapsed),
        ]
        rows.extend(_hook_rows(self.sampler, "log_summary_rows", self.logger))
        self.logger.warning(format_two_column_log(f"{self.method} Sampler Summary ->", rows))

    def _stop_reason(self, outcome: Any) -> str:
        status = str(getattr(outcome, "status", "") or "")
        if status == "interrupted":
            return "interrupted (continue with --resume)"
        hook = getattr(self.sampler, "log_stop_reason", None)
        if callable(hook):
            try:
                reason = hook()
            except Exception:
                self.logger.exception("%s Sampler could not report its stop reason", self.method)
                reason = None
            if reason:
                return str(reason)
        return {
            "success": "all samples evaluated",
            "partial_failure": "all samples evaluated, some failed",
            "failed": "all samples failed",
            "error": "stopped by an error",
        }.get(status, status or "unknown")


__all__ = ["SamplerLifecycleLog", "sampler_logger"]
