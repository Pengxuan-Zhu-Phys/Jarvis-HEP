#!/usr/bin/env python3
"""Per-Sample logging primitives for Jarvis-HEP V2 (V1 contract, reimplemented)."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from datetime import datetime
import logging
import os
import threading
from typing import Any, Deque, Iterable, Mapping


DEFAULT_BUFFER_MAX_EVENTS = 2048

_LEVEL_NUMBERS = {
    "CRITICAL": logging.CRITICAL,
    "ERROR": logging.ERROR,
    "WARNING": logging.WARNING,
    "WARN": logging.WARNING,
    "INFO": logging.INFO,
    "DEBUG": logging.DEBUG,
    "NOTSET": logging.NOTSET,
}


def _level_number(level: Any) -> int:
    if isinstance(level, int):
        return int(level)
    return _LEVEL_NUMBERS.get(str(level).strip().upper(), logging.INFO)


def _forward_sample_event(
    *,
    extra: Mapping[str, Any] | dict[str, Any],
    module: str,
    level: str,
    message: str,
    raw: bool = False,
) -> None:
    """Forward an eligible sample event to the process terminal logger.

    Sample files remain the source of complete per-sample detail.  The terminal
    path is opt-in and independently thresholded so ordinary runs show only
    sample errors while ``check`` can expose debug output.
    """
    to_console = extra.get("sample_log_to_console", extra.get("to_console", False))
    if not bool(to_console) or bool(extra.get("sample_log_silence", extra.get("log_silence", False))):
        return

    try:
        from jarvishep2.logging.toplevel import (
            console_logging_enabled,
            get_jarvis_logger,
        )

        if not console_logging_enabled():
            return
    except Exception:
        return

    threshold = _level_number(extra.get("sample_console_level", extra.get("console_level", "ERROR")))
    level_no = _level_number(level)
    if level_no < threshold:
        return

    worker_id = extra.get("_sample_worker_id", extra.get("worker_id"))
    try:
        worker_id = int(worker_id) if worker_id is not None else None
    except (TypeError, ValueError):
        worker_id = None
    try:
        logger = get_jarvis_logger("worker", worker_id=worker_id).bind(
            jarvis_module=module
        )
        # Preserve subprocess stdout/stderr passthrough semantics.  The
        # SampleLogger backend already uses ``raw`` to omit the sample header;
        # the terminal forwarding path must carry it through as well.
        if raw:
            logger = logger.bind(raw=True)
        logger.log(level_no, message)
    except Exception:
        # Logging must never make a sample calculation fail.
        return


@dataclass(frozen=True)
class SampleLogEvent:
    timestamp: str
    module: str
    level: str
    message: str
    raw: bool


class _SampleLogBackend:
    """Sample-local file sink with exact legacy sample-log formatting."""

    def __init__(self, path: str, *, time_provider=None) -> None:
        self.path = os.path.abspath(path)
        self._time_provider = time_provider or datetime.now
        self._lock = threading.Lock()
        self._closed = False
        self._has_written = False
        self._last_ended_newline = False
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        self._handle = open(self.path, "a", encoding="utf-8")

    def _timestamp(self) -> str:
        from jarvishep2.logging.style import sample_style

        style = sample_style()
        fmt = str(style.get("date_format") or "%m-%d %H:%M:%S.%f")
        stamp = self._time_provider().strftime(fmt)
        frac = int(style.get("frac_digits") or 0)
        # Historical: "%f" is 6 digits; keep first frac_digits (default 3 → ms).
        if "%f" in fmt and frac > 0 and "." in stamp:
            head, tail = stamp.rsplit(".", 1)
            stamp = f"{head}.{tail[:frac]}"
        return stamp

    def write(
        self,
        *,
        module: str,
        level: str,
        message: str,
        raw: bool,
        timestamp: str | None = None,
    ) -> None:
        if self._closed:
            return

        from jarvishep2.logging.style import sample_style

        style = sample_style()
        text = str(message)
        if not raw:
            prefix = "\n"
            if (
                bool(style.get("double_newline_gap", True))
                and self._has_written
                and not self._last_ended_newline
            ):
                prefix = "\n\n"
            ts = timestamp if timestamp is not None else self._timestamp()
            bullet = str(style.get("bullet") or "·•·")
            head_tmpl = str(
                style.get("head")
                or "\n{bullet} {module} \n\t-> {timestamp} - [{level}] >>> \n"
            )
            # head template already starts with \n; combine with gap prefix.
            head = head_tmpl.format(
                bullet=bullet,
                module=module,
                timestamp=ts,
                level=level,
                message="",
            )
            if prefix == "\n\n" and head.startswith("\n"):
                head = "\n" + head
            elif prefix == "\n" and not self._has_written and head.startswith("\n"):
                pass
            text = f"{head}{text}"
        else:
            if self._has_written and not self._last_ended_newline and not text.startswith("\n"):
                text = "\n" + text
            if not text.endswith("\n"):
                text += "\n"

        with self._lock:
            if self._closed:
                return
            self._handle.write(text)
            self._handle.flush()
            self._has_written = True
            self._last_ended_newline = text.endswith("\n")

    def close(self) -> None:
        with self._lock:
            if self._closed:
                return
            self._closed = True
            self._handle.close()


def replay_sample_log_events(
    path: str,
    events: Iterable[SampleLogEvent],
    *,
    time_provider=None,
) -> None:
    backend = _SampleLogBackend(path, time_provider=time_provider)
    try:
        for event in events:
            backend.write(
                module=event.module,
                level=event.level,
                message=event.message,
                raw=event.raw,
                timestamp=event.timestamp,
            )
    finally:
        backend.close()


class SampleLogger:
    """Sample-local file logger with optional thresholded terminal forwarding."""

    def __init__(self, backend: _SampleLogBackend, *, extra: dict[str, Any] | None = None) -> None:
        self._backend = backend
        self._extra = dict(extra or {})
        self._options = (None, None, None, None, None, None, None, None, self._extra)

    @classmethod
    def open(
        cls,
        path: str,
        *,
        module: str,
        time_provider=None,
        extra: dict[str, Any] | None = None,
    ) -> SampleLogger:
        merged = dict(extra or {})
        merged.setdefault("module", module)
        backend = _SampleLogBackend(path, time_provider=time_provider)
        return cls(backend, extra=merged)

    def bind(self, **extra: Any) -> SampleLogger:
        merged = dict(self._extra)
        merged.update(extra)
        return type(self)(self._backend, extra=merged)

    def close(self) -> None:
        self._backend.close()

    def log(self, level: Any, message: Any, *args: Any, **kwargs: Any) -> None:
        module = str(self._extra.get("module", "No module"))
        raw = "raw" in self._extra
        rendered = self._render_message(message, *args, **kwargs)
        normalized_level = self._normalize_level(level)
        self._backend.write(
            module=module,
            level=normalized_level,
            message=rendered,
            raw=raw,
        )
        _forward_sample_event(
            extra=self._extra,
            module=module,
            level=normalized_level,
            message=rendered,
            raw=raw,
        )

    def debug(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("DEBUG", message, *args, **kwargs)

    def info(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("INFO", message, *args, **kwargs)

    def warning(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("WARNING", message, *args, **kwargs)

    def error(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("ERROR", message, *args, **kwargs)

    def critical(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("CRITICAL", message, *args, **kwargs)

    @staticmethod
    def _normalize_level(level: Any) -> str:
        text = str(level).strip().upper()
        return text or "INFO"

    @staticmethod
    def _render_message(message: Any, *args: Any, **kwargs: Any) -> str:
        text = str(message)
        if args or kwargs:
            return text.format(*args, **kwargs)
        return text


class BufferedSampleLogger:
    """In-memory sample logger used before sample artifacts are materialized."""

    def __init__(
        self,
        *,
        extra: dict[str, Any] | None = None,
        events: Deque[SampleLogEvent] | None = None,
        max_events: int = DEFAULT_BUFFER_MAX_EVENTS,
        time_provider=None,
    ) -> None:
        self._extra = dict(extra or {})
        self._max_events = max(1, int(max_events))
        self._events = events if events is not None else deque(maxlen=self._max_events)
        self._time_provider = time_provider or datetime.now
        self._closed = False
        self._options = (None, None, None, None, None, None, None, None, self._extra)

    @property
    def events(self) -> tuple[SampleLogEvent, ...]:
        return tuple(self._events)

    @property
    def event_count(self) -> int:
        return len(self._events)

    def bind(self, **extra: Any) -> BufferedSampleLogger:
        merged = dict(self._extra)
        merged.update(extra)
        return type(self)(
            extra=merged,
            events=self._events,
            max_events=self._max_events,
            time_provider=self._time_provider,
        )

    def close(self) -> None:
        self._closed = True

    def discard(self) -> None:
        self._events.clear()
        self._closed = True

    def replay_to(self, path: str) -> None:
        replay_sample_log_events(path, self.events, time_provider=self._time_provider)

    def log(self, level: Any, message: Any, *args: Any, **kwargs: Any) -> None:
        if self._closed:
            return
        module = str(self._extra.get("module", "No module"))
        raw = "raw" in self._extra
        rendered = SampleLogger._render_message(message, *args, **kwargs)
        normalized_level = SampleLogger._normalize_level(level)
        timestamp = self._time_provider().strftime("%m-%d %H:%M:%S.%f")[:-3]
        self._events.append(
            SampleLogEvent(
                timestamp=timestamp,
                module=module,
                level=normalized_level,
                message=rendered,
                raw=raw,
            )
        )
        _forward_sample_event(
            extra=self._extra,
            module=module,
            level=normalized_level,
            message=rendered,
            raw=raw,
        )

    def debug(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("DEBUG", message, *args, **kwargs)

    def info(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("INFO", message, *args, **kwargs)

    def warning(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("WARNING", message, *args, **kwargs)

    def error(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("ERROR", message, *args, **kwargs)

    def critical(self, message: Any, *args: Any, **kwargs: Any) -> None:
        self.log("CRITICAL", message, *args, **kwargs)
