#!/usr/bin/env python3
"""Per-Sample logging primitives for Jarvis-HEP V2 (V1 contract, reimplemented)."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from datetime import datetime
import os
import threading
from typing import Any, Deque, Iterable


DEFAULT_BUFFER_MAX_EVENTS = 2048


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
        return self._time_provider().strftime("%m-%d %H:%M:%S.%f")[:-3]

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

        text = str(message)
        if not raw:
            prefix = "\n"
            if self._has_written and not self._last_ended_newline:
                prefix = "\n\n"
            ts = timestamp if timestamp is not None else self._timestamp()
            text = (
                f"{prefix}·•· {module} \n"
                f"\t-> {ts} - [{level}] >>> \n"
                f"{text}"
            )
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
    """Custom sample-local logging adapter used only for per-sample file logging."""

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
        self._backend.write(
            module=module,
            level=self._normalize_level(level),
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
        timestamp = self._time_provider().strftime("%m-%d %H:%M:%S.%f")[:-3]
        self._events.append(
            SampleLogEvent(
                timestamp=timestamp,
                module=module,
                level=SampleLogger._normalize_level(level),
                message=rendered,
                raw=raw,
            )
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