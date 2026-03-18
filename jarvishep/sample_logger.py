#!/usr/bin/env python3
from __future__ import annotations

from datetime import datetime
import os
import threading
from typing import Any


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

    def write(self, *, module: str, level: str, message: str, raw: bool) -> None:
        if self._closed:
            return

        text = str(message)
        if not raw:
            prefix = "\n"
            if self._has_written and not self._last_ended_newline:
                prefix = "\n\n"
            text = (
                f"{prefix}·•· {module} \n"
                f"\t-> {self._timestamp()} - [{level}] >>> \n"
                f"{text}"
            )

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


class SampleLogger:
    """Thin compatibility logger used only for per-sample file logging."""

    def __init__(self, backend: _SampleLogBackend, *, extra: dict[str, Any] | None = None) -> None:
        self._backend = backend
        self._extra = dict(extra or {})
        # Keep a minimal loguru-like shape for code that introspects bound extras.
        self._options = (None, None, None, None, None, None, None, None, self._extra)

    @classmethod
    def open(
        cls,
        path: str,
        *,
        module: str,
        time_provider=None,
        extra: dict[str, Any] | None = None,
    ) -> "SampleLogger":
        merged = dict(extra or {})
        merged.setdefault("module", module)
        backend = _SampleLogBackend(path, time_provider=time_provider)
        return cls(backend, extra=merged)

    def bind(self, **extra: Any) -> "SampleLogger":
        merged = dict(self._extra)
        merged.update(extra)
        return type(self)(self._backend, extra=merged)

    def close(self) -> None:
        self._backend.close()

    def log(self, level: Any, message: Any, *args: Any, **kwargs: Any) -> None:
        if "module" not in self._extra:
            module = "No module"
        else:
            module = str(self._extra.get("module"))

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
        if text:
            return text
        return "INFO"

    @staticmethod
    def _render_message(message: Any, *args: Any, **kwargs: Any) -> str:
        text = str(message)
        if args or kwargs:
            return text.format(*args, **kwargs)
        return text
