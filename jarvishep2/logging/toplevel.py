#!/usr/bin/env python3
"""Top-level Jarvis logging over stdlib logging (no loguru)."""

from __future__ import annotations

import atexit
import logging
import os
from logging.handlers import QueueHandler, QueueListener, RotatingFileHandler
import queue
import sys
from typing import Any

JARVIS_HEP_LOG_DOMAIN = "jarvis_hep"
_DEFAULT_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s | %(message)s"
_DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

_RECORD_SKIP_KEYS = frozenset(
    {
        "args",
        "asctime",
        "created",
        "exc_info",
        "exc_text",
        "filename",
        "funcName",
        "levelname",
        "levelno",
        "lineno",
        "message",
        "module",
        "msecs",
        "msg",
        "name",
        "pathname",
        "process",
        "processName",
        "relativeCreated",
        "stack_info",
        "taskName",
        "thread",
        "threadName",
    }
)

_state: dict[str, Any] = {
    "configured": False,
    "listener": None,
    "log_queue": None,
    "log_path": None,
}


class JarvisContextFormatter(logging.Formatter):
    """Formatter that appends bound context as key=value pairs."""

    def format(self, record: logging.LogRecord) -> str:
        base = super().format(record)
        context = format_record_context(record)
        if context:
            return f"{base} {context}"
        return base


class JarvisLoggerAdapter(logging.LoggerAdapter):
    """loguru-like adapter with .bind() over stdlib logging."""

    def process(self, msg: Any, kwargs: dict[str, Any]) -> tuple[Any, dict[str, Any]]:
        extra = dict(kwargs.get("extra") or {})
        for key, value in self.extra.items():
            extra.setdefault(key, value)
        kwargs["extra"] = extra
        return msg, kwargs

    def bind(self, **ctx: Any) -> JarvisLoggerAdapter:
        merged = dict(self.extra)
        merged.update(ctx)
        return type(self)(self.logger, merged)


def format_record_context(record: logging.LogRecord) -> str:
    """Render non-standard record fields as sorted key=value context."""
    parts: list[str] = []
    for key in sorted(record.__dict__):
        if key in _RECORD_SKIP_KEYS or key.startswith("_"):
            continue
        value = record.__dict__[key]
        if value is None:
            continue
        parts.append(f"{key}={value}")
    return " ".join(parts)


def _resolve_level(level: str | int) -> int:
    if isinstance(level, int):
        return level
    name = str(level).strip().upper()
    return getattr(logging, name, logging.INFO)


def _make_console_handler(*, level: int) -> logging.Handler:
    stream = logging.StreamHandler(sys.stderr)
    stream.setLevel(level)
    try:
        import colorlog

        class _JarvisColoredFormatter(colorlog.ColoredFormatter):
            def format(self, record: logging.LogRecord) -> str:
                base = super().format(record)
                context = format_record_context(record)
                return f"{base} {context}" if context else base

        stream.setFormatter(
            _JarvisColoredFormatter(
                "%(log_color)s%(asctime)s %(levelname)s %(name)s | %(message)s%(reset)s",
                datefmt=_DEFAULT_DATE_FORMAT,
                log_colors={
                    "DEBUG": "cyan",
                    "INFO": "green",
                    "WARNING": "yellow",
                    "ERROR": "red",
                    "CRITICAL": "bold_red",
                },
            )
        )
    except ImportError:
        stream.setFormatter(
            JarvisContextFormatter(_DEFAULT_LOG_FORMAT, datefmt=_DEFAULT_DATE_FORMAT)
        )
    return stream


def _make_file_handler(
    *,
    log_path: str,
    level: int,
    max_bytes: int,
    backup_count: int,
) -> RotatingFileHandler:
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    handler = RotatingFileHandler(
        log_path,
        maxBytes=max_bytes,
        backupCount=backup_count,
        encoding="utf-8",
    )
    handler.setLevel(level)
    handler.setFormatter(JarvisContextFormatter(_DEFAULT_LOG_FORMAT, datefmt=_DEFAULT_DATE_FORMAT))
    return handler


def shutdown_jarvis_logging() -> None:
    """Stop the queue listener and drain pending records."""
    listener = _state.get("listener")
    if listener is not None:
        listener.stop()
    _state["listener"] = None
    _state["log_queue"] = None
    _state["configured"] = False


def setup_jarvis_logging(
    *,
    level: str | int = "INFO",
    console: bool = True,
    role: str = "jarvis",
    log_dir: str = "logs",
    max_bytes: int = 100 * 2**20,
    backup_count: int = 7,
    use_queue: bool = True,
) -> str:
    """Configure process-level Jarvis logging once per process."""
    shutdown_jarvis_logging()

    resolved_level = _resolve_level(level)

    log_root = os.path.abspath(log_dir)
    os.makedirs(log_root, exist_ok=True)
    log_path = os.path.join(log_root, f"jarvis_{role}_{os.getpid()}.log")

    logger = logging.getLogger(JARVIS_HEP_LOG_DOMAIN)
    logger.handlers.clear()
    logger.setLevel(resolved_level)
    logger.propagate = False

    sink_handlers: list[logging.Handler] = []
    if console:
        sink_handlers.append(_make_console_handler(level=resolved_level))
    sink_handlers.append(
        _make_file_handler(
            log_path=log_path,
            level=resolved_level,
            max_bytes=max_bytes,
            backup_count=backup_count,
        )
    )

    if use_queue and sink_handlers:
        log_queue: queue.Queue[logging.LogRecord] = queue.Queue(-1)
        queue_handler = QueueHandler(log_queue)
        queue_handler.setLevel(resolved_level)
        logger.addHandler(queue_handler)

        listener = QueueListener(log_queue, *sink_handlers, respect_handler_level=True)
        listener.start()
        _state["listener"] = listener
        _state["log_queue"] = log_queue
    else:
        for handler in sink_handlers:
            logger.addHandler(handler)

    _state["configured"] = True
    _state["log_path"] = log_path
    atexit.register(shutdown_jarvis_logging)
    return log_path


def get_jarvis_logger(name: str = "jarvis") -> JarvisLoggerAdapter:
    """Return a bound-capable adapter over the Jarvis log domain."""
    qualified = name if name.startswith(JARVIS_HEP_LOG_DOMAIN) else f"{JARVIS_HEP_LOG_DOMAIN}.{name}"
    return JarvisLoggerAdapter(logging.getLogger(qualified), {})