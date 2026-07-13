#!/usr/bin/env python3
"""Top-level Jarvis logging over stdlib logging (no loguru).

QueueListener architecture is preserved (D0.3). Rendering matches the V1 visual
contract (D12.2): ``·•·`` / ``Ϡ`` prefixes, module label, ``MM-DD HH:mm:ss.SSS``,
and raw message passthrough.
"""

from __future__ import annotations

import atexit
import logging
import os
import queue
import sys
from datetime import datetime
from logging.handlers import QueueHandler, QueueListener, RotatingFileHandler
from typing import Any

JARVIS_HEP_LOG_DOMAIN = "jarvis_hep"
_DEFAULT_DATE_FORMAT = "%m-%d %H:%M:%S"

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
        # Bound context / V1 presentation keys (not dumped as trailing kv).
        "Jarvis",
        "to_console",
        "_log_domain",
        "raw",
        "colorize",
        "jarvis_module",
    }
)

_state: dict[str, Any] = {
    "configured": False,
    "listener": None,
    "log_queue": None,
    "log_path": None,
}


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


class JarvisContextFormatter(logging.Formatter):
    """V1-style top-level formatter with optional trailing key=value context."""

    def __init__(self, *, colorize: bool = False) -> None:
        super().__init__()
        self.colorize = bool(colorize)

    def formatTime(self, record: logging.LogRecord, datefmt: str | None = None) -> str:
        created = datetime.fromtimestamp(record.created)
        stamp = created.strftime(_DEFAULT_DATE_FORMAT)
        millis = int(record.msecs)
        return f"{stamp}.{millis:03d}"

    def _module_label(self, record: logging.LogRecord) -> str:
        # ``module`` is a reserved LogRecord attribute (source file stem). The
        # V1 presentation label is carried as ``jarvis_module`` instead.
        bound = record.__dict__.get("jarvis_module")
        if isinstance(bound, str) and bound.strip():
            return bound.strip()
        name = str(record.name or "")
        if name.startswith(f"{JARVIS_HEP_LOG_DOMAIN}."):
            return name[len(JARVIS_HEP_LOG_DOMAIN) + 1 :] or name
        return name or "Jarvis-HEP"

    def format(self, record: logging.LogRecord) -> str:
        if bool(getattr(record, "raw", False)):
            message = record.getMessage()
            if record.exc_info:
                if not record.exc_text:
                    record.exc_text = self.formatException(record.exc_info)
                if record.exc_text:
                    message = f"{message}\n{record.exc_text}" if message else record.exc_text
            return message if message.endswith("\n") else f"{message}\n"

        module = self._module_label(record)
        timestamp = self.formatTime(record)
        level = record.levelname
        message = record.getMessage()
        if record.exc_info:
            if not record.exc_text:
                record.exc_text = self.formatException(record.exc_info)
            if record.exc_text:
                message = f"{message}\n{record.exc_text}" if message else record.exc_text

        bullet = "Ϡ" if module == "Jarvis-HEP.hdf5-Writter" else "·•·"
        if self.colorize and sys.stderr.isatty():
            # Approximate V1 loguru colors without loguru.
            cyan, green, reset = "\033[36m", "\033[32m", "\033[0m"
            level_colors = {
                "DEBUG": "\033[36m",
                "INFO": "\033[32m",
                "WARNING": "\033[33m",
                "ERROR": "\033[31m",
                "CRITICAL": "\033[1;31m",
            }
            level_color = level_colors.get(level, "")
            head = (
                f"\n{bullet} {cyan}{module}{reset} \n"
                f"\t-> {green}{timestamp}{reset} - "
                f"[{level_color}{level}{reset}] >>> \n"
            )
            body = f"{level_color}{message}{reset}" if level_color else message
        else:
            head = f"\n{bullet} {module} \n\t-> {timestamp} - [{level}] >>> \n"
            body = message

        context = format_record_context(record)
        if context:
            body = f"{body} {context}" if body else context
        return f"{head}{body}"


class JarvisLoggerAdapter(logging.LoggerAdapter):
    """loguru-like adapter with .bind() over stdlib logging."""

    @staticmethod
    def _normalize_extra(extra: dict[str, Any]) -> dict[str, Any]:
        """Map V1 ``module=`` presentation labels onto non-reserved ``jarvis_module``."""
        normalized = dict(extra)
        if "module" in normalized and "jarvis_module" not in normalized:
            normalized["jarvis_module"] = normalized.pop("module")
        elif "module" in normalized:
            normalized.pop("module")
        return normalized

    def process(self, msg: Any, kwargs: dict[str, Any]) -> tuple[Any, dict[str, Any]]:
        extra = self._normalize_extra(dict(kwargs.get("extra") or {}))
        for key, value in self.extra.items():
            extra.setdefault(key, value)
        kwargs["extra"] = self._normalize_extra(extra)
        return msg, kwargs

    def bind(self, **ctx: Any) -> JarvisLoggerAdapter:
        merged = dict(self.extra)
        merged.update(self._normalize_extra(ctx))
        return type(self)(self.logger, merged)


def _resolve_level(level: str | int) -> int:
    if isinstance(level, int):
        return level
    name = str(level).strip().upper()
    return getattr(logging, name, logging.INFO)


def _make_console_handler(*, level: int) -> logging.Handler:
    stream = logging.StreamHandler(sys.stderr)
    stream.setLevel(level)
    stream.setFormatter(JarvisContextFormatter(colorize=True))
    return stream


def _make_file_handler(
    *,
    log_path: str,
    level: int,
    max_bytes: int,
    backup_count: int,
) -> RotatingFileHandler:
    parent = os.path.dirname(log_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    handler = RotatingFileHandler(
        log_path,
        maxBytes=max_bytes,
        backupCount=backup_count,
        encoding="utf-8",
    )
    handler.setLevel(level)
    handler.setFormatter(JarvisContextFormatter(colorize=False))
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
    log_path: str | None = None,
    max_bytes: int = 5 * 2**20,
    backup_count: int = 7,
    use_queue: bool = True,
) -> str:
    """Configure process-level Jarvis logging once per process.

    File layout (D12.2):
    - If ``log_path`` is given, use it (control process prefers
      ``<task_root>/logs/<scan>/<scan>.log``).
    - Otherwise ``<log_dir>/jarvis_{role}_{pid}.log`` (Workers / ad-hoc).
    """
    shutdown_jarvis_logging()

    resolved_level = _resolve_level(level)

    if log_path:
        resolved_path = os.path.abspath(str(log_path))
        parent = os.path.dirname(resolved_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
    else:
        log_root = os.path.abspath(log_dir)
        os.makedirs(log_root, exist_ok=True)
        resolved_path = os.path.join(log_root, f"jarvis_{role}_{os.getpid()}.log")

    logger = logging.getLogger(JARVIS_HEP_LOG_DOMAIN)
    logger.handlers.clear()
    logger.setLevel(resolved_level)
    logger.propagate = False

    sink_handlers: list[logging.Handler] = []
    if console:
        sink_handlers.append(_make_console_handler(level=resolved_level))
    sink_handlers.append(
        _make_file_handler(
            log_path=resolved_path,
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
    _state["log_path"] = resolved_path
    atexit.register(shutdown_jarvis_logging)
    return resolved_path


def get_jarvis_logger(name: str = "jarvis") -> JarvisLoggerAdapter:
    """Return a bound-capable adapter over the Jarvis log domain."""
    qualified = (
        name
        if name.startswith(JARVIS_HEP_LOG_DOMAIN)
        else f"{JARVIS_HEP_LOG_DOMAIN}.{name}"
    )
    return JarvisLoggerAdapter(
        logging.getLogger(qualified),
        {"jarvis_module": name},
    )


__all__ = [
    "JARVIS_HEP_LOG_DOMAIN",
    "JarvisContextFormatter",
    "JarvisLoggerAdapter",
    "format_record_context",
    "get_jarvis_logger",
    "setup_jarvis_logging",
    "shutdown_jarvis_logging",
]
