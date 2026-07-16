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
    "scan_logs_dir": None,
    "component": None,
}

# Canonical files under ``logs/<scan>/`` (one process → one primary file).
COMPONENT_LOG_BASENAME = {
    "core": "core.log",
    "control": "core.log",
    "jarvis": "core.log",
    "factory": "factory.log",  # same process as core unless setup alone
    "sampler": "sampler.log",
    "archiver": "archiver.log",
    "worker": "worker.log",  # overridden by worker_id → worker-00.log
}

# V1 presentation labels (``·•· <module>``).
COMPONENT_MODULE_LABEL = {
    "core": "Jarvis-HEP",
    "control": "Jarvis-HEP",
    "jarvis": "Jarvis-HEP",
    "factory": "Factory",
    "sampler": "Sampler",
    "archiver": "Archiver",
    "worker": "Worker",
}


def scan_logs_dir(task_root: str, scan_name: str) -> str:
    """Return ``<task_root>/logs/<scan>/`` (absolute)."""
    root = os.path.abspath(os.path.expanduser(str(task_root or ".")))
    scan = str(scan_name or "scan").strip() or "scan"
    return os.path.join(root, "logs", scan)


def component_log_path(
    scan_logs: str,
    component: str,
    *,
    worker_id: int | None = None,
) -> str:
    """Resolve a component log file under the scan logs directory."""
    base = os.path.abspath(str(scan_logs))
    role = str(component or "jarvis").strip().lower() or "jarvis"
    if role == "worker" and worker_id is not None:
        return os.path.join(base, f"worker-{int(worker_id):02d}.log")
    name = COMPONENT_LOG_BASENAME.get(role, f"{role}.log")
    return os.path.join(base, name)


def component_module_label(component: str, *, worker_id: int | None = None) -> str:
    """Default ``·•·`` module label for a component."""
    role = str(component or "jarvis").strip().lower() or "jarvis"
    if role == "worker" and worker_id is not None:
        return f"Worker-{int(worker_id)}"
    return COMPONENT_MODULE_LABEL.get(role, role)


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
    _state["log_path"] = None
    _state["scan_logs_dir"] = None
    _state["component"] = None


def setup_jarvis_logging(
    *,
    level: str | int = "INFO",
    console: bool = True,
    role: str = "jarvis",
    component: str | None = None,
    log_dir: str = "logs",
    scan_logs_dir: str | None = None,
    log_path: str | None = None,
    worker_id: int | None = None,
    max_bytes: int = 5 * 2**20,
    backup_count: int = 7,
    use_queue: bool = True,
) -> str:
    """Configure process-level Jarvis logging once per process.

    Preferred layout (scan-scoped, by component)::

        logs/<scan>/core.log
        logs/<scan>/worker-00.log
        logs/<scan>/archiver.log

    Resolution order for the file path:

    1. Explicit ``log_path``
    2. ``scan_logs_dir`` + ``component``/``role`` (+ ``worker_id`` for workers)
    3. Legacy fallback: ``<log_dir>/jarvis_{role}_{pid}.log``
    """
    shutdown_jarvis_logging()

    resolved_level = _resolve_level(level)
    comp = str(component or role or "jarvis").strip().lower() or "jarvis"
    scan_dir = str(scan_logs_dir or "").strip() or None

    if log_path:
        resolved_path = os.path.abspath(str(log_path))
    elif scan_dir:
        resolved_path = component_log_path(scan_dir, comp, worker_id=worker_id)
    else:
        log_root = os.path.abspath(log_dir)
        os.makedirs(log_root, exist_ok=True)
        resolved_path = os.path.join(log_root, f"jarvis_{comp}_{os.getpid()}.log")

    parent = os.path.dirname(resolved_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

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
    _state["scan_logs_dir"] = os.path.dirname(resolved_path)
    _state["component"] = comp
    atexit.register(shutdown_jarvis_logging)
    return resolved_path


def get_jarvis_logger(
    name: str = "jarvis",
    *,
    module: str | None = None,
    worker_id: int | None = None,
) -> JarvisLoggerAdapter:
    """Return a bound-capable adapter over the Jarvis log domain.

    ``module`` sets the V1 ``·•·`` label (via ``jarvis_module``). When omitted,
    a sensible default is derived from ``name`` / ``worker_id``.
    """
    qualified = (
        name
        if name.startswith(JARVIS_HEP_LOG_DOMAIN)
        else f"{JARVIS_HEP_LOG_DOMAIN}.{name}"
    )
    # Strip domain for component heuristic.
    short = name
    if short.startswith(f"{JARVIS_HEP_LOG_DOMAIN}."):
        short = short[len(JARVIS_HEP_LOG_DOMAIN) + 1 :]
    root = short.split(".", 1)[0] if short else "jarvis"
    label = module
    if label is None:
        if root.startswith("sampler"):
            label = f"Sampler:{short.split('.', 1)[-1]}" if "." in short else "Sampler"
        else:
            label = component_module_label(root, worker_id=worker_id)
    return JarvisLoggerAdapter(
        logging.getLogger(qualified),
        {"jarvis_module": label},
    )


__all__ = [
    "COMPONENT_LOG_BASENAME",
    "COMPONENT_MODULE_LABEL",
    "JARVIS_HEP_LOG_DOMAIN",
    "JarvisContextFormatter",
    "JarvisLoggerAdapter",
    "component_log_path",
    "component_module_label",
    "format_record_context",
    "get_jarvis_logger",
    "scan_logs_dir",
    "setup_jarvis_logging",
    "shutdown_jarvis_logging",
]
