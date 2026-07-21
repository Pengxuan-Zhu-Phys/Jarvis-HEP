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
from typing import Any, Mapping

from jarvishep2.logging.style import process_style

JARVIS_HEP_LOG_DOMAIN = "jarvis_hep"

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
    "log_paths": None,
    "scan_logs_dir": None,
    "component": None,
}

# Canonical files under ``logs/<scan>/``.
COMPONENT_LOG_BASENAME = {
    "core": "core.log",
    "control": "core.log",
    "jarvis": "core.log",
    "factory": "factory.log",
    "sampler": "sampler.log",
    "archiver": "archiver.log",
    "datarecorder": "datarecorder.log",
    "worker": "worker.log",  # overridden by worker_id → worker-00.log
}

# Control-process multi-sink files (same OS process, separate RotatingFileHandlers).
CONTROL_MULTI_SINK_COMPONENTS: tuple[str, ...] = (
    "core",
    "factory",
    "sampler",
    "archiver",
    "datarecorder",
)

# Archiver process: pack lifecycle + DATABASE DataRecorder.
ARCHIVER_MULTI_SINK_COMPONENTS: tuple[str, ...] = (
    "archiver",
    "datarecorder",
)

# Process presentation labels (``·•· <module>``). Always ``Jarvis-HEP`` + dot segments.
# Examples: Jarvis-HEP, Jarvis-HEP.Factory, Jarvis-HEP.Sampler.Dynesty.Inner
JARVIS_HEP_LABEL_ROOT = "Jarvis-HEP"

COMPONENT_MODULE_LABEL = {
    "core": "Jarvis-HEP",
    "control": "Jarvis-HEP",
    "jarvis": "Jarvis-HEP",
    "factory": "Jarvis-HEP.Factory",
    "sampler": "Jarvis-HEP.Sampler",
    "archiver": "Jarvis-HEP.Archiver",
    "datarecorder": "Jarvis-HEP.DataRecorder",
    "worker": "Jarvis-HEP.Worker",
}

DATARECORDER_MODULE_LABEL = "Jarvis-HEP.DataRecorder"


# Canonical casing for well-known segments (logger short names → display).
_LABEL_SEGMENT_CANON: dict[str, str] = {
    "sampler": "Sampler",
    "factory": "Factory",
    "archiver": "Archiver",
    "worker": "Worker",
    "datarecorder": "DataRecorder",
    "dynesty": "Dynesty",
    "multinest": "MultiNest",
    "inner": "Inner",
    "pool": "Pool",
    "monitor": "Monitor",
    "watchdog": "Watchdog",
}


def _canon_label_segment(seg: str) -> str:
    key = str(seg or "").strip()
    if not key:
        return ""
    mapped = _LABEL_SEGMENT_CANON.get(key.lower())
    if mapped:
        return mapped
    # Digit worker ids stay zero-padded as given.
    if key.isdigit():
        return key
    # Title-case pure lower/upper tokens; keep MixedCase (e.g. MultiNest).
    if key.islower() or key.isupper():
        return key[:1].upper() + key[1:].lower() if len(key) > 1 else key.upper()
    return key


def jarvis_module_label(*parts: str) -> str:
    """Build a presentation label ``Jarvis-HEP[.Seg...]`` (dots only, no colons).

    Segments are cleaned of ``:`` / spaces; a leading ``Jarvis-HEP`` is not doubled.
    """
    segments: list[str] = []
    for raw in parts:
        text = str(raw or "").strip()
        if not text:
            continue
        text = text.replace(":", ".").replace(" ", ".")
        while ".." in text:
            text = text.replace("..", ".")
        text = text.strip(".")
        if not text:
            continue
        if text == JARVIS_HEP_LABEL_ROOT:
            continue
        if text.startswith(f"{JARVIS_HEP_LABEL_ROOT}."):
            text = text[len(JARVIS_HEP_LABEL_ROOT) + 1 :]
        for seg in text.split("."):
            seg = _canon_label_segment(seg)
            if not seg or seg == JARVIS_HEP_LABEL_ROOT:
                continue
            segments.append(seg)
    if not segments:
        return JARVIS_HEP_LABEL_ROOT
    return JARVIS_HEP_LABEL_ROOT + "." + ".".join(segments)


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
    """Default ``·•·`` module label for a component (``Jarvis-HEP.*``)."""
    role = str(component or "jarvis").strip().lower() or "jarvis"
    if role == "worker" and worker_id is not None:
        return jarvis_module_label("Worker", f"{int(worker_id):02d}")
    if role in COMPONENT_MODULE_LABEL:
        return COMPONENT_MODULE_LABEL[role]
    return jarvis_module_label(role)


def resolve_record_component(record: logging.LogRecord) -> str:
    """Map a LogRecord to a scan log component (``core`` / ``factory`` / ``sampler`` / …).

    Routing keys (first match wins):
    1. Logger name under ``jarvis_hep.<component>...``
    2. Presentation label ``jarvis_module`` (Dynesty.Inner → sampler, Factory → factory)
    3. Fallback ``core`` so nothing is dropped.
    """
    name = str(getattr(record, "name", "") or "")
    if name.startswith(f"{JARVIS_HEP_LOG_DOMAIN}."):
        short = name[len(JARVIS_HEP_LOG_DOMAIN) + 1 :]
    else:
        short = name
    root = short.split(".", 1)[0].strip().lower() if short else ""
    if root in {"factory", "sampler", "archiver", "worker", "datarecorder"}:
        return root
    if root in {"core", "control", "jarvis", "client"}:
        return "core"

    mod = str(record.__dict__.get("jarvis_module") or "").strip()
    if mod:
        # Accept both new Jarvis-HEP.* labels and legacy short names.
        if "DataRecorder" in mod or "hdf5" in mod.lower():
            return "datarecorder"
        if (
            ".Sampler" in mod
            or mod.startswith("Sampler")
            or "Dynesty" in mod
            or "MultiNest" in mod
            or mod.endswith(".Inner")
            or ".Pool" in mod
        ):
            return "sampler"
        if ".Factory" in mod or mod == "Factory" or mod.endswith(".Factory"):
            return "factory"
        if ".Archiver" in mod or mod == "Archiver" or mod.endswith(".Archiver"):
            return "archiver"
        if ".Worker" in mod or mod.startswith("Worker"):
            return "worker"
    return "core"


class ComponentFilter(logging.Filter):
    """Allow only records whose resolved component is in *accept*."""

    def __init__(self, accept: set[str] | frozenset[str] | tuple[str, ...]) -> None:
        super().__init__()
        self.accept = frozenset(str(item) for item in accept)

    def filter(self, record: logging.LogRecord) -> bool:
        return resolve_record_component(record) in self.accept


class CoreResidualFilter(logging.Filter):
    """Residual sink: accept *primary* plus any component without its own file.

    Used for ``core.log`` (control) and ``archiver.log`` (archiver process) so
    unclaimed components are never dropped.
    """

    def __init__(
        self,
        dedicated: set[str] | frozenset[str] | tuple[str, ...],
        *,
        primary: str = "core",
    ) -> None:
        super().__init__()
        self.primary = str(primary or "core")
        self.dedicated = frozenset(
            str(item) for item in dedicated if item != self.primary
        )

    def filter(self, record: logging.LogRecord) -> bool:
        resolved = resolve_record_component(record)
        return resolved == self.primary or resolved not in self.dedicated


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
    """Process-log formatter driven by ``card/logging.yaml`` (default = V1 look)."""

    def __init__(
        self,
        *,
        colorize: bool = False,
        style: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__()
        self.colorize = bool(colorize)
        self._style = dict(style or process_style())

    def formatTime(self, record: logging.LogRecord, datefmt: str | None = None) -> str:
        created = datetime.fromtimestamp(record.created)
        fmt = str(self._style.get("date_format") or "%m-%d %H:%M:%S")
        stamp = created.strftime(fmt)
        digits = int(self._style.get("millis_digits") or 0)
        if digits > 0:
            millis = int(record.msecs)
            stamp = f"{stamp}.{millis:03d}"
        return stamp

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

        # DataRecorder (DATABASE / samples.hdf5) uses the special Ϡ bullet.
        data_mod = str(
            self._style.get("hdf5_writer_module") or DATARECORDER_MODULE_LABEL
        )
        default_bullet = str(self._style.get("bullet") or "·•·")
        data_bullet = str(self._style.get("hdf5_writer_bullet") or "Ϡ")
        bullet = (
            data_bullet
            if module == data_mod or "DataRecorder" in module
            else default_bullet
        )
        head_tmpl = str(
            self._style.get("head")
            or "\n{bullet} {module} \n\t-> {timestamp} - [{level}] >>> \n"
        )

        if self.colorize and sys.stderr.isatty():
            cyan, green, reset = "\033[36m", "\033[32m", "\033[0m"
            level_colors = {
                "DEBUG": "\033[36m",
                "INFO": "\033[32m",
                "WARNING": "\033[33m",
                "ERROR": "\033[31m",
                "CRITICAL": "\033[1;31m",
            }
            level_color = level_colors.get(level, "")
            # Colorize module/timestamp/level inside the configured head template.
            head = head_tmpl.format(
                bullet=bullet,
                module=f"{cyan}{module}{reset}",
                timestamp=f"{green}{timestamp}{reset}",
                level=f"{level_color}{level}{reset}" if level_color else level,
                message="",
            )
            body = f"{level_color}{message}{reset}" if level_color else message
        else:
            head = head_tmpl.format(
                bullet=bullet,
                module=module,
                timestamp=timestamp,
                level=level,
                message="",
            )
            body = message

        if bool(self._style.get("append_context", True)):
            context = format_record_context(record)
            if context:
                body = f"{body} {context}" if body else context
        return f"{head}{body}"


class JarvisLoggerAdapter(logging.LoggerAdapter):
    """loguru-like adapter with .bind() over stdlib logging."""

    @staticmethod
    def _normalize_extra(extra: dict[str, Any]) -> dict[str, Any]:
        """Map ``module=`` onto ``jarvis_module`` as ``Jarvis-HEP.*`` (dots only)."""
        normalized = dict(extra)
        if "module" in normalized and "jarvis_module" not in normalized:
            raw = normalized.pop("module")
            text = str(raw or "").strip()
            if text == JARVIS_HEP_LABEL_ROOT:
                normalized["jarvis_module"] = JARVIS_HEP_LABEL_ROOT
            elif text:
                normalized["jarvis_module"] = jarvis_module_label(text)
            else:
                normalized["jarvis_module"] = JARVIS_HEP_LABEL_ROOT
        elif "module" in normalized:
            normalized.pop("module")
        # Normalize existing jarvis_module to the canonical stamp when possible.
        if "jarvis_module" in normalized and normalized["jarvis_module"] is not None:
            jm = str(normalized["jarvis_module"]).strip()
            if jm == JARVIS_HEP_LABEL_ROOT:
                normalized["jarvis_module"] = JARVIS_HEP_LABEL_ROOT
            elif jm and not jm.startswith(f"{JARVIS_HEP_LABEL_ROOT}."):
                # Leave sample labels (Sample@uuid...) alone.
                if jm.startswith("Sample@") or " (Likelihood)" in jm or " (Operas" in jm:
                    pass
                else:
                    normalized["jarvis_module"] = jarvis_module_label(jm)
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


def _make_console_handler(*, level: int, style: Mapping[str, Any] | None = None) -> logging.Handler:
    stream = logging.StreamHandler(sys.stderr)
    stream.setLevel(level)
    stream.setFormatter(JarvisContextFormatter(colorize=True, style=style))
    return stream


def _make_file_handler(
    *,
    log_path: str,
    level: int,
    max_bytes: int,
    backup_count: int,
    style: Mapping[str, Any] | None = None,
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
    handler.setFormatter(JarvisContextFormatter(colorize=False, style=style))
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
    _state["log_paths"] = None
    _state["scan_logs_dir"] = None
    _state["component"] = None


def setup_jarvis_logging(
    *,
    level: str | int = "INFO",
    console: bool = True,
    console_level: str | int | None = None,
    silence: bool = False,
    role: str = "jarvis",
    component: str | None = None,
    log_dir: str = "logs",
    scan_logs_dir: str | None = None,
    log_path: str | None = None,
    style_path: str | None = None,
    worker_id: int | None = None,
    max_bytes: int = 5 * 2**20,
    backup_count: int = 7,
    use_queue: bool = True,
    multi_sink: bool | None = None,
) -> str:
    """Configure process-level Jarvis logging once per process.

    Preferred layout (scan-scoped, by component)::

        logs/<scan>/core.log
        logs/<scan>/factory.log
        logs/<scan>/sampler.log
        logs/<scan>/archiver.log
        logs/<scan>/worker-00.log

    **Control process** (role/component ``core``) with ``scan_logs_dir`` opens a
    multi-sink by default: Factory / Sampler / Archiver / DataRecorder lines go
    to dedicated files; residual control lines stay in ``core.log``.
    **Archiver process** multi-sinks ``archiver.log`` + ``datarecorder.log``.
    Worker processes keep a single ``worker-NN.log``.

    Console:
      - default on at INFO (``console_level``)
      - ``silence=True`` disables screen output
      - file level follows ``level`` (default INFO)

    Style templates load from ``card/logging.yaml`` (or ``style_path``).
    Returns the primary log path (``core.log`` for multi-sink control).
    """
    from jarvishep2.logging.style import clear_logging_style_cache, process_style

    shutdown_jarvis_logging()
    clear_logging_style_cache()
    style = process_style(style_path)

    file_level = _resolve_level(level)
    # Console defaults to INFO; silence wins over console=True.
    use_console = bool(console) and not bool(silence)
    cons_level = _resolve_level(
        console_level if console_level is not None else level
    )
    # Root logger must admit the lowest handler threshold.
    root_level = min(file_level, cons_level) if use_console else file_level

    comp = str(component or role or "jarvis").strip().lower() or "jarvis"
    scan_dir = str(scan_logs_dir or "").strip() or None

    # Multi-sink when we have a scan logs directory (control or archiver process).
    # Explicit log_path forces single-file (tests / specialized callers).
    if multi_sink is None:
        multi_sink = bool(
            scan_dir
            and worker_id is None
            and comp in {"core", "control", "jarvis", "archiver"}
            and not log_path
        )
    else:
        multi_sink = bool(multi_sink)

    log_paths: dict[str, str] = {}
    residual_primary = "core"
    if multi_sink and scan_dir:
        os.makedirs(os.path.abspath(scan_dir), exist_ok=True)
        if comp == "archiver":
            sink_names = ARCHIVER_MULTI_SINK_COMPONENTS
            residual_primary = "archiver"
        else:
            sink_names = CONTROL_MULTI_SINK_COMPONENTS
            residual_primary = "core"
        for name in sink_names:
            log_paths[name] = component_log_path(scan_dir, name)
        resolved_path = log_paths.get(residual_primary) or next(iter(log_paths.values()))
    elif log_path:
        resolved_path = os.path.abspath(str(log_path))
        log_paths[comp if comp != "control" else "core"] = resolved_path
    elif scan_dir:
        resolved_path = component_log_path(scan_dir, comp, worker_id=worker_id)
        key = "worker" if comp == "worker" else (comp if comp != "control" else "core")
        log_paths[key] = resolved_path
    else:
        log_root = os.path.abspath(log_dir)
        os.makedirs(log_root, exist_ok=True)
        resolved_path = os.path.join(log_root, f"jarvis_{comp}_{os.getpid()}.log")
        log_paths[comp] = resolved_path

    for path in log_paths.values():
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)

    logger = logging.getLogger(JARVIS_HEP_LOG_DOMAIN)
    logger.handlers.clear()
    logger.setLevel(root_level)
    logger.propagate = False

    sink_handlers: list[logging.Handler] = []
    if use_console:
        sink_handlers.append(_make_console_handler(level=cons_level, style=style))

    if multi_sink and scan_dir:
        dedicated = set(log_paths)
        for name, path in log_paths.items():
            handler = _make_file_handler(
                log_path=path,
                level=file_level,
                max_bytes=max_bytes,
                backup_count=backup_count,
                style=style,
            )
            if name == residual_primary:
                # Residual bucket: primary + any component without a dedicated sink.
                handler.addFilter(
                    CoreResidualFilter(dedicated, primary=residual_primary)
                )
            else:
                handler.addFilter(ComponentFilter({name}))
            sink_handlers.append(handler)
    else:
        sink_handlers.append(
            _make_file_handler(
                log_path=resolved_path,
                level=file_level,
                max_bytes=max_bytes,
                backup_count=backup_count,
                style=style,
            )
        )

    if use_queue and sink_handlers:
        log_queue: queue.Queue[logging.LogRecord] = queue.Queue(-1)
        queue_handler = QueueHandler(log_queue)
        queue_handler.setLevel(root_level)
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
    _state["log_paths"] = dict(log_paths)
    _state["scan_logs_dir"] = (
        os.path.abspath(scan_dir)
        if scan_dir
        else os.path.dirname(resolved_path)
    )
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

    ``module`` sets the ``·•·`` presentation label (via ``jarvis_module``).
    Prefer ``Jarvis-HEP.*`` with dots only. When omitted, a label is derived
    from ``name`` / ``worker_id`` (e.g. ``sampler.dynesty`` →
    ``Jarvis-HEP.Sampler.Dynesty``).
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
    if module is not None:
        label = jarvis_module_label(module) if module else JARVIS_HEP_LABEL_ROOT
        # Keep exact "Jarvis-HEP" when caller passed the root stamp.
        if str(module).strip() == JARVIS_HEP_LABEL_ROOT:
            label = JARVIS_HEP_LABEL_ROOT
    elif root == "sampler" or short.startswith("sampler."):
        rest = short.split(".", 1)[1] if "." in short else ""
        label = (
            jarvis_module_label("Sampler", *rest.split("."))
            if rest
            else jarvis_module_label("Sampler")
        )
    elif root == "factory" or short.startswith("factory."):
        label = jarvis_module_label("Factory")
    else:
        label = component_module_label(root, worker_id=worker_id)
    return JarvisLoggerAdapter(
        logging.getLogger(qualified),
        {"jarvis_module": label},
    )


__all__ = [
    "ARCHIVER_MULTI_SINK_COMPONENTS",
    "COMPONENT_LOG_BASENAME",
    "COMPONENT_MODULE_LABEL",
    "CONTROL_MULTI_SINK_COMPONENTS",
    "ComponentFilter",
    "CoreResidualFilter",
    "DATARECORDER_MODULE_LABEL",
    "JARVIS_HEP_LABEL_ROOT",
    "JARVIS_HEP_LOG_DOMAIN",
    "JarvisContextFormatter",
    "JarvisLoggerAdapter",
    "component_log_path",
    "component_module_label",
    "format_record_context",
    "get_jarvis_logger",
    "jarvis_module_label",
    "resolve_record_component",
    "scan_logs_dir",
    "setup_jarvis_logging",
    "shutdown_jarvis_logging",
]
