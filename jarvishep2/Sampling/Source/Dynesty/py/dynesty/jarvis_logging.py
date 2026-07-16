#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Jarvis logger bridge for the vendored dynesty package (D13.5).

All dynesty progress / warning / info output should go through
:func:`get_dynesty_logger` so formatting matches the V1 visual contract
(``·•· Dynesty.Inner`` via :mod:`jarvishep2.logging`).
"""

from __future__ import annotations

import warnings
from typing import Any

# Process-local logger injected by DynestySampler / tests.
_DYNESTY_LOGGER: Any | None = None


def set_dynesty_logger(logger: Any | None) -> None:
    """Install the Jarvis logger used by all dynesty internals."""
    global _DYNESTY_LOGGER
    _DYNESTY_LOGGER = logger


def get_dynesty_logger() -> Any:
    """Return the active logger (injected or a default Jarvis adapter)."""
    if _DYNESTY_LOGGER is not None:
        return _DYNESTY_LOGGER
    try:
        from jarvishep2.logging import get_jarvis_logger

        return get_jarvis_logger(
            "sampler.dynesty.inner",
            module="Dynesty.Inner",
        )
    except Exception:
        # Extremely defensive fallback for bare unit imports.
        import logging

        return logging.getLogger("jarvishep2.dynesty")


def bind_inner(logger: Any | None = None) -> Any:
    """Return a child logger labeled ``Dynesty.Inner`` (V1 presentation)."""
    base = logger if logger is not None else get_dynesty_logger()
    binder = getattr(base, "bind", None)
    if callable(binder):
        try:
            return binder(module="Dynesty.Inner", jarvis_module="Dynesty.Inner")
        except Exception:
            try:
                return binder(module="Dynesty.Inner")
            except Exception:
                pass
    return base


def emit_info(message: str, *, logger: Any | None = None) -> None:
    log = bind_inner(logger)
    try:
        log.info(str(message))
    except Exception:
        pass


def emit_warning(
    message: str,
    category: type = UserWarning,
    *,
    logger: Any | None = None,
) -> None:
    """Route a dynesty warning through Jarvis logger (V1 ``_emit_warning``)."""
    log = bind_inner(logger)
    text = str(message)
    if category is not UserWarning and category is not None:
        try:
            name = getattr(category, "__name__", str(category))
            text = f"{name}: {text}"
        except Exception:
            pass
    try:
        log.warning(text)
        return
    except Exception:
        pass
    try:
        warnings.warn(str(message), category or UserWarning, stacklevel=3)
    except Exception:
        pass


def emit_progress(
    print_func,
    logger_obj,
    *args,
    **kwargs,
) -> Any:
    """Call progress printer with optional ``logger=`` (V1 ``_emit_progress``)."""
    if logger_obj is None:
        logger_obj = get_dynesty_logger()
    try:
        return print_func(*args, logger=logger_obj, **kwargs)
    except TypeError:
        return print_func(*args, **kwargs)


def resolve_progress_title(logger: Any | None = None) -> str:
    """Title for progress blocks (``Dynesty`` / ``MultiNest``)."""
    title = "Dynesty"
    log = logger if logger is not None else get_dynesty_logger()
    module = None
    try:
        extra = getattr(log, "extra", None)
        if isinstance(extra, dict):
            module = str(
                extra.get("jarvis_module")
                or extra.get("module")
                or ""
            )
    except Exception:
        module = None
    if not module:
        try:
            # loguru-style bound logger
            options = getattr(log, "_options", None)
            if isinstance(options, tuple) and len(options) >= 9:
                extra = options[8]
                if isinstance(extra, dict):
                    module = str(extra.get("module", "") or extra.get("jarvis_module", ""))
        except Exception:
            module = None
    if module:
        if "MultiNest" in module:
            return "MultiNest"
        if "Dynesty" in module:
            return "Dynesty"
    return title


def format_progress_block(fn_args, rate=None, logger=None) -> str:
    """V1-style multi-line progress block for Jarvis logger."""
    metrics = getattr(fn_args, "metrics", None) or {}
    rows = [("iter", fn_args.niter)]
    if rate is not None:
        try:
            import numpy as np

            if np.isfinite(rate):
                rows.append(("rate(it/s)", "{:.2f}".format(float(rate))))
        except Exception:
            pass
    if metrics.get("add_live_it") is not None:
        rows.append(("add_live", "+{:d}".format(int(metrics["add_live_it"]))))
    if metrics.get("batch") is not None:
        rows.append(("batch", int(metrics["batch"])))
    if "bound" in metrics:
        rows.append(("bound", metrics["bound"]))
    if "nc" in metrics:
        rows.append(("nc", metrics["nc"]))
    if "ncall" in metrics:
        rows.append(("ncall", metrics["ncall"]))
    if "eff" in metrics:
        rows.append(("eff(%)", "{:6.3f}".format(float(metrics["eff"]))))
    if all(k in metrics for k in ("logl_min", "loglstar", "logl_max")):
        rows.append(
            (
                "loglstar",
                "{:6.3f} < {:6.3f} < {:6.3f}".format(
                    metrics["logl_min"],
                    metrics["loglstar"],
                    metrics["logl_max"],
                ),
            )
        )
    if "logz" in metrics and "logzerr" in metrics:
        rows.append(
            (
                "logz",
                "{:6.3f} +/- {:6.3f}".format(metrics["logz"], metrics["logzerr"]),
            )
        )
    if metrics.get("dlogz") is not None and "delta_logz" in metrics:
        rows.append(
            (
                "dlogz",
                "{:6.3f} > {:6.3f}".format(metrics["delta_logz"], metrics["dlogz"]),
            )
        )
    elif metrics.get("stop_val") is not None:
        rows.append(("stop", "{:6.3f}".format(metrics["stop_val"])))

    title = resolve_progress_title(logger)
    return "{} Progress ->\n".format(title) + "\n".join(
        "\t{:<10} -> {}".format(name, value) for name, value in rows
    )


def install_warnings_bridge() -> None:
    """Route Python ``warnings.warn`` from dynesty modules to Jarvis logger.

    Safe to call multiple times. Only intercepts messages whose stack includes
    ``dynesty`` paths so unrelated library warnings stay intact.
    """
    if getattr(install_warnings_bridge, "_installed", False):
        return

    original_showwarning = warnings.showwarning

    def _showwarning(message, category, filename, lineno, file=None, line=None):
        path = str(filename or "")
        if "dynesty" in path.replace("\\", "/").lower():
            emit_warning(str(message), category)
            return
        return original_showwarning(message, category, filename, lineno, file=file, line=line)

    warnings.showwarning = _showwarning  # type: ignore[assignment]
    install_warnings_bridge._installed = True  # type: ignore[attr-defined]


__all__ = [
    "bind_inner",
    "emit_info",
    "emit_progress",
    "emit_warning",
    "format_progress_block",
    "get_dynesty_logger",
    "install_warnings_bridge",
    "resolve_progress_title",
    "set_dynesty_logger",
]
