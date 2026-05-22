#!/usr/bin/env python3
from __future__ import annotations

from importlib import import_module
from typing import Any


Direction = str


class IOAdapter:
    """Jarvis-HEP calculator IO handler factory."""

    format_name: str
    direction: Direction

    def create_handler(self, **kwargs):
        raise NotImplementedError("IO adapters must create a runtime handler.")


class LegacyHandlerAdapter(IOAdapter):
    """Adapter wrapper around legacy Jarvis-HEP InputFile/OutputFile classes."""

    def __init__(
        self,
        format_name: str,
        direction: Direction,
        handler_module: str,
        handler_class: str,
        log_label: str,
    ):
        self.format_name = str(format_name)
        self.direction = _normalize_direction(direction)
        self.handler_module = str(handler_module)
        self.handler_class = str(handler_class)
        self.log_label = str(log_label)

    def _handler_type(self):
        module = import_module(self.handler_module)
        return getattr(module, self.handler_class)

    def create_handler(self, **kwargs):
        handler_type = self._handler_type()
        handler = handler_type(
            kwargs["name"],
            kwargs["path"],
            kwargs["file_type"],
            kwargs["variables"],
            kwargs["save"],
            kwargs["logger"],
            kwargs["PackID"],
            kwargs["sample_uuid"],
            kwargs["sample_save_dir"],
            kwargs["module"],
            kwargs["funcs"],
            kwargs.get("io_manager"),
        )
        _apply_handler_metadata(handler, kwargs)
        return handler


class PortalHandlerAdapter(IOAdapter):
    """Factory for Portal-backed calculator IO handlers.

    Format ownership and dispatch stay in Jarvis-Portal. This HEP-side adapter
    only records the resolved Portal adapter so unsupported-format errors and
    available-format listings come from the Portal registry.
    """

    def __init__(self, portal_adapter: Any, direction: Direction):
        self.portal_adapter = portal_adapter
        self.format_name = str(portal_adapter.format_name)
        self.direction = _normalize_direction(direction)
        self.log_label = f"Portal:{self.format_name}"

    def create_handler(self, **kwargs):
        if self.direction == "input":
            from jarvishep.IOs.portal import PortalInputFile as handler_type
        elif self.direction == "output":
            from jarvishep.IOs.portal import PortalOutputFile as handler_type
        else:
            raise ValueError(f"Unsupported Portal handler direction '{self.direction}'.")

        handler = handler_type(
            kwargs["name"],
            kwargs["path"],
            kwargs["file_type"],
            kwargs["variables"],
            kwargs["save"],
            kwargs["logger"],
            kwargs["PackID"],
            kwargs["sample_uuid"],
            kwargs["sample_save_dir"],
            kwargs["module"],
            kwargs["funcs"],
            kwargs.get("io_manager"),
        )
        handler.portal_adapter = self.portal_adapter
        _apply_handler_metadata(handler, kwargs)
        return handler


class IORegistry:
    def __init__(self):
        self._legacy_adapters: dict[tuple[str, Direction], IOAdapter] = {}

    def register(
        self,
        format_name: str,
        adapter: IOAdapter,
        direction: Direction | None = None,
        *,
        override: bool = False,
    ) -> None:
        fmt = _normalize_format(format_name)
        adapter_direction = _normalize_direction(direction or adapter.direction)
        key = (fmt, adapter_direction)
        if key in self._legacy_adapters and not override:
            raise ValueError(
                f"Calculator IO adapter for format '{format_name}' "
                f"and direction '{adapter_direction}' is already registered."
            )
        self._legacy_adapters[key] = adapter

    def get(self, format_name: str, direction: Direction) -> IOAdapter:
        fmt = _normalize_format(format_name)
        wanted_direction = _normalize_direction(direction)
        adapter = self._legacy_adapters.get((fmt, wanted_direction))
        if adapter is not None:
            return adapter

        both_adapter = self._legacy_adapters.get((fmt, "both"))
        if both_adapter is not None:
            return both_adapter

        portal_adapter = _portal_get(format_name, wanted_direction)
        return PortalHandlerAdapter(portal_adapter, wanted_direction)

    def available_formats(self, direction: Direction | None = None) -> list[str]:
        normalized_direction = _normalize_direction(direction) if direction is not None else None
        names = set(_portal_available_formats(normalized_direction))
        for (_, adapter_direction), adapter in self._legacy_adapters.items():
            if normalized_direction is None or adapter_direction in {normalized_direction, "both"}:
                names.add(adapter.format_name)
        return sorted(names, key=str.lower)


def _apply_handler_metadata(handler, kwargs) -> None:
    for key in ("header", "columns", "comment", "spec"):
        if key in kwargs and kwargs[key] is not None:
            setattr(handler, key, kwargs[key])


def _portal_get(format_name: str, direction: Direction):
    try:
        from jarvis_portal import get as portal_get
    except ImportError as exc:
        raise ImportError(
            "Portal-backed calculator IO requires Jarvis-HEP-Portal. "
            "Install it with `pip install Jarvis-HEP-Portal`."
        ) from exc
    return portal_get(format_name, direction)


def _portal_available_formats(direction: Direction | None = None) -> list[str]:
    try:
        from jarvis_portal import available_formats as portal_available_formats
    except ImportError:
        return []
    return portal_available_formats(direction)


def _normalize_format(format_name: str) -> str:
    text = str(format_name or "").strip()
    if not text:
        raise ValueError("Calculator IO format name must not be empty.")
    return text.casefold()


def _normalize_direction(direction: Direction) -> Direction:
    text = str(direction or "").strip().casefold()
    if text not in {"input", "output", "both"}:
        raise ValueError(f"Unsupported calculator IO direction '{direction}'.")
    return text


_DEFAULT_REGISTRY = IORegistry()
_BUILTINS_REGISTERED = False


def register(
    format_name: str,
    adapter: IOAdapter,
    direction: Direction | None = None,
    *,
    override: bool = False,
) -> None:
    _ensure_legacy_adapters()
    _DEFAULT_REGISTRY.register(format_name, adapter, direction, override=override)


def get(format_name: str, direction: Direction) -> IOAdapter:
    _ensure_legacy_adapters()
    return _DEFAULT_REGISTRY.get(format_name, direction)


def available_formats(direction: Direction | None = None) -> list[str]:
    _ensure_legacy_adapters()
    return _DEFAULT_REGISTRY.available_formats(direction)


def _ensure_legacy_adapters() -> None:
    global _BUILTINS_REGISTERED
    if _BUILTINS_REGISTERED:
        return

    builtins = [
        LegacyHandlerAdapter(
            "SLHA",
            "input",
            "jarvishep.IOs.Input",
            "SLHAInputFile",
            "SLHAInputFile",
        ),
        LegacyHandlerAdapter(
            "SLHA",
            "output",
            "jarvishep.IOs.Output",
            "SLHAOutputFile",
            "SLHAOutputFile",
        ),
        LegacyHandlerAdapter(
            "xSLHA",
            "output",
            "jarvishep.IOs.Output",
            "xSLHAOutputFile",
            "xSLHAOutputFile",
        ),
        LegacyHandlerAdapter(
            "File",
            "output",
            "jarvishep.IOs.Output",
            "FileOutput",
            "FileOutput",
        ),
    ]
    for adapter in builtins:
        _DEFAULT_REGISTRY.register(adapter.format_name, adapter, adapter.direction)
    _BUILTINS_REGISTERED = True
