"""Compatibility import aliases for the D25.10 package layout.

Canonical packages:

* ``jarvishep2.sampling`` (legacy ``jarvishep2.Sampling``)
* ``jarvishep2.calculators`` (legacy ``jarvishep2.Module``)
* ``jarvishep2.cardload.contracts`` (legacy ``jarvishep2.contracts``)

Single-file modules keep a shim at the old path (see ``jarvishep2/client.py``
etc.). This finder only covers the renamed *packages*, so
``from jarvishep2.Sampling.bridson import Bridson`` still resolves.
"""

from __future__ import annotations

import importlib
import sys
from importlib.abc import Loader, MetaPathFinder
from importlib.machinery import ModuleSpec
from types import ModuleType

_PACKAGE_ALIASES: tuple[tuple[str, str], ...] = (
    ("jarvishep2.Sampling", "jarvishep2.sampling"),
    ("jarvishep2.Module", "jarvishep2.calculators"),
    ("jarvishep2.contracts", "jarvishep2.cardload.contracts"),
)


def canonical_module_name(fullname: str) -> str | None:
    """Return the canonical module name for a legacy alias, or None."""
    for old, new in _PACKAGE_ALIASES:
        if fullname == old:
            return new
        prefix = old + "."
        if fullname.startswith(prefix):
            return new + fullname[len(old) :]
    return None


class _AliasLoader(Loader):
    """Reuse an already-executed canonical module object under a legacy name."""

    def __init__(self, module: ModuleType) -> None:
        self._module = module

    def create_module(self, spec: ModuleSpec) -> ModuleType:
        return self._module

    def exec_module(self, module: ModuleType) -> None:
        return None


class _LegacyPackageFinder(MetaPathFinder):
    """Map ``jarvishep2.Sampling.*`` (and friends) onto the canonical packages."""

    def find_spec(
        self,
        fullname: str,
        path: object | None,
        target: ModuleType | None = None,
    ) -> ModuleSpec | None:
        canonical = canonical_module_name(fullname)
        if canonical is None:
            return None
        module = importlib.import_module(canonical)
        spec = ModuleSpec(
            fullname,
            _AliasLoader(module),
            origin=getattr(module, "__file__", None),
            is_package=hasattr(module, "__path__"),
        )
        if hasattr(module, "__path__"):
            spec.submodule_search_locations = list(module.__path__)
        return spec


def install_legacy_paths() -> None:
    """Install the package alias finder once on ``sys.meta_path``."""
    if any(isinstance(item, _LegacyPackageFinder) for item in sys.meta_path):
        return
    sys.meta_path.insert(0, _LegacyPackageFinder())


__all__ = ["canonical_module_name", "install_legacy_paths"]
