# JarvisPLOT.layers registry
from __future__ import annotations
from typing import Dict, Type

registry: Dict[str, Type] = {}
LAYER_REGISTRY = registry  # public alias expected by Figure


def get(name: str):
    return registry.get(name)


def register_layer(cls: Type) -> None:
    name = getattr(cls, "TYPE", None)
    if not name:
        return
    registry[name] = cls


def _safe_import(module_path: str, class_name: str):
    try:
        mod = __import__(module_path, fromlist=[class_name])
        return getattr(mod, class_name, None)
    except Exception:
        return None

# ---- Attempt to import and register built-ins ----
# scatter (baseline)
_cls = _safe_import("JarvisPLOT.layers.scatter", "ScatterLayer")
if _cls:
    register_layer(_cls)

# triangle_scatter (ternary plot)
_cls = _safe_import("JarvisPLOT.layers.triangle_scatter", "TriangleScatterLayer")
if _cls:
    register_layer(_cls)

# colorbar (optional)
_cls = _safe_import("JarvisPLOT.layers.colorbar", "ColorbarLayer")
if _cls:
    register_layer(_cls)

# contour (optional)
_cls = _safe_import("JarvisPLOT.layers.contour", "ContourLayer")
if _cls:
    register_layer(_cls)

# other optional layers if present (no error if missing)
for mod, clsname in [
    ("JarvisPLOT.layers.line", "LineLayer"),
    ("JarvisPLOT.layers.band", "BandLayer"),
    ("JarvisPLOT.layers.ring_guides", "RingGuidesLayer"),
]:
    _cls = _safe_import(mod, clsname)
    if _cls:
        register_layer(_cls)

__all__ = ["registry", "LAYER_REGISTRY", "register_layer", "get"]