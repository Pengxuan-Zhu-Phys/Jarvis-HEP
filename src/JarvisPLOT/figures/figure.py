# JarvisPLOT Figure module
from __future__ import annotations
import matplotlib.pyplot as plt
from dataclasses import dataclass
import os, json
from typing import Dict, Any, List, Optional
from ..utils import deep_merge
from ..layers import registry as LAYER_REGISTRY

# ------------------ Spec dataclasses ------------------

@dataclass
class AxesSpec:
    name: str
    preset: Optional[str] = None
    cfg: Optional[Dict[str, Any]] = None

@dataclass
class LayerSpec:
    name: Optional[str]
    type: str
    axes: str
    cfg: Dict[str, Any]

@dataclass
class FigureSpec:
    name: str
    size: Optional[List[float]]
    axes: List[AxesSpec]
    layers: List[LayerSpec]

# ------------------ Runtime classes ------------------

@dataclass
class AxesRuntime:
    name: str
    ax: plt.Axes

class Figure:
    def __init__(self, mpl_fig, axes_map: Dict[str, AxesRuntime], spec: FigureSpec):
        self.fig = mpl_fig
        self.axes_map = axes_map
        self.spec = spec
        self.mappables = {}

    @classmethod
    def build(cls, spec: FigureSpec, style):
        # Prefer a resolved profile (attached by Engine) for figure defaults & axes
        resolved_profile = getattr(spec, "_resolved_profile", None)

        # Collect layer types to select per-layer defaults from the style bundle
        layer_types = [ls.type for ls in spec.layers]

        # ---- Figure size/dpi selection order ----
        # 1) profile.figure
        # 2) style.figure_defaults_for(layer_types)
        # 3) style.figure_defaults
        if resolved_profile and isinstance(resolved_profile, dict):
            pf_fig = resolved_profile.get("figure", {}) or {}
        else:
            pf_fig = {}

        fd = getattr(style, "figure_defaults_for", None)
        if callable(fd):
            fdef = style.figure_defaults_for(layer_types)
        else:
            fdef = getattr(style, "figure_defaults", {}) or {}

        size = spec.size or pf_fig.get("size") or fdef.get("size", [6, 4.5])
        dpi  = pf_fig.get("dpi", fdef.get("dpi", 150))
        fig = plt.figure(figsize=size, dpi=dpi)

        # ---- Axes creation ----
        axes_map: Dict[str, AxesRuntime] = {}

        if (not spec.axes) and resolved_profile and isinstance(resolved_profile, dict):
            # Use style helper to expand profile-defined axes (preset + override merged)
            expanded = style.expand_profile_axes(resolved_profile) if hasattr(style, "expand_profile_axes") else []
            axes_specs = [AxesSpec(name=ax.get("name", "main"), preset=None, cfg=ax) for ax in expanded]
        else:
            axes_specs = spec.axes

        for ax_spec in axes_specs:
            # base config from preset file (if specified) else from style presets
            if ax_spec.preset:
                # Load preset JSON relative to the style bundle root
                preset_path = os.path.join(style.root, ax_spec.preset)
                try:
                    with open(preset_path, "r", encoding="utf-8") as f:
                        base = json.load(f)
                except Exception:
                    base = {}
            else:
                apf = getattr(style, "axes_preset_for", None)
                if callable(apf):
                    base = style.axes_preset_for(ax_spec.name, layer_types)
                else:
                    base = style.axes_preset(ax_spec.name)

            merged = deep_merge(base, ax_spec.cfg or {})
            ax = _create_axes(fig, merged)
            axes_map[ax_spec.name] = AxesRuntime(ax_spec.name, ax)

        return cls(fig, axes_map, spec)
      
    def render(self, ctx):
        # Expose this figure's mappables to layers via context
        setattr(ctx, "_mappables", self.mappables)
        for layer_spec in self.spec.layers:
            layer_cls = LAYER_REGISTRY.get(layer_spec.type)
            if not layer_cls:
                available = ", ".join(sorted(LAYER_REGISTRY.keys()))
                raise RuntimeError(f"Unknown layer type '{layer_spec.type}'. Available: {available}")
            layer = layer_cls(layer_spec)

            # Resolve axes: if the requested axes isn't found, fall back to the first available axes
            ax_name = getattr(layer_spec, "axes", None)
            if ax_name not in self.axes_map:
                if self.axes_map:
                    fallback_name = next(iter(self.axes_map))
                    print(f"[JarvisPLOT] Warning: axes '{ax_name}' not found; using '{fallback_name}' instead.")
                    ax = self.axes_map[fallback_name].ax
                else:
                    raise RuntimeError("No axes available in this figure to render the layer.")
            else:
                ax = self.axes_map[ax_name].ax

            mappable = layer.render(ax, ctx)
            if layer_spec.name and mappable is not None:
                self.mappables[layer_spec.name] = mappable

    def save(self, out_spec):
        import os
        os.makedirs(out_spec.dir, exist_ok=True)
        paths = []
        for fmt in out_spec.formats:
            p = os.path.join(out_spec.dir, f"{self.spec.name}.{fmt}")
            self.fig.savefig(p, dpi=out_spec.dpi)
            paths.append(p)
        return paths

    def close(self):
        plt.close(self.fig)

# ------------------ Helper ------------------

def _create_axes(fig, cfg):
    rect = cfg.get("rect", [0.12,0.12,0.72,0.76])
    proj = cfg.get("projection", "cartesian")
    if proj == "polar":
        ax = fig.add_axes(rect, projection="polar")
    else:
        ax = fig.add_axes(rect)
    if cfg.get("grid"):
        ax.grid(True, which="both", alpha=0.3)
    if fc := cfg.get("facecolor"):
        ax.set_facecolor(fc)
    return ax