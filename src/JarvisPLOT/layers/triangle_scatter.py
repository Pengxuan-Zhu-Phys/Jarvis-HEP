#!/usr/bin/env python3

# File: /Users/p.zhu/Workshop/Jarvis-HEP/src/JarvisPLOT/layers/triangle_scatter.py
from __future__ import annotations
from typing import Tuple, Optional
import numpy as np
from .base import LayerBase

# --- helpers -----------------------------------------------------------------

def _barycentric_to_xy(a: np.ndarray, b: np.ndarray, c: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Map barycentric coords (a,b,c), a+b+c=1, to 2D coords in an equilateral triangle.
    If c is None, use c = 1 - a - b.
    Triangle vertices: V1=(0,0), V2=(1,0), V3=(0.5, sqrt(3)/2).
    """
    if c is None:
        c = 1.0 - a - b
    # Vertex coordinates
    x1, y1 = 0.0, 0.0
    x2, y2 = 1.0, 0.0
    x3, y3 = 0.5, np.sqrt(3.0)/2.0
    x = a * x1 + b * x2 + c * x3
    y = a * y1 + b * y2 + c * y3
    return x, y


def _draw_triangle_frame(ax, linewidth=1.0, color="#444", alpha=1.0):
    import matplotlib.pyplot as plt
    verts = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3.0)/2.0], [0.0, 0.0]])
    ax.plot(verts[:,0], verts[:,1], linewidth=linewidth, color=color, alpha=alpha)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, np.sqrt(3.0)/2.0 + 0.05)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])


def _draw_grid(ax, ticks=(0.2, 0.4, 0.6, 0.8), color="#bbb", linewidth=0.5, alpha=0.6):
    # Grid lines parallel to each edge for reference
    h = np.sqrt(3.0)/2.0
    for t in ticks:
        # Lines parallel to base (V1-V2): y = h*(1 - t)
        y = h * (1.0 - t)
        x0 = (1.0 - t) * 0.5
        x1 = 1.0 - x0
        ax.plot([x0, x1], [y, y], color=color, linewidth=linewidth, alpha=alpha)
        # Lines parallel to edge (V2-V3)
        # paramization via two points at given t-level of a and b; draw using barycentric grid idea
        # a = t → line between points where a=t and b ranges in [0, 1-t]
        a = t
        b0, b1 = 0.0, 1.0 - a
        xA0, yA0 = _barycentric_to_xy(np.full(1,a), np.full(1,b0))
        xA1, yA1 = _barycentric_to_xy(np.full(1,a), np.full(1,b1))
        ax.plot([xA0[0], xA1[0]], [yA0[0], yA1[0]], color=color, linewidth=linewidth, alpha=alpha)
        # b = t → line between points where b=t and a ranges in [0, 1-t]
        b = t
        a0, a1 = 0.0, 1.0 - b
        xB0, yB0 = _barycentric_to_xy(np.full(1,a0), np.full(1,b))
        xB1, yB1 = _barycentric_to_xy(np.full(1,a1), np.full(1,b))
        ax.plot([xB0[0], xB1[0]], [yB0[0], yB1[0]], color=color, linewidth=linewidth, alpha=alpha)

# --- layer --------------------------------------------------------------------

class TriangleScatterLayer(LayerBase):
    """Scatter in ternary (triangle) space. Config keys:
    - a, b, c: expressions or arrays; if c missing, use 1-a-b
    - mask_outside: bool (default True) → drop points where a<0 or b<0 or c<0
    - draw_frame: bool (default True)
    - draw_grid: bool (default False)
    - grid_ticks: list of floats in (0,1)
    - label_a/b/c: optional corner labels
    - any other kwargs passed to ax.scatter (s, alpha, cmap, c, norm, etc.)
    """
    TYPE = "triangle_scatter"

    def render(self, ax, ctx):
        cfg = self.spec.cfg
        ns  = ctx.namespace
        # Resolve a, b, (optional) c
        a = eval(cfg["a"], ns, {}) if isinstance(cfg.get("a"), str) else cfg.get("a")
        b = eval(cfg["b"], ns, {}) if isinstance(cfg.get("b"), str) else cfg.get("b")
        c = None
        if "c" in cfg:
            c = eval(cfg["c"], ns, {}) if isinstance(cfg.get("c"), str) else cfg.get("c")
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        c = np.asarray(1.0 - a - b, dtype=float) if c is None else np.asarray(c, dtype=float)

        # Optional row-wise normalization if a+b+c != 1
        if cfg.get("auto_normalize", True):
            s = a + b + c
            bad = ~(np.isfinite(s))
            s[bad] = 1.0
            # Only normalize rows where sum is not ~1
            need = np.abs(s - 1.0) > 1e-8
            if np.any(need):
                a = np.where(need, a / s, a)
                b = np.where(need, b / s, b)
                c = np.where(need, c / s, c)

        mask_outside = bool(cfg.get("mask_outside", True))
        if mask_outside:
            m = (a >= 0) & (b >= 0) & (c >= 0)
            kept = int(np.count_nonzero(m))
            if kept == 0:
                print("[triangle_scatter] All points masked out (a,b,c have negatives). Disable mask_outside to view.")
            a, b, c = a[m], b[m], c[m]

        if a.size == 0:
            return None

        x, y = _barycentric_to_xy(a, b, c)

        # triangle adornments
        if cfg.get("draw_frame", True):
            _draw_triangle_frame(ax,
                                 linewidth=cfg.get("frame_linewidth", 1.0),
                                 color=cfg.get("frame_color", "#444"),
                                 alpha=cfg.get("frame_alpha", 1.0))
        if cfg.get("draw_grid", False):
            _draw_grid(ax,
                       ticks=cfg.get("grid_ticks", (0.2, 0.4, 0.6, 0.8)),
                       color=cfg.get("grid_color", "#bbb"),
                       linewidth=cfg.get("grid_linewidth", 0.5),
                       alpha=cfg.get("grid_alpha", 0.6))

        # Merge style: bundle preset → YAML cfg
        preset = ctx.style.layer_preset("triangle_scatter") or {}
        # Merge order: preset < profile.layer defaults < YAML cfg
        prof_defaults = {}
        if hasattr(ctx.style, "layer_defaults_from_profile"):
            prof_defaults = ctx.style.layer_defaults_from_profile(getattr(ctx, "profile", None), "triangle_scatter") or {}

        allowed_skip = {"type","axes","a","b","c","name","draw_frame","draw_grid","mask_outside",
                        "frame_linewidth","frame_color","frame_alpha","grid_ticks","grid_color","grid_linewidth","grid_alpha"}
        style = {k:v for k,v in preset.items() if k not in allowed_skip}
        # profile defaults (if any)
        for k, v in prof_defaults.items():
            if k not in allowed_skip:
                style.setdefault(k, v)
        # YAML overrides last
        for k, v in cfg.items():
            if k not in allowed_skip:
                style[k] = v

        mappable = None
        color_arr = None
        if "c" in style:
            color_arr = style.pop("c")
            if isinstance(color_arr, str):
                color_arr = eval(color_arr, ns, {})
        if color_arr is not None:
            sc = ax.scatter(x, y, c=color_arr, **style)
            mappable = sc
        else:
            sc = ax.scatter(x, y, **style)

        # optional corner labels
        if any(k in cfg for k in ("label_a","label_b","label_c")):
            h = np.sqrt(3.0)/2.0
            if cfg.get("label_a"):
                ax.text(0.0, 0.0, str(cfg["label_a"]), ha="right", va="top")
            if cfg.get("label_b"):
                ax.text(1.0, 0.0, str(cfg["label_b"]), ha="left", va="top")
            if cfg.get("label_c"):
                ax.text(0.5, h, str(cfg["label_c"]), ha="center", va="bottom")

        return mappable
