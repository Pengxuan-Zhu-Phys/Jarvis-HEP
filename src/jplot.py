from __future__ import annotations

# #!/usr/bin/env python3
# """
# JPlot — a lightweight, standalone plotting engine for JarvisPLOT-style YAML jobs.

# Goals
# -----
# - **Independent** from Jarvis-HEP internals (no Base, no BudingPLOT dependency).
# - **Axes-first** composition: a Figure owns named axes (free-form layout via rects).
# - **Layer registry** with a few core layers: scatter, contour, ring_guides, colorbar.
# - **Safe expression** evaluation against loaded datasources (numpy/pandas namespace).
# - **JSON presets + YAML overrides ready** (hooks provided; optional in minimal demo).

# Usage
# -----
# python -m jplot path/to/plot.yaml

# YAML (minimal sketch)
# ---------------------
# figures:
#   - name: demo
#     size: [8, 5]
#     dpi: 150
#     axes:
#       - { name: main, rect: [0.1, 0.12, 0.65, 0.8], projection: cartesian }
#       - { name: cbar, rect: [0.80, 0.12, 0.03, 0.8], projection: cartesian, type: colorbar, from_layer: scatter1, orientation: vertical }
#     datasources:
#       - { name: df, type: csv, path: ./data.csv }
#     layers:
#       - { name: scatter1, type: scatter, axes: main, x: "df.x", y: "df.y", c: "df.L", cmap: viridis, s: 10, alpha: 0.8 }
#       - { name: contourL, type: contour, axes: main, x: "grid.X", y: "grid.Y", z: "grid.L", levels: 12 }
#       - { type: ring_guides, axes: main, rings: [ {center: [0,0], radius: 1.0} ] }
#       - { type: colorbar, axes: cbar, from_layer: scatter1, orientation: vertical }

# This file is intentionally small but structured to be extended.
# """
# from __future__ import annotations

# import os
# import sys
# import json
# import math
# import types
# import typing as T

# import numpy as np
# import pandas as pd
# import yaml
# import matplotlib
# import matplotlib.pyplot as plt
# from matplotlib.colors import Normalize, LogNorm

# # -------------------------
# # Utilities
# # -------------------------

# SAFE_GLOBALS_BASE = {
#     "np": np,
#     "numpy": np,
#     "pd": pd,
#     "math": math,
# }


# def deep_merge(base: dict, *overrides: dict) -> dict:
#     """Recursive dict merge (rightmost wins)."""
#     import copy
#     out = copy.deepcopy(base)
#     for ov in overrides:
#         if not ov:
#             continue
#         for k, v in ov.items():
#             if isinstance(v, dict) and isinstance(out.get(k), dict):
#                 out[k] = deep_merge(out[k], v)
#             else:
#                 out[k] = v
#     return out


# # -------------------------
# # Data environment
# # -------------------------
# class DataEnv:
#     """Loads named datasources and provides an eval() namespace.

#     Supports: CSV (pandas.read_csv). `yaml_dir` is used to resolve relative paths.
#     """

#     def __init__(self, dataspecs: T.List[dict] | None, yaml_dir: str | None = None):
#         self.frames: dict[str, pd.DataFrame] = {}
#         self.yaml_dir = yaml_dir or os.getcwd()
#         if dataspecs:
#             for ds in dataspecs:
#                 name = ds["name"]
#                 typ = ds.get("type", "csv").lower()
#                 path = ds.get("path")
#                 if path is None:
#                     raise ValueError(f"Datasource '{name}' is missing a 'path'")
#                 if not os.path.isabs(path):
#                     path = os.path.join(self.yaml_dir, path)
#                 if typ == "csv":
#                     self.frames[name] = pd.read_csv(path)
#                 else:
#                     raise ValueError(f"Unsupported datasource type: {typ}")

#     def namespace(self) -> dict:
#         ns = dict(SAFE_GLOBALS_BASE)
#         # expose each dataframe under its name; columns accessible as ns['df']['col'] via eval strings like "df.x"
#         for k, v in self.frames.items():
#             ns[k] = v
#         return ns

#     def eval_expr(self, expr: str) -> np.ndarray | pd.Series:
#         ns = self.namespace()
#         try:
#             val = eval(expr, ns, {})
#         except Exception as e:
#             raise RuntimeError(f"Eval failed for expression '{expr}': {e}")
#         # Accept pandas Series/ndarray/scalars
#         if isinstance(val, (pd.Series, np.ndarray, list, tuple)):
#             return np.asarray(val)
#         return np.asarray([val])


# # -------------------------
# # Axes factory
# # -------------------------
# class AxesFactory:
#     @staticmethod
#     def create(fig: plt.Figure, spec: dict) -> plt.Axes:
#         rect = spec.get("rect", [0.1, 0.1, 0.8, 0.8])
#         proj = spec.get("projection", "cartesian")
#         if proj == "polar":
#             ax = fig.add_axes(rect, projection="polar")
#             polar_cfg = spec.get("polar", {})
#             if "theta_zero_location" in polar_cfg:
#                 ax.set_theta_zero_location(polar_cfg["theta_zero_location"])  # e.g. 'N'
#             if "theta_direction" in polar_cfg:
#                 ax.set_theta_direction(polar_cfg["theta_direction"])  # 1 or -1
#         elif proj == "ternary":
#             # Placeholder: try mpltern if available, else cartesian fallback with a note
#             try:
#                 import mpltern  # noqa: F401
#                 ax = fig.add_axes(rect, projection="ternary")
#             except Exception:
#                 ax = fig.add_axes(rect)
#                 ax.text(0.5, 0.5, "[ternary backend missing]", ha="center", va="center", transform=ax.transAxes, fontsize=9, color="#999")
#         else:
#             ax = fig.add_axes(rect)
#         # common titles/labels
#         if title := spec.get("title"):
#             ax.set_title(title)
#         if xlabel := spec.get("xlabel"):
#             ax.set_xlabel(xlabel)
#         if ylabel := spec.get("ylabel"):
#             ax.set_ylabel(ylabel)
#         if face := spec.get("facecolor"):
#             ax.set_facecolor(face)
#         if spec.get("grid", False):
#             ax.grid(True, which="both", alpha=0.3)
#         return ax


# # -------------------------
# # Layer registry
# # -------------------------
# class Layer:
#     def __init__(self, name: str | None, cfg: dict):
#         self.name = name
#         self.cfg = cfg
#         self.mappable = None  # for colorbar

#     def render(self, ax: plt.Axes, data: DataEnv, context: dict):
#         raise NotImplementedError


# class ScatterLayer(Layer):
#     def render(self, ax: plt.Axes, data: DataEnv, context: dict):
#         # evaluate coordinates
#         x = data.eval_expr(self.cfg["x"])  # e.g., "df.x"
#         y = data.eval_expr(self.cfg["y"])  # e.g., "df.y"
#         style = {k: v for k, v in self.cfg.items() if k not in {"type", "axes", "x", "y", "name"}}
#         cexpr = self.cfg.get("c")
#         if cexpr is not None:
#             c = data.eval_expr(cexpr)
#             # color scaling
#             vmin = self.cfg.get("vmin")
#             vmax = self.cfg.get("vmax")
#             norm = None
#             if self.cfg.get("scale") == "log":
#                 vmin = vmin if vmin is not None else max(np.min(c[c>0]) if np.any(c>0) else 1e-6, 1e-12)
#                 vmax = vmax if vmax is not None else np.max(c)
#                 norm = LogNorm(vmin=vmin, vmax=vmax)
#             # filter style so we don't pass a duplicate keyword for color 'c'
#             style_filtered = {k: v for k, v in style.items() if k not in {"c", "vmin", "vmax", "scale"}}
#             sc = ax.scatter(x, y, c=c, norm=norm, vmin=vmin, vmax=vmax, **style_filtered)
#             self.mappable = sc
#         else:
#             sc = ax.scatter(x, y, **style)
#             self.mappable = sc
#         context["_last_mappable"] = self.mappable


# class ContourLayer(Layer):
#     def render(self, ax: plt.Axes, data: DataEnv, context: dict):
#         X = data.eval_expr(self.cfg["x"])  # can be 1D unique grid or flattened
#         Y = data.eval_expr(self.cfg["y"])  # same
#         Z = data.eval_expr(self.cfg["z"])  # same length as X/Y
#         # attempt to reshape if grid dims are provided
#         nx = self.cfg.get("nx")
#         ny = self.cfg.get("ny")
#         if nx and ny:
#             X = X.reshape(ny, nx)
#             Y = Y.reshape(ny, nx)
#             Z = Z.reshape(ny, nx)
#         levels = self.cfg.get("levels", 10)
#         cs = ax.contour(X, Y, Z, levels=levels, linewidths=self.cfg.get("linewidth", 0.8), alpha=self.cfg.get("alpha", 0.9), cmap=self.cfg.get("cmap"))
#         if self.cfg.get("clabel"):
#             ax.clabel(cs, inline=True, fontsize=self.cfg.get("clabel_fontsize", 8))
#         self.mappable = cs
#         context["_last_mappable"] = self.mappable


# class RingGuidesLayer(Layer):
#     def render(self, ax: plt.Axes, data: DataEnv, context: dict):
#         rings = self.cfg.get("rings", [])
#         import matplotlib.patches as mpatches
#         style = {k: v for k, v in self.cfg.items() if k not in {"type", "axes", "rings", "name"}}
#         for r in rings:
#             (cx, cy) = r.get("center", [0.0, 0.0])
#             rad = r.get("radius", 1.0)
#             circ = mpatches.Circle((cx, cy), rad, fill=False, **style)
#             ax.add_patch(circ)


# class ColorbarLayer(Layer):
#     def render(self, ax: plt.Axes, data: DataEnv, context: dict):
#         from_layer = self.cfg.get("from_layer")
#         orientation = self.cfg.get("orientation", "vertical")
#         mappable = None
#         if from_layer:
#             mappable = context.get("mappables", {}).get(from_layer)
#         if mappable is None:
#             mappable = context.get("_last_mappable")
#         if mappable is None:
#             ax.text(0.5, 0.5, "[no mappable]", transform=ax.transAxes, ha="center", va="center", fontsize=9, color="#999")
#             return
#         cb = plt.colorbar(mappable, cax=ax, orientation=orientation)
#         if label := self.cfg.get("label"):
#             if orientation == "vertical":
#                 ax.set_ylabel(label)
#             else:
#                 ax.set_xlabel(label)


# LAYER_REGISTRY: dict[str, T.Type[Layer]] = {
#     "scatter": ScatterLayer,
#     "contour": ContourLayer,
#     "ring_guides": RingGuidesLayer,
#     "colorbar": ColorbarLayer,
# }


# # -------------------------
# # Figure renderer
# # -------------------------
# class JPlotEngine:
#     def __init__(self, cfg: dict):
#         self.cfg = cfg
#         self.yaml_dir = os.path.dirname(os.path.abspath(cfg.get("_yaml_path", ".")))

#     def run(self):
#         figures = self.cfg.get("figures", [])
#         out_dir_cfg = self.cfg.get("output", {}).get("dir", "./IMAGE")
#         # Resolve output directory relative to the YAML file location if not absolute
#         out_dir = out_dir_cfg if os.path.isabs(out_dir_cfg) else os.path.abspath(os.path.join(self.yaml_dir, out_dir_cfg))
#         dpi = self.cfg.get("output", {}).get("dpi", 150)
#         fmts = self.cfg.get("output", {}).get("formats", ["png"])
#         os.makedirs(out_dir, exist_ok=True)

#         for fig_spec in figures:
#             name = fig_spec.get("name", "figure")
#             size = fig_spec.get("size", [8, 5])
#             fig = plt.figure(figsize=size, dpi=dpi)

#             # datasources per figure (optional)
#             env = DataEnv(fig_spec.get("datasources"), yaml_dir=self.yaml_dir)

#             # create axes by name
#             axes_map: dict[str, plt.Axes] = {}
#             for ax_spec in fig_spec.get("axes", []):
#                 ax = AxesFactory.create(fig, ax_spec)
#                 axes_map[ax_spec["name"]] = ax

#             # context stores mappables for colorbars, etc.
#             context: dict = {"mappables": {}}

#             # render layers in order
#             for layer_spec in fig_spec.get("layers", []):
#                 ltype = layer_spec["type"]
#                 lname = layer_spec.get("name")
#                 target = layer_spec.get("axes")
#                 layer_cls = LAYER_REGISTRY.get(ltype)
#                 if layer_cls is None:
#                     raise ValueError(f"Unknown layer type: {ltype}")
#                 layer = layer_cls(lname, layer_spec)
#                 ax = axes_map[target] if target else plt.gca()
#                 layer.render(ax, env, context)
#                 if lname and getattr(layer, "mappable", None) is not None:
#                     context["mappables"][lname] = layer.mappable

#             # save outputs
#             saved_paths = []
#             for fmt in fmts:
#                 path = os.path.join(out_dir, f"{name}.{fmt}")
#                 fig.savefig(path, dpi=dpi)
#                 saved_paths.append(path)
#             plt.close(fig)
#             print("Saved figure(s):")
#             for p in saved_paths:
#                 print("  ", os.path.abspath(p))


# # -------------------------
# # CLI
# # -------------------------

# def load_yaml(path: str) -> dict:
#     with open(path, "r") as f:
#         cfg = yaml.safe_load(f)
#     # remember where this YAML came from for relative path resolution
#     cfg["_yaml_path"] = os.path.abspath(path)
#     return cfg


# def main(argv: list[str] | None = None):
#     argv = sys.argv[1:] if argv is None else argv
#     if not argv:
#         print("Usage: python -m jplot <plot.yaml>")
#         return 2
#     cfg = load_yaml(argv[0])
#     JPlotEngine(cfg).run()
#     return 0


# if __name__ == "__main__":
#     raise SystemExit(main())
