#!/usr/bin/env python3

import json
import os
from itertools import combinations

import numpy as np
import shapely as sp
import yaml

from jarvishep.base import Base


def create_round_square(ct, rd=1, nn=3.5, frame=None):
    if nn <= 0:
        raise ValueError("nn must be > 0")

    def _signed_power(values, power):
        # Stable for fractional power on negative inputs: sign(x) * |x|**p
        return np.sign(values) * np.power(np.abs(values), power)

    tt = np.linspace(0.0, 2.0 * np.pi, 314)
    power = 2.0 / float(nn)
    xx = _signed_power(np.cos(tt), power)
    yy = _signed_power(np.sin(tt), power)
    xx = rd * xx + ct[0]
    yy = rd * yy + ct[1]
    cdn = np.stack((xx, yy), axis=-1)
    rs = sp.Polygon(cdn)
    if frame is not None:
        rs = rs.intersection(frame)
    return rs


def draw_logo_in_square(ax):
    ax.axis("off")
    ax.set_xlim(-1.02, 1.02)
    ax.set_ylim(-1.02, 1.02)

    outframe = create_round_square((0.0, 0.0), rd=1.0, nn=4.0)
    x_coords = [point[0] for point in outframe.exterior.coords]
    y_coords = [point[1] for point in outframe.exterior.coords]
    ax.fill(x_coords, y_coords, fc="#34495E", ec=None)

    inner = create_round_square((0.0, 0.0), rd=0.76, nn=4.0)
    x_inner = [point[0] for point in inner.exterior.coords]
    y_inner = [point[1] for point in inner.exterior.coords]
    ax.fill(x_inner, y_inner, fc="#00fdff", ec=None, alpha=0.85)

    core = create_round_square((0.0, 0.0), rd=0.35, nn=4.0)
    x_core = [point[0] for point in core.exterior.coords]
    y_core = [point[1] for point in core.exterior.coords]
    ax.fill(x_core, y_core, fc="#072346", ec=None)


class CustomDumper(yaml.Dumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)
        if len(self.indents) == 1:
            super().write_line_break()


class JarvisPLOT(Base):
    def __init__(self):
        super().__init__()
        self.info = {}

    def get_plot_config_from_Jarvis(self, info, scan_yaml):
        """Emit a JarvisPLOT YAML for a Jarvis-HEP run."""
        self.info.update(info)
        self.emit_jplot(scan_yaml)

    def emit_jplot(self, scan_yaml):
        """Emit a JarvisPLOT YAML for current sampling setup."""
        method = None
        try:
            method = scan_yaml.config.get("Sampling", {}).get("Method", None)
        except Exception:
            method = None
        include_dynesty = method in {"Dynesty", "MultiNest"}

        out_yaml_path = self.info["plot"]["config"]
        out_yaml_dir = os.path.abspath(os.path.dirname(out_yaml_path))
        os.makedirs(out_yaml_dir, exist_ok=True)

        plots_dir = os.path.join(out_yaml_dir, ".")
        os.makedirs(plots_dir, exist_ok=True)

        def _rel(path_str):
            try:
                return os.path.relpath(os.path.abspath(path_str), start=out_yaml_dir)
            except Exception:
                return path_str

        if method == "MultiNest":
            dynesty_csv = os.path.join(self.info["sample"]["task_result_dir"], "DATABASE", "multinest_result.csv")
        else:
            dynesty_csv = os.path.join(self.info["sample"]["task_result_dir"], "DATABASE", "dynesty_result.csv")
        has_dynesty_dataset = include_dynesty and os.path.exists(dynesty_csv)

        out_csv = None
        try:
            if not self.info["db"].get("out_csv", False):
                with open(self.info["db"]["info"], "r") as f1:
                    dbinfo = json.loads(f1.read())
                    out_csv = dbinfo.get("converted", None)
                    self.info["db"]["out_csv"] = out_csv
            else:
                out_csv = self.info["db"]["out_csv"]
        except Exception:
            out_csv = None

        if isinstance(out_csv, str):
            out_csv_list = [out_csv]
        elif isinstance(out_csv, list):
            out_csv_list = out_csv
        else:
            out_csv_list = []

        datasets = []
        if has_dynesty_dataset:
            datasets.append({
                "name": "dynesty",
                "path": _rel(dynesty_csv),
                "type": "csv",
            })

        if out_csv_list:
            for path_item in out_csv_list:
                if not path_item:
                    continue
                base = os.path.splitext(os.path.basename(str(path_item)))[0]
                name = base.replace("-", "_").replace(".", "_")
                if not name.startswith("df"):
                    name = f"df_{name}"
                datasets.append({
                    "name": name,
                    "path": _rel(str(path_item)),
                    "type": "csv",
                })
        else:
            datasets.append({"name": "df", "path": "samples.csv", "type": "csv"})

        scan_vars = []
        try:
            scan_vars = scan_yaml.config.get("Sampling", {}).get("Variables", [])
        except Exception:
            scan_vars = []

        var_names = []
        var_meta = {}
        try:
            for var in scan_vars:
                if not isinstance(var, dict) or not var.get("name", None):
                    continue
                name = str(var["name"])
                var_names.append(name)

                vmin, vmax = None, None
                try:
                    params = (var.get("distribution", {}) or {}).get("parameters", {}) or {}
                    vmin = params.get("min", None)
                    vmax = params.get("max", None)
                except Exception:
                    vmin, vmax = None, None

                scale = "linear"
                try:
                    dtype = str((var.get("distribution", {}) or {}).get("type", "")).strip().lower()
                    if "log" in dtype:
                        scale = "log"
                except Exception:
                    scale = "linear"

                var_meta[name] = {
                    "lim": [vmin, vmax],
                    "scale": scale,
                }
        except Exception:
            var_names = []
            var_meta = {}

        if len(var_names) == 0:
            var_names = ["x", "y"]
        elif len(var_names) == 1:
            var_names = [var_names[0], var_names[0]]

        sample_sources = [ds["name"] for ds in datasets if ds["name"] != "dynesty"]
        if not sample_sources:
            sample_sources = ["df"]

        figures = []

        for xx, yy in combinations(var_names, 2):
            safe_x = str(xx).replace("-", "_").replace(".", "_").replace("/", "_")
            safe_y = str(yy).replace("-", "_").replace(".", "_").replace("/", "_")
            fig_name = f"scatter_{safe_x}__{safe_y}"

            figures.append(
                {
                    "name": fig_name,
                    "enable": True,
                    "style": ["a4paper_2x1", "rectcmap"],
                    "frame": {
                        "ax": {
                            "labels": {"x": str(xx), "y": str(yy)},
                            "xlim": (var_meta.get(str(xx), {}) or {}).get("lim", [None, None]),
                            "ylim": (var_meta.get(str(yy), {}) or {}).get("lim", [None, None]),
                            "xscale": (var_meta.get(str(xx), {}) or {}).get("scale", "linear"),
                            "yscale": (var_meta.get(str(yy), {}) or {}).get("scale", "linear"),
                        },
                        "axc": {
                            "label": {"ylabel": "LogL"},
                        },
                    },
                    "layers": [
                        {
                            "name": "scatter",
                            "data": [{"source": src} for src in sample_sources],
                            "axes": "ax",
                            "method": "scatter",
                            "coordinates": {
                                "x": {"expr": str(xx)},
                                "y": {"expr": str(yy)},
                                "c": {"expr": "LogL"},
                            },
                            "style": {
                                "marker": ".",
                                "s": 2,
                                "cmap": "jarvis_rainbow2_r",
                                "norm": "linear",
                                "zorder": 1,
                            },
                        }
                    ],
                }
            )

        if has_dynesty_dataset:
            figures.append(self._build_dynesty_diagnostics_figure())

        jplot_yaml = {
            "DataSet": datasets,
            "Figures": figures,
            "output": {
                "dir": _rel(plots_dir),
                "dpi": 200,
                "formats": ["png"],
            },
            "project": {
                "name": self.info.get("project_name", "Jarvis-HEP"),
            },
            "version": 0.3,
            "Functions": [],
        }

        with open(out_yaml_path, "w") as file:
            yaml.dump(
                jplot_yaml,
                file,
                Dumper=CustomDumper,
                default_flow_style=False,
                allow_unicode=True,
                sort_keys=False,
            )

        self.logger.warning(f"JarvisPLOT YAML generated: {out_yaml_path}")
        self.logger.warning(f"Render with external Jarvis-PLOT CLI: jarvisplot {out_yaml_path}")

    def _build_dynesty_diagnostics_figure(self):
        return {
            "name": "runplot",
            "enable": True,
            "style": ["a4paper_2x1", "rect_5x1"],
            "frame": {
                "ax0": {"labels": {"x": "$-\\ln(X)$", "y": "$N_{\\mathrm{live}}$"}},
                "ax1": {"labels": {"x": "$-\\ln(X)$", "y": "Likelihood"}},
                "ax2": {"labels": {"x": "$-\\ln(X)$", "y": "Importance\nweight PDF"}},
                "ax3": {"labels": {"x": "$-\\ln(X)$", "y": "Evidence"}},
                "ax4": {"labels": {"x": "$-\\ln(X)$", "y": "Iters"}},
            },
            "layers": [
                {
                    "name": "nlive",
                    "data": [{"source": "dynesty"}],
                    "axes": "ax0",
                    "method": "scatter",
                    "coordinates": {"x": {"expr": "-log_PriorVolume"}, "y": {"expr": "samples_nlive"}},
                    "style": {"alpha": 0.7, "zorder": 1},
                },
                {
                    "name": "loglike",
                    "data": [{"source": "dynesty"}],
                    "axes": "ax1",
                    "method": "scatter",
                    "coordinates": {
                        "x": {"expr": "-log_PriorVolume"},
                        "y": {"expr": "np.exp(log_Like) / np.exp(np.max(log_Like))"},
                    },
                    "style": {"alpha": 0.7, "zorder": 1},
                },
                {
                    "name": "Importance\nweight PDF",
                    "data": [{"source": "dynesty"}],
                    "axes": "ax2",
                    "method": "scatter",
                    "coordinates": {
                        "x": {"expr": "-log_PriorVolume"},
                        "y": {"expr": "np.exp(log_weight) / np.exp(np.max(log_weight))"},
                    },
                    "style": {"alpha": 0.7, "zorder": 1},
                },
                {
                    "name": "logevidence",
                    "data": [{"source": "dynesty"}],
                    "axes": "ax3",
                    "method": "scatter",
                    "coordinates": {"x": {"expr": "-log_PriorVolume"}, "y": {"expr": "np.exp(log_Evidence)"}},
                    "style": {"alpha": 0.7, "zorder": 1},
                },
                {
                    "name": "iters",
                    "data": [{"source": "dynesty"}],
                    "axes": "ax4",
                    "method": "scatter",
                    "coordinates": {"x": {"expr": "-log_PriorVolume"}, "y": {"expr": "samples_it"}},
                    "style": {"alpha": 0.7, "zorder": 1},
                },
            ],
        }
