#!/usr/bin/env python3
"""Emit JarvisPLOT scene YAML from scan outputs (D11.5).

HEP owns DATABASE / levelset → scene mapping only. Rendering stays in JarvisPLOT.
"""

from __future__ import annotations

import json
import os
from collections.abc import Mapping, Sequence
from typing import Any

import yaml


def _write_yaml(path: str, payload: Mapping[str, Any]) -> str:
    output = os.path.abspath(str(path))
    parent = os.path.dirname(output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(output, "w", encoding="utf-8") as handle:
        yaml.safe_dump(dict(payload), handle, sort_keys=False, allow_unicode=True)
    return output


def emit_levelset_overlay_scene(
    levelset_path: str,
    *,
    output_yaml: str,
    title: str | None = None,
) -> str | None:
    """Build a 2D polyline overlay scene from ``levelset.json`` if polylines exist."""
    path = os.path.abspath(str(levelset_path))
    if not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, Mapping):
        return None
    polylines = payload.get("polylines_x") or payload.get("polylines_u") or []
    if not polylines:
        return None

    series: list[dict[str, Any]] = []
    for index, poly in enumerate(polylines):
        if not isinstance(poly, Sequence) or not poly:
            continue
        xs: list[float] = []
        ys: list[float] = []
        for point in poly:
            if not isinstance(point, Sequence) or len(point) < 2:
                continue
            xs.append(float(point[0]))
            ys.append(float(point[1]))
        if len(xs) < 2:
            continue
        series.append(
            {
                "name": f"contour_{index}",
                "x": xs,
                "y": ys,
                "style": "line",
            }
        )
    if not series:
        return None

    scene = {
        "schema": "jarvisplot.scene/v1",
        "scene_type": "overlay",
        "scene_id": "levelset_overlay",
        "metadata": {
            "producer": "Jarvis-HEP-V2",
            "source": path,
            "title": title or "Adaptive level-set contour",
            "target_expression": payload.get("target_expression"),
            "target_value": payload.get("target_value"),
        },
        "figure": {
            "xlabel": "x",
            "ylabel": "y",
            "title": title or "Adaptive level-set contour",
        },
        "series": series,
    }
    return _write_yaml(output_yaml, scene)


def emit_scan_scatter_scene_from_hdf5(
    db_path: str,
    *,
    output_yaml: str,
    x_key: str = "x",
    y_key: str = "y",
    color_key: str = "LogL",
    limit: int = 5000,
    title: str | None = None,
) -> str | None:
    """Build a simple scatter scene from DATABASE/samples.hdf5 observables."""
    path = os.path.abspath(str(db_path))
    if not os.path.isfile(path):
        return None
    try:
        from jarvishep2.database import SimpleHDF5Writer
    except Exception:
        return None
    try:
        records = SimpleHDF5Writer(path).read_records()
    except Exception:
        return None
    if not records:
        return None

    xs: list[float] = []
    ys: list[float] = []
    colors: list[float] = []
    for row in records[: max(1, int(limit))]:
        if not isinstance(row, Mapping):
            continue
        if x_key not in row or y_key not in row:
            continue
        try:
            xs.append(float(row[x_key]))
            ys.append(float(row[y_key]))
            if color_key in row:
                colors.append(float(row[color_key]))
        except (TypeError, ValueError):
            continue
    if len(xs) < 1:
        return None

    series: dict[str, Any] = {
        "name": "samples",
        "x": xs,
        "y": ys,
        "style": "scatter",
    }
    if len(colors) == len(xs):
        series["c"] = colors
        series["colorbar_label"] = color_key

    scene = {
        "schema": "jarvisplot.scene/v1",
        "scene_type": "scatter",
        "scene_id": "scan_scatter",
        "metadata": {
            "producer": "Jarvis-HEP-V2",
            "source": path,
            "title": title or "Scan samples",
            "n_points": len(xs),
        },
        "figure": {
            "xlabel": x_key,
            "ylabel": y_key,
            "title": title or "Scan samples",
        },
        "series": [series],
    }
    return _write_yaml(output_yaml, scene)


def emit_plot_scenes_from_run(
    task_result_dir: str,
    *,
    scan_name: str = "scan",
    images_dir: str | None = None,
) -> dict[str, str]:
    """Emit available plot scene YAMLs under ``images/`` for a finished run."""
    root = os.path.abspath(str(task_result_dir))
    out_dir = os.path.abspath(images_dir or os.path.join(root, "images"))
    os.makedirs(out_dir, exist_ok=True)
    written: dict[str, str] = {}

    levelset = os.path.join(root, "levelset.json")
    overlay = emit_levelset_overlay_scene(
        levelset,
        output_yaml=os.path.join(out_dir, f"{scan_name}_levelset.yaml"),
        title=f"{scan_name} level-set",
    )
    if overlay:
        written["levelset_overlay"] = overlay

    db_path = os.path.join(root, "DATABASE", "samples.hdf5")
    scatter = emit_scan_scatter_scene_from_hdf5(
        db_path,
        output_yaml=os.path.join(out_dir, f"{scan_name}_scatter.yaml"),
        title=f"{scan_name} samples",
    )
    if scatter:
        written["scatter"] = scatter
    return written


__all__ = [
    "emit_levelset_overlay_scene",
    "emit_plot_scenes_from_run",
    "emit_scan_scatter_scene_from_hdf5",
]
