#!/usr/bin/env python3
"""Emit JarvisPLOT scene / jplot YAML from scan outputs (D11.5 / D10.3).

HEP owns DATABASE / levelset → plot input mapping only. Rendering stays in JarvisPLOT
via ``Jarvis2 plot`` or optional post-run auto-render (no PLOT package fork).
"""

from __future__ import annotations

import csv
import json
import os
from collections.abc import Mapping, Sequence
from typing import Any

import yaml


def _write_yaml(path: str, payload: Mapping[str, Any] | list[Any]) -> str:
    output = os.path.abspath(str(path))
    parent = os.path.dirname(output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(output, "w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False, allow_unicode=True)
    return output


def _write_csv(path: str, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]) -> str:
    output = os.path.abspath(str(path))
    parent = os.path.dirname(output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})
    return output


def _load_levelset(levelset_path: str) -> Mapping[str, Any] | None:
    path = os.path.abspath(str(levelset_path))
    if not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, Mapping) else None


def _polylines_physical(payload: Mapping[str, Any]) -> list[list[list[float]]]:
    raw = payload.get("polylines_x") or payload.get("polylines_u") or []
    out: list[list[list[float]]] = []
    for poly in raw:
        if not isinstance(poly, Sequence) or not poly:
            continue
        pts: list[list[float]] = []
        for point in poly:
            if isinstance(point, Mapping):
                # physical map may be dict {name: value}
                vals = list(point.values())
                if len(vals) >= 2:
                    try:
                        pts.append([float(vals[0]), float(vals[1])])
                    except (TypeError, ValueError):
                        continue
                continue
            if not isinstance(point, Sequence) or len(point) < 2:
                continue
            try:
                pts.append([float(point[0]), float(point[1])])
            except (TypeError, ValueError):
                continue
        if len(pts) >= 2:
            out.append(pts)
    return out


def emit_levelset_overlay_scene(
    levelset_path: str,
    *,
    output_yaml: str,
    title: str | None = None,
) -> str | None:
    """Build a lightweight intermediate overlay scene (series form) from levelset.json."""
    payload = _load_levelset(levelset_path)
    if payload is None:
        return None
    polylines = _polylines_physical(payload)
    if not polylines:
        return None

    series: list[dict[str, Any]] = []
    for index, poly in enumerate(polylines):
        xs = [pt[0] for pt in poly]
        ys = [pt[1] for pt in poly]
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
            "source": os.path.abspath(str(levelset_path)),
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


def _read_hdf5_records(db_path: str, *, limit: int | None = None) -> list[dict[str, Any]]:
    path = os.path.abspath(str(db_path))
    if not os.path.isfile(path):
        return []
    try:
        from jarvishep2.database import SimpleHDF5Writer
    except Exception:
        return []
    try:
        records = SimpleHDF5Writer(path).read_records()
    except Exception:
        return []
    out: list[dict[str, Any]] = []
    if limit is None:
        selected = records
    else:
        selected = records[: max(1, int(limit))]
    for row in selected:
        if isinstance(row, Mapping):
            out.append(dict(row))
    return out


def _csv_cell(value: Any) -> Any:
    """Flatten a record value for CSV (scalars as-is; nested → JSON text)."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return value
    if isinstance(value, (str, int, float)):
        return value
    if hasattr(value, "item") and callable(value.item):
        try:
            return value.item()
        except Exception:
            pass
    try:
        return json.dumps(value, ensure_ascii=False, default=str)
    except Exception:
        return str(value)


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
    """Build a simple intermediate scatter scene from DATABASE/samples.hdf5."""
    records = _read_hdf5_records(db_path, limit=limit)
    if not records:
        return None

    xs: list[float] = []
    ys: list[float] = []
    colors: list[float] = []
    for row in records:
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
            "source": os.path.abspath(str(db_path)),
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


def export_samples_csv_from_hdf5(
    db_path: str,
    *,
    output_csv: str,
    x_key: str = "x",
    y_key: str = "y",
    color_key: str = "LogL",
    limit: int | None = None,
) -> str | None:
    """Export DATABASE records to CSV — **full** observables by default.

    Every key present in any HDF5 record becomes a column. Nested values
    (lists/dicts) are JSON-encoded as cell text. Preferred columns
    ``x`` / ``y`` / ``LogL`` / ``uuid`` are ordered first when present for
    jplot convenience; no row is dropped for missing plot keys.

    ``limit=None`` (default) exports all rows; pass a positive int to cap.
    """
    records = _read_hdf5_records(db_path, limit=limit)
    if not records:
        return None

    preferred = [x_key, y_key, color_key, "uuid"]
    fieldnames: list[str] = []
    seen: set[str] = set()
    for key in preferred:
        if any(key in row for row in records):
            fieldnames.append(key)
            seen.add(key)
    # Stable: first-seen order across records, then any late keys.
    for row in records:
        for key in row.keys():
            text = str(key)
            if text not in seen:
                seen.add(text)
                fieldnames.append(text)
    if not fieldnames:
        return None

    rows: list[dict[str, Any]] = []
    for row in records:
        rows.append({key: _csv_cell(row.get(key)) for key in fieldnames})
    return _write_csv(output_csv, rows, fieldnames)


def emit_jplot_scan_levelset_yaml(
    task_result_dir: str,
    *,
    scan_name: str = "scan",
    images_dir: str | None = None,
    project_root: str | None = None,
    x_key: str = "x",
    y_key: str = "y",
    color_key: str = "LogL",
    limit: int = 5000,
) -> str | None:
    """Emit a **stock jplot** YAML: sample scatter + level-set polyline overlay.

    - Plot YAML → ``<project>/images/<scan>/``
    - Samples CSV → next to HDF5 under ``DATABASE/samples.csv`` (full columns)

    Consumable by ``Jarvis2 plot path.yaml`` / ``jplot`` without a V2 adapter.
    """
    from jarvishep2.base import (
        infer_project_root_from_task_result_dir,
        project_images_dir,
    )

    root = os.path.abspath(str(task_result_dir))
    proj = os.path.abspath(
        str(project_root or infer_project_root_from_task_result_dir(root))
    )
    out_dir = os.path.abspath(
        images_dir or project_images_dir(project_root=proj, scan_name=scan_name)
    )
    os.makedirs(out_dir, exist_ok=True)

    db_path = os.path.join(root, "DATABASE", "samples.hdf5")
    # Full Worker evaluation archive (every logL call), not nested dead points.
    # Nested clean CSV is DATABASE/dynesty_result.csv (written by DynestySampler).
    # Caller must wait for Archiver catch-up before this snapshot.
    samples_csv = os.path.join(root, "DATABASE", "samples.csv")
    csv_path = export_samples_csv_from_hdf5(
        db_path,
        output_csv=samples_csv,
        x_key=x_key,
        y_key=y_key,
        color_key=color_key,
        limit=None,
    )

    levelset_path = os.path.join(root, "levelset.json")
    payload = _load_levelset(levelset_path)
    polylines = _polylines_physical(payload) if payload else []

    if not csv_path and not polylines:
        return None

    datasets: list[dict[str, Any]] = []
    layers: list[dict[str, Any]] = []

    if csv_path:
        # Absolute path so jplot works from images/<scan>/ regardless of cwd.
        datasets.append(
            {
                "name": "samples",
                "path": os.path.abspath(csv_path),
                "type": "csv",
            }
        )
        layers.append(
            {
                "name": "scatter",
                "data": [{"source": "samples"}],
                "axes": "ax",
                "method": "scatter",
                "coordinates": {
                    "x": {"expr": x_key},
                    "y": {"expr": y_key},
                    "c": {"expr": color_key},
                },
                "style": {
                    "marker": ".",
                    "s": 6,
                    "zorder": 1,
                },
            }
        )

    for index, poly in enumerate(polylines):
        xs = [pt[0] for pt in poly]
        ys = [pt[1] for pt in poly]
        layers.append(
            {
                "name": f"levelset_{index}",
                "axes": "ax",
                "method": "plot",
                "coordinates": {
                    "x": xs,
                    "y": ys,
                },
                "style": {
                    "linewidth": 1.2,
                    "linestyle": "-",
                    "color": "black",
                    "zorder": 3,
                },
            }
        )

    if not layers:
        return None

    title = f"{scan_name} samples"
    if polylines:
        target = None
        if payload is not None:
            target = payload.get("target_expression")
            tval = payload.get("target_value")
            if target is not None and tval is not None:
                title = f"{scan_name}: {target} = {tval}"
            elif target is not None:
                title = f"{scan_name}: {target}"
        else:
            title = f"{scan_name} + level-set"

    jplot_doc: dict[str, Any] = {
        "DataSet": datasets,
        "Figures": [
            {
                "name": f"{scan_name}_levelset",
                "enable": True,
                "style": ["a4paper_2x1"],
                "frame": {
                    "ax": {
                        "labels": {
                            "x": f"${x_key}$",
                            "y": f"${y_key}$",
                        },
                    },
                    "axc": {
                        "label": {
                            "ylabel": f"${color_key}$",
                        },
                    },
                },
                "layers": layers,
            }
        ],
        # Non-jplot metadata (ignored by plot engine if unknown)
        "_jarvishep2": {
            "producer": "Jarvis-HEP-V2",
            "scan_name": scan_name,
            "levelset": os.path.basename(levelset_path) if payload else None,
            "title": title,
        },
    }
    yaml_path = os.path.join(out_dir, f"{scan_name}_levelset_jplot.yaml")
    return _write_yaml(yaml_path, jplot_doc)


def emit_jplot_dynesty_runplot_yaml(
    task_result_dir: str,
    *,
    scan_name: str = "scan",
    images_dir: str | None = None,
    project_root: str | None = None,
    dynesty_csv: str | None = None,
    source_label: str = "dynesty",
    yaml_basename: str | None = None,
) -> str | None:
    """Emit stock jplot YAML for nested-sampling diagnostics (V1 / Jarvis-PLOT).

    Requires a nested-result CSV (``dynesty_result.csv`` / ``multinest_result.csv``
    or an explicit ``dynesty_csv``). Dataset name **must** be ``dynesty`` so the
    ``dynesty_runplot`` card picks it up without a custom layer block (see
    Jarvis-PLOT DYNESTY_RUNPLOT spec). MultiNest uses the same figure card with
    path pointing at ``multinest_result.csv`` (V1 plot.py parity).

    Layout:
    - CSV: ``<task>/DATABASE/<nested>_result.csv``
    - YAML: ``<project>/images/<scan>/<scan>_*_runplot.yaml``
    """
    from jarvishep2.base import (
        infer_project_root_from_task_result_dir,
        project_images_dir,
    )

    root = os.path.abspath(str(task_result_dir))
    proj = os.path.abspath(
        str(project_root or infer_project_root_from_task_result_dir(root))
    )
    out_dir = os.path.abspath(
        images_dir or project_images_dir(project_root=proj, scan_name=scan_name)
    )
    os.makedirs(out_dir, exist_ok=True)

    csv_path = os.path.abspath(
        str(dynesty_csv or os.path.join(root, "DATABASE", "dynesty_result.csv"))
    )
    if not os.path.isfile(csv_path):
        return None

    # Prefer a readable source basename for metadata / default yaml name.
    csv_base = os.path.basename(csv_path)
    label = str(source_label or "dynesty")
    jplot_doc: dict[str, Any] = {
        "version": 0.3,
        "project": {
            "name": scan_name or "Jarvis-HEP",
        },
        "DataSet": [
            {
                # Style card binds to name "dynesty" regardless of Method.
                "name": "dynesty",
                "path": csv_path,
                "type": "csv",
            }
        ],
        "Figures": [
            {
                "name": "dynesty_runplot",
                "enable": True,
                "style": ["a4paper_2x1", "dynesty_runplot"],
            }
        ],
        "output": {
            "dir": out_dir,
            "dpi": 200,
            "formats": ["png"],
        },
        "Functions": [],
        "_jarvishep2": {
            "producer": "Jarvis-HEP-V2",
            "scan_name": scan_name,
            "source": csv_base,
            "method": label,
            "title": f"{scan_name} {label} diagnostics",
        },
    }
    base = yaml_basename or f"{scan_name}_{label}_runplot.yaml"
    if not base.endswith((".yaml", ".yml")):
        base = f"{base}.yaml"
    yaml_path = os.path.join(out_dir, base)
    return _write_yaml(yaml_path, jplot_doc)


def emit_plot_scenes_from_run(
    task_result_dir: str,
    *,
    scan_name: str = "scan",
    images_dir: str | None = None,
    project_root: str | None = None,
    x_key: str = "x",
    y_key: str = "y",
    color_key: str = "LogL",
    auto_render: bool = False,
) -> dict[str, str]:
    """Emit plot inputs for a finished run.

    Layout:
    - ``<project>/images/<scan>/`` — jplot YAML + intermediate scene YAML / PNG
    - ``<scan>/DATABASE/samples.csv`` — full observables CSV next to ``samples.hdf5``
    - ``<scan>/DATABASE/dynesty_result.csv`` + dynesty_runplot jplot (when present)

    When ``auto_render`` is true and JarvisPLOT is installed, also render the
    jplot YAML to PNG via the stock plot bridge.
    """
    from jarvishep2.base import (
        infer_project_root_from_task_result_dir,
        project_images_dir,
    )

    root = os.path.abspath(str(task_result_dir))
    proj = os.path.abspath(
        str(project_root or infer_project_root_from_task_result_dir(root))
    )
    out_dir = os.path.abspath(
        images_dir or project_images_dir(project_root=proj, scan_name=scan_name)
    )
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
        x_key=x_key,
        y_key=y_key,
        color_key=color_key,
        title=f"{scan_name} samples",
    )
    if scatter:
        written["scatter"] = scatter

    jplot_yaml = emit_jplot_scan_levelset_yaml(
        root,
        scan_name=scan_name,
        images_dir=out_dir,
        project_root=proj,
        x_key=x_key,
        y_key=y_key,
        color_key=color_key,
    )
    if jplot_yaml:
        written["jplot_levelset"] = jplot_yaml
        # Also record the DATABASE CSV if it was written.
        samples_csv = os.path.join(root, "DATABASE", "samples.csv")
        if os.path.isfile(samples_csv):
            written["samples_csv"] = samples_csv
        if auto_render:
            png = _try_render_jplot(jplot_yaml)
            if png:
                written["jplot_levelset_png"] = png

    # Nested-sampling diagnostics (Jarvis-PLOT dynesty_runplot).
    # Dynesty → dynesty_result.csv; MultiNest → multinest_result.csv (same card).
    nested_csv_specs = (
        ("dynesty", "dynesty_result.csv", "jplot_dynesty_runplot", "jplot_dynesty_runplot_png"),
        (
            "multinest",
            "multinest_result.csv",
            "jplot_multinest_runplot",
            "jplot_multinest_runplot_png",
        ),
    )
    for label, csv_name, jplot_key, png_key in nested_csv_specs:
        nested_csv = os.path.join(root, "DATABASE", csv_name)
        if not os.path.isfile(nested_csv):
            continue
        written[f"{label}_result_csv"] = nested_csv
        jplot_nested = emit_jplot_dynesty_runplot_yaml(
            root,
            scan_name=scan_name,
            images_dir=out_dir,
            project_root=proj,
            dynesty_csv=nested_csv,
            source_label=label,
        )
        if jplot_nested:
            written[jplot_key] = jplot_nested
            if auto_render:
                png = _try_render_jplot(jplot_nested)
                if png:
                    written[png_key] = png

    return written


def _try_render_jplot(plot_yaml: str) -> str | None:
    """Best-effort stock JarvisPLOT render; never raises."""
    try:
        from jarvishep2.plot_bridge import run_plot
    except Exception:
        return None
    try:
        code = int(run_plot(plot_yaml))
    except Exception:
        return None
    if code != 0:
        return None
    # Common jplot output locations relative to yaml dir / cwd
    base = os.path.splitext(os.path.abspath(plot_yaml))[0]
    candidates = [
        base + ".png",
        os.path.join(os.path.dirname(base), "images", os.path.basename(base) + ".png"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    # Search sibling images dirs produced by jplot project layout
    parent = os.path.dirname(os.path.abspath(plot_yaml))
    for root, _dirs, files in os.walk(parent):
        for name in files:
            if name.endswith(".png") and os.path.basename(base) in name:
                return os.path.join(root, name)
    return None


__all__ = [
    "emit_jplot_dynesty_runplot_yaml",
    "emit_jplot_scan_levelset_yaml",
    "emit_levelset_overlay_scene",
    "emit_plot_scenes_from_run",
    "emit_scan_scatter_scene_from_hdf5",
    "export_samples_csv_from_hdf5",
]
