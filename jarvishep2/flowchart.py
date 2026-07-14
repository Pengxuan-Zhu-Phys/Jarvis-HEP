#!/usr/bin/env python3
"""Workflow flowchart scene export from V2 execution plans (D11.5 / D12.3 precursor).

Produces a ``jarvisplot.scene/v1`` flowchart scene that JarvisPLOT can render.
This is a plan-based exporter (layers + module nodes), not a full V1 IO-port graph.
"""

from __future__ import annotations

import json
import os
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.sample import ExecutionStep
from jarvishep2.workflow import build_execution_plan, group_by_layer

FLOWCHART_SCHEMA = "jarvisplot.scene/v1"
FLOWCHART_SCENE_TYPE = "flowchart"
FLOWCHART_SCENE_ID = "workflow_main"


def _role_for_step_type(step_type: str) -> str:
    normalized = str(step_type or "").strip().lower()
    if normalized == "calculator":
        return "calculator"
    if normalized == "opera":
        return "operas"
    if normalized == "likelihood":
        return "likelihood"
    return normalized or "module"


def build_flowchart_scene(
    steps: Sequence[ExecutionStep] | Sequence[Mapping[str, Any]],
    *,
    scan_name: str = "scan",
    project_name: str = "",
    producer: str = "Jarvis-HEP-V2",
) -> dict[str, Any]:
    """Build a flowchart scene dict from an execution plan."""
    normalized: list[ExecutionStep] = []
    for item in steps:
        if isinstance(item, ExecutionStep):
            normalized.append(item)
        elif isinstance(item, Mapping):
            normalized.append(
                ExecutionStep(
                    type=str(item.get("type", "module")),
                    name=str(item.get("name", "Module")),
                    layer=int(item.get("layer", 0)),
                    params=dict(item.get("params") or {}),
                )
            )
    layers_steps = group_by_layer(normalized)

    nodes: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    layer_entries: list[dict[str, Any]] = []
    prev_layer_node_ids: list[str] = []

    # Synthetic parameter source
    param_id = "Parameters"
    nodes.append(
        {
            "id": param_id,
            "kind": "source",
            "role": "parameter_source",
            "label": "Parameters",
            "layer": "layer_0",
            "out_ports": [{"id": "out", "role": "parameter"}],
        }
    )
    layer_entries.append(
        {
            "id": "layer_0",
            "index": 0,
            "label": "Parameters",
            "nodes": [param_id],
        }
    )
    prev_layer_node_ids = [param_id]

    for layer_steps in layers_steps:
        if not layer_steps:
            continue
        layer_index = int(layer_steps[0].layer) + 1  # shift after Parameters
        layer_id = f"layer_{layer_index}"
        layer_node_ids: list[str] = []
        labels: list[str] = []
        for step in layer_steps:
            node_id = str(step.name)
            labels.append(node_id)
            layer_node_ids.append(node_id)
            nodes.append(
                {
                    "id": node_id,
                    "kind": "module",
                    "role": _role_for_step_type(step.type),
                    "label": node_id,
                    "layer": layer_id,
                    "in_ports": [{"id": "in", "role": "data"}],
                    "out_ports": [{"id": "out", "role": "data"}],
                    "metadata": {"module_type": step.type},
                }
            )
            for source_id in prev_layer_node_ids:
                edges.append(
                    {
                        "source": {"node": source_id, "port": "out"},
                        "target": {"node": node_id, "port": "in"},
                        "role": "dataflow",
                    }
                )
        layer_entries.append(
            {
                "id": layer_id,
                "index": layer_index,
                "label": ", ".join(labels) if labels else f"Layer {layer_index}",
                "nodes": layer_node_ids,
            }
        )
        prev_layer_node_ids = layer_node_ids

    return {
        "schema": FLOWCHART_SCHEMA,
        "scene_type": FLOWCHART_SCENE_TYPE,
        "scene_id": FLOWCHART_SCENE_ID,
        "metadata": {
            "producer": producer,
            "scan_name": str(scan_name or "scan"),
            "project_name": str(project_name or ""),
        },
        "layers": layer_entries,
        "nodes": nodes,
        "edges": edges,
    }


def build_flowchart_scene_from_config(
    config: Mapping[str, Any],
    *,
    include_likelihood: bool = True,
) -> dict[str, Any]:
    """Build a flowchart scene from a task config's Calculators/Operas blocks."""
    calculators = (config.get("Calculators") or {}).get("Modules") or []
    operas = (config.get("Operas") or {}).get("Modules") or []
    calc_list = [dict(item) for item in calculators if isinstance(item, Mapping)]
    opera_list = [dict(item) for item in operas if isinstance(item, Mapping)]
    steps = build_execution_plan(
        calculator_modules=calc_list,
        opera_modules=opera_list,
        include_likelihood=include_likelihood and bool(calc_list or opera_list),
    )
    scan = config.get("Scan") if isinstance(config.get("Scan"), Mapping) else {}
    scan_name = str(
        (scan or {}).get("name") or config.get("scan_name") or "scan"
    )
    project_name = str(config.get("project_name") or "")
    return build_flowchart_scene(
        steps,
        scan_name=scan_name,
        project_name=project_name,
    )


def export_flowchart_semantics(
    scene: Mapping[str, Any],
    save_path: str,
) -> str:
    """Write flowchart scene JSON; return absolute path."""
    output_path = os.path.abspath(str(save_path))
    parent = os.path.dirname(output_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(dict(scene), handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return output_path


def render_flowchart_png(
    scene: Mapping[str, Any],
    output_path: str,
) -> str | None:
    """Render a flowchart scene via JarvisPLOT if available. Returns PNG path or None."""
    try:
        from jarvisplot import render_flowchart
    except ImportError:
        return None
    out = os.path.abspath(str(output_path))
    parent = os.path.dirname(out)
    if parent:
        os.makedirs(parent, exist_ok=True)
    rendered = render_flowchart(dict(scene), output_path=out)
    return str(rendered or out)


__all__ = [
    "FLOWCHART_SCHEMA",
    "FLOWCHART_SCENE_ID",
    "FLOWCHART_SCENE_TYPE",
    "build_flowchart_scene",
    "build_flowchart_scene_from_config",
    "export_flowchart_semantics",
    "render_flowchart_png",
]
