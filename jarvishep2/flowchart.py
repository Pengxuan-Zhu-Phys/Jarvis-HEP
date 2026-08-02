#!/usr/bin/env python3
"""Workflow flowchart scene export (D12.3 / V1 semantic parity).

Builds a ``jarvisplot.scene/v1`` flowchart scene from a task YAML config so that
stock JarvisPLOT can render it with no V2-specific adapter:

  jarvisplot.render_flowchart(scene)  # or ``jplot flowchart flowchart.json``

The graph shape matches V1 ``Workflow.build_flowchart_semantics``: Parameters
source, variable nodes, file I/O nodes, module ports, bridges, selectionflow.
"""

from __future__ import annotations

import json
import os
import re
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.calculator_modes import expand_calculator_modes, mode_info
from jarvishep2.sample import ExecutionStep
from jarvishep2.versioning import get_runtime_version
from jarvishep2.workflow import resolve_module_layers

FLOWCHART_SCHEMA = "jarvisplot.scene/v1"
FLOWCHART_SCENE_TYPE = "flowchart"
FLOWCHART_SCENE_ID = "workflow_main"

_EDGE_ROLES = frozenset(
    {"parameterflow", "dataflow", "fileflow", "bridgeflow", "selectionflow"}
)
_KNOWN_CONSTANTS = frozenset(
    {
        "Pi",
        "pi",
        "E",
        "e",
        "I",
        "True",
        "False",
        "None",
        "nan",
        "inf",
        "Infinity",
    }
)


# ---------------------------------------------------------------------------
# JSON / port helpers (V1 Workflow static methods)
# ---------------------------------------------------------------------------


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, set):
        return [_jsonable(item) for item in sorted(value, key=str)]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    try:
        import numpy as np

        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
    except Exception:
        pass
    return str(value)


def _layer_name(layer_index: int) -> str:
    return f"layer_{int(layer_index)}"


def _edge_role(role: str) -> str:
    normalized = str(role).strip().lower()
    return normalized if normalized in _EDGE_ROLES else "dataflow"


def _port(port_id: Any, role: str, label: Any = None, metadata: Any = None) -> dict[str, Any]:
    port: dict[str, Any] = {"id": str(port_id), "role": str(role)}
    if label is not None and str(label) != str(port_id):
        port["label"] = str(label)
    if metadata:
        port["metadata"] = _jsonable(metadata)
    return port


def _unique_names(values: Any) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values or []:
        text = str(value)
        if not text or text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    return ordered


def _expression_free_names(expression: str) -> list[str]:
    """Best-effort free-symbol names (V1 Dump ``inc`` semantics)."""
    text = str(expression or "").strip()
    if not text:
        return []
    try:
        from jarvishep2.expression import ExpressionContext

        compiled = ExpressionContext(maxsize=32).compile(text)
        return list(compiled.variable_names)
    except Exception:
        pass
    # Fallback: identifier scan, drop common constants / bare numbers.
    names: list[str] = []
    for token in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", text):
        if token in _KNOWN_CONSTANTS:
            continue
        names.append(token)
    return _unique_names(names)


def _input_source_names(spec: Any) -> list[str]:
    """Names of upstream variables that feed this input port (V1 parity).

    For Operas/calculator ports like ``{name: x, expression: "xx * Pi"}``, the
    **port** is ``x`` but the **sources** are free symbols of the expression
    (``xx``), not the port name. Falling back to ``name`` wrongly invents
    ``var::x`` and disconnects ``var::xx`` (see EggBox_Bridson_05 vs V2).
    """
    source_names: list[str] = []
    if isinstance(spec, str):
        return [str(spec)]
    if not isinstance(spec, Mapping):
        return [str(spec)]

    # Top-level expression on the port itself (typical Operas input).
    if "expression" in spec and str(spec.get("expression") or "").strip():
        free = _expression_free_names(str(spec["expression"]))
        if free:
            source_names.extend(free)

    variables = spec.get("variables")
    if isinstance(variables, Mapping):
        for item in variables.values():
            if isinstance(item, Mapping) and "inc" in item:
                source_names.extend(str(val) for val in item.get("inc") or [])
            elif isinstance(item, Mapping) and item.get("name"):
                if "expression" in item:
                    source_names.extend(_expression_free_names(str(item["expression"])))
                else:
                    source_names.append(str(item["name"]))
            elif isinstance(item, str):
                source_names.append(str(item))
    elif isinstance(variables, list):
        for item in variables:
            if isinstance(item, Mapping):
                if "expression" in item:
                    free = _expression_free_names(str(item["expression"]))
                    if free:
                        source_names.extend(free)
                    elif item.get("name"):
                        source_names.append(str(item["name"]))
                elif item.get("name"):
                    source_names.append(str(item["name"]))
            elif isinstance(item, str):
                source_names.append(str(item))

    # V1 calculator analyze_config rewrites Dump vars into dict with ``inc``.
    actions = spec.get("actions")
    if isinstance(actions, list):
        for action in actions:
            if not isinstance(action, Mapping):
                continue
            action_type = str(action.get("type") or "").strip()
            if action_type == "File" and action.get("source"):
                source_names.append(str(action["source"]))
            vars_list = action.get("variables")
            if not isinstance(vars_list, list):
                continue
            for var in vars_list:
                if not isinstance(var, Mapping):
                    continue
                if "expression" in var:
                    source_names.extend(_expression_free_names(str(var["expression"])))
                elif var.get("name"):
                    source_names.append(str(var["name"]))

    if not source_names and spec.get("_inc"):
        source_names.extend(str(val) for val in spec.get("_inc") or [])
    if not source_names and spec.get("entry"):
        source_names.append(str(spec["entry"]))
    # Only if there is still no expression-derived source: use the port name
    # (plain passthrough ports without expression).
    if not source_names and spec.get("name") and "expression" not in spec:
        source_names.append(str(spec["name"]))
    return _unique_names(source_names)


def _output_produced_names(spec: Any) -> list[str]:
    output_names: list[str] = []
    if isinstance(spec, str):
        return [str(spec)]
    if not isinstance(spec, Mapping):
        return [str(spec)]

    variables = spec.get("variables")
    if isinstance(variables, Mapping):
        output_names.extend(str(key) for key in variables.keys())
    elif isinstance(variables, list):
        for item in variables:
            if isinstance(item, Mapping) and item.get("name"):
                output_names.append(str(item["name"]))
            elif isinstance(item, str):
                output_names.append(str(item))

    if not output_names and spec.get("name"):
        output_names.append(str(spec["name"]))
    return _unique_names(output_names)


def _spec_is_file(spec: Any, module_type: str, *, direction: str) -> bool:
    if isinstance(spec, str) or not isinstance(spec, Mapping):
        return False
    spec_type = str(spec.get("type", "") or "").strip().lower()
    if spec.get("path"):
        return True
    if spec_type in {"file", "slha", "xslha", "json", "csv", "tsv", "dat", "wolfram"}:
        return True
    if direction == "input":
        return (
            module_type == "Calculator"
            and bool(spec.get("variables") or spec.get("actions"))
            and not spec.get("expression")
        )
    return bool(spec.get("variables"))


def _module_role(module_type: str) -> tuple[str, str]:
    normalized = str(module_type or "Module").strip().lower().replace(" ", "_")
    if normalized == "parameter":
        return "source", "parameter_source"
    if normalized in {"operas", "opera"}:
        return "module", "operas"
    if normalized == "calculator":
        return "module", "calculator"
    if normalized == "likelihood":
        return "module", "likelihood"
    return "module", normalized or "module"


# ---------------------------------------------------------------------------
# Graph mutation helpers
# ---------------------------------------------------------------------------


def _add_node(
    nodes_by_id: dict[str, dict[str, Any]],
    layer_nodes: dict[int, list[str]],
    node_id: str,
    *,
    kind: str,
    role: str,
    layer: int,
    label: str,
    metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    node = nodes_by_id.get(node_id)
    if node is None:
        node = {
            "id": str(node_id),
            "kind": str(kind),
            "role": str(role),
            "layer_index": int(layer),
            "label": str(label),
            "in_ports": [],
            "out_ports": [],
        }
        if metadata:
            node["metadata"] = _jsonable(dict(metadata))
        nodes_by_id[node_id] = node
        layer_nodes.setdefault(int(layer), []).append(node_id)
        return node

    if int(layer) < int(node.get("layer_index", layer)):
        node["layer_index"] = int(layer)
    if metadata:
        existing = node.get("metadata", {})
        if not isinstance(existing, dict):
            existing = {}
        merged = dict(existing)
        merged.update(_jsonable(dict(metadata)))
        node["metadata"] = merged
    return node


def _add_port(node: dict[str, Any], port_direction: str, port: Mapping[str, Any]) -> None:
    key = f"{port_direction}:{port['id']}"
    port_keys = node.setdefault("_port_keys", set())
    if key in port_keys:
        return
    node[f"{port_direction}_ports"].append(dict(port))
    port_keys.add(key)


def _add_edge(
    edges: list[dict[str, Any]],
    edge_keys: set[tuple[Any, ...]],
    source_node: str,
    source_port: str,
    target_node: str,
    target_port: str,
    role: str,
    metadata: Mapping[str, Any] | None = None,
) -> None:
    role = _edge_role(role)
    key = (source_node, source_port, target_node, target_port, role)
    if key in edge_keys:
        return
    edge: dict[str, Any] = {
        "source": {"node": str(source_node), "port": str(source_port)},
        "target": {"node": str(target_node), "port": str(target_port)},
        "role": str(role),
    }
    if metadata:
        edge["metadata"] = _jsonable(dict(metadata))
    edges.append(edge)
    edge_keys.add(key)


def _bridge_chain(
    nodes_by_id: dict[str, dict[str, Any]],
    layer_nodes: dict[int, list[str]],
    edges: list[dict[str, Any]],
    edge_keys: set[tuple[Any, ...]],
    source_node_id: str,
    source_port_id: str,
    source_layer: int,
    target_node_id: str,
    target_port_id: str,
    target_layer: int,
    role: str,
    bridge_label: str,
    metadata: Mapping[str, Any] | None = None,
) -> None:
    if int(target_layer) - int(source_layer) <= 1:
        _add_edge(
            edges,
            edge_keys,
            source_node_id,
            source_port_id,
            target_node_id,
            target_port_id,
            role,
            metadata=metadata,
        )
        return

    prev_node = source_node_id
    prev_port = source_port_id
    for layer in range(int(source_layer) + 1, int(target_layer)):
        bridge_id = f"bridge::{bridge_label}::L{layer}"
        bridge_node = _add_node(
            nodes_by_id,
            layer_nodes,
            bridge_id,
            kind="bridge",
            role="bridge_relay",
            layer=layer,
            label=bridge_label,
            metadata={
                "source_layer": int(source_layer),
                "target_layer": int(target_layer),
                "relay_for": bridge_label,
            },
        )
        _add_port(bridge_node, "in", _port("in", "bridge_in"))
        _add_port(bridge_node, "out", _port("out", "bridge_out"))
        _add_edge(
            edges,
            edge_keys,
            prev_node,
            prev_port,
            bridge_id,
            "in",
            "bridgeflow",
            metadata=metadata,
        )
        prev_node = bridge_id
        prev_port = "out"

    _add_edge(
        edges,
        edge_keys,
        prev_node,
        prev_port,
        target_node_id,
        target_port_id,
        role,
        metadata=metadata,
    )


# ---------------------------------------------------------------------------
# Module specs from task YAML (V1 analyze_config parity without live modules)
# ---------------------------------------------------------------------------


def _sampling_parameter_names(config: Mapping[str, Any]) -> list[str]:
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    variables = sampling.get("Variables") if isinstance(sampling, Mapping) else None
    names: list[str] = []
    if isinstance(variables, list):
        for item in variables:
            if isinstance(item, Mapping) and item.get("name"):
                names.append(str(item["name"]))
            elif isinstance(item, str):
                names.append(str(item))
    return _unique_names(names)


def _selection_spec(module: Mapping[str, Any]) -> dict[str, Any] | None:
    expression = module.get("selection")
    if expression is None or str(expression).strip() == "":
        return None
    text = str(expression)
    return {
        "expression": text,
        "variables": _expression_free_names(text),
    }


def _calculator_io_lists(module: Mapping[str, Any]) -> tuple[list[Any], list[Any]]:
    """Resolve calculator input/output lists from V1 execution.layout or flat keys."""
    execution = module.get("execution") if isinstance(module.get("execution"), Mapping) else {}
    raw_inputs = list(execution.get("input") or module.get("input") or [])
    raw_outputs = list(execution.get("output") or module.get("output") or [])
    return raw_inputs, raw_outputs


def _build_io_specs(
    raw_items: Sequence[Any],
    *,
    module_type: str,
    direction: str,
) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    for index, item in enumerate(raw_items or []):
        if isinstance(item, str):
            name = item
            is_file = False
            source_or_produced = [item]
            metadata: dict[str, Any] = {}
        elif isinstance(item, Mapping):
            name = str(item.get("name", f"{direction}_{index}"))
            is_file = _spec_is_file(item, module_type, direction=direction)
            if direction == "input":
                source_or_produced = _input_source_names(item)
            else:
                source_or_produced = _output_produced_names(item)
            metadata = {}
            file_type = str(item.get("type", "") or "").strip()
            if file_type:
                metadata["file_type"] = file_type
            if item.get("entry"):
                metadata["entry"] = str(item["entry"])
            if item.get("expression"):
                metadata["expression"] = str(item["expression"])
        else:
            continue
        entry: dict[str, Any] = {
            "id": str(name),
            "label": str(name),
            "kind": "file" if is_file else "variable",
            "metadata": metadata,
        }
        if direction == "input":
            entry["source_names"] = _unique_names(source_or_produced)
        else:
            entry["produced_names"] = _unique_names(source_or_produced)
        specs.append(entry)
    return specs


def _module_specs_from_yaml(
    module: Mapping[str, Any],
    *,
    module_type: str,
) -> dict[str, Any]:
    if module_type == "Calculator":
        raw_inputs, raw_outputs = _calculator_io_lists(module)
    else:
        raw_inputs = list(module.get("input") or module.get("inputs") or [])
        raw_outputs = list(module.get("output") or module.get("outputs") or [])
        # Operas may use mapping form for inputs/outputs.
        if isinstance(module.get("inputs"), Mapping):
            raw_inputs = list(module["inputs"].keys())
        if isinstance(module.get("outputs"), Mapping):
            raw_outputs = list(module["outputs"].keys())

    required = [
        str(val).strip()
        for val in (module.get("required_modules") or [])
        if str(val).strip()
    ]
    if not required and module_type in {"Calculator", "Operas", "Likelihood"}:
        required = ["Parameters"]

    kind, role = _module_role(module_type)
    return {
        "name": str(module.get("name") or module_type),
        "type": module_type,
        "kind": kind,
        "role": role,
        "inputs": _build_io_specs(raw_inputs, module_type=module_type, direction="input"),
        "outputs": _build_io_specs(raw_outputs, module_type=module_type, direction="output"),
        "selection": _selection_spec(module),
        "required_modules": required,
        "operator": str(module.get("operator") or "") if module.get("operator") else None,
        "call_mode": str(module.get("call_mode") or "") if module.get("call_mode") else None,
    }


def _likelihood_module_from_config(config: Mapping[str, Any]) -> dict[str, Any] | None:
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    entries = []
    if isinstance(sampling, Mapping):
        raw = sampling.get("LogLikelihood") or sampling.get("loglikelihood") or []
        if isinstance(raw, list):
            entries = [item for item in raw if isinstance(item, Mapping)]
    if not entries:
        return None
    inputs: list[dict[str, Any]] = []
    outputs: list[dict[str, Any]] = []
    source_names: list[str] = []
    for item in entries:
        name = str(item.get("name") or "LogL")
        expr = str(item.get("expression") or "")
        free = _expression_free_names(expr) if expr else []
        source_names.extend(free)
        outputs.append(
            {
                "id": name,
                "label": name,
                "kind": "variable",
                "produced_names": [name],
                "metadata": {"expression": expr} if expr else {},
            }
        )
    for var in _unique_names(source_names):
        inputs.append(
            {
                "id": var,
                "label": var,
                "kind": "variable",
                "source_names": [var],
                "metadata": {},
            }
        )
    return {
        "name": "LogLikelihood",
        "type": "Likelihood",
        "kind": "module",
        "role": "likelihood",
        "inputs": inputs,
        "outputs": outputs,
        "selection": None,
        "required_modules": ["Parameters"],
        "operator": None,
        "call_mode": None,
    }


def _collect_modules_and_layers(
    config: Mapping[str, Any],
    *,
    include_likelihood: bool,
) -> tuple[dict[str, dict[str, Any]], dict[int, list[str]]]:
    """Return module_specs_by_name and layer_id → module name list (1-based layers)."""
    calculators = (config.get("Calculators") or {}).get("Modules") or []
    operas = (config.get("Operas") or {}).get("Modules") or []
    # Match the execution plan: a multi-mode calculator is one physical
    # PackID pool but several logical workflow steps with independent I/O.
    calc_list = expand_calculator_modes(
        [dict(item) for item in calculators if isinstance(item, Mapping)]
    )
    opera_list = [dict(item) for item in operas if isinstance(item, Mapping)]

    calc_layers = resolve_module_layers(calc_list)
    opera_layers = resolve_module_layers(opera_list)
    max_calc = max(calc_layers.values()) if calc_layers else -1
    opera_base = max_calc + 1 if calc_layers else 0

    specs: dict[str, dict[str, Any]] = {}
    layer_map: dict[int, list[str]] = {}

    # V1: Parameters occupy layer 1; first workflow modules start at layer 2.
    for module in calc_list:
        name = str(module.get("name") or "").strip()
        if not name:
            continue
        plan_layer = int(calc_layers.get(name, 0))
        layer_id = plan_layer + 2
        spec = _module_specs_from_yaml(module, module_type="Calculator")
        shared_mode = mode_info(module)
        if shared_mode is not None:
            parent, mode = shared_mode
            # Same-layer sibling modes are logically independent but cannot
            # execute concurrently: they acquire one shared physical pool.
            spec["shared_runtime"] = {
                "parent": parent,
                "mode": mode,
                "pool": parent,
                "execution_constraint": "serialized_with_sibling_modes",
            }
            spec["display_label"] = (
                f"{name}\n(shared PackID: {parent}; serialized)"
            )
        specs[name] = spec
        layer_map.setdefault(layer_id, []).append(name)

    for module in opera_list:
        name = str(module.get("name") or "").strip()
        if not name:
            continue
        plan_layer = int(opera_layers.get(name, 0))
        layer_id = opera_base + plan_layer + 2
        specs[name] = _module_specs_from_yaml(module, module_type="Operas")
        layer_map.setdefault(layer_id, []).append(name)

    last_module_layer = max(layer_map.keys()) if layer_map else 1

    # Nuisance optimize step (D13.4) — drawn when Sampling.Nuisance / Nuisance present.
    from jarvishep2.Module.nuisance import extract_nuisance_config

    nuisance_block = extract_nuisance_config(config)
    if isinstance(nuisance_block, Mapping) and nuisance_block:
        nuis_layer = last_module_layer + 1
        method = str(nuisance_block.get("Method") or "Profile1D")
        specs["NuisanceOptimize"] = {
            "name": "NuisanceOptimize",
            "type": "Nuisance",
            "kind": "module",
            "role": "nuisance",
            "inputs": [],
            "outputs": [
                {
                    "id": "NuisanceLogL",
                    "label": "NuisanceLogL",
                    "kind": "variable",
                    "produced_names": ["NuisanceLogL"],
                    "metadata": {"method": method},
                }
            ],
            "selection": None,
            "required_modules": ["Parameters"],
            "operator": None,
            "call_mode": None,
        }
        layer_map.setdefault(nuis_layer, []).append("NuisanceOptimize")
        last_module_layer = nuis_layer

    # Likelihood is optional on the flowchart (off by default).
    if include_likelihood and (calc_list or opera_list or nuisance_block):
        likelihood = _likelihood_module_from_config(config)
        if likelihood is None:
            likelihood = {
                "name": "LogLikelihood",
                "type": "Likelihood",
                "kind": "module",
                "role": "likelihood",
                "inputs": [],
                "outputs": [
                    {
                        "id": "LogL",
                        "label": "LogL",
                        "kind": "variable",
                        "produced_names": ["LogL"],
                        "metadata": {},
                    }
                ],
                "selection": None,
                "required_modules": ["Parameters"],
                "operator": None,
                "call_mode": None,
            }
        ll_layer = last_module_layer + 1
        specs["LogLikelihood"] = likelihood
        layer_map.setdefault(ll_layer, []).append("LogLikelihood")

    return specs, layer_map


# ---------------------------------------------------------------------------
# Public builders
# ---------------------------------------------------------------------------


def build_flowchart_scene_from_config(
    config: Mapping[str, Any],
    *,
    include_likelihood: bool = False,
) -> dict[str, Any]:
    """Build a full V1-compatible flowchart scene from a task config mapping.

    Likelihood is **not** drawn by default (``include_likelihood=False``): it is
    a scoring step, not a workflow module users usually want on the flowchart.
    """
    module_specs, module_layers = _collect_modules_and_layers(
        config, include_likelihood=include_likelihood
    )

    nodes_by_id: dict[str, dict[str, Any]] = {}
    layer_nodes: dict[int, list[str]] = {}
    edges: list[dict[str, Any]] = []
    edge_keys: set[tuple[Any, ...]] = set()
    variable_nodes: dict[str, str] = {}

    file_input_names = {
        source_name
        for specs in module_specs.values()
        for spec in specs.get("inputs") or []
        if spec.get("kind") == "file"
        for source_name in (spec.get("source_names") or [])
    }

    def ensure_variable_node(variable_name: str, layer: int, role: str = "observable", metadata=None):
        node_id = f"var::{variable_name}"
        _add_node(
            nodes_by_id,
            layer_nodes,
            node_id,
            kind="variable",
            role=role,
            layer=layer,
            label=variable_name,
            metadata=metadata or {},
        )
        node = nodes_by_id[node_id]
        _add_port(node, "in", _port("in", role))
        _add_port(node, "out", _port("out", role))
        variable_nodes[variable_name] = node_id
        return node_id

    def get_variable_node(variable_name: str, layer: int, role: str = "observable", metadata=None):
        node_id = variable_nodes.get(variable_name)
        if node_id is not None:
            node = nodes_by_id[node_id]
            if layer < int(node["layer_index"]):
                node["layer_index"] = int(layer)
            return node_id
        return ensure_variable_node(variable_name, layer, role=role, metadata=metadata)

    def get_producer_layer(variable_name: str, fallback_layer: int) -> int:
        node_id = variable_nodes.get(variable_name)
        if node_id and node_id in nodes_by_id:
            return int(nodes_by_id[node_id]["layer_index"])
        return int(fallback_layer)

    scan = config.get("Scan") if isinstance(config.get("Scan"), Mapping) else {}
    scan_name = str((scan or {}).get("name") or config.get("scan_name") or "scan")
    project_name = str(config.get("project_name") or "")
    workflow_label = str(config.get("workflow_name") or scan_name or "Workflow")

    scene_metadata: dict[str, Any] = {
        "producer": "Jarvis-HEP",
        "producer_version": get_runtime_version(),
        "workflow_name": workflow_label,
        "scan_name": scan_name,
    }
    if project_name:
        scene_metadata["project_name"] = project_name

    # --- Parameters (layer 1) ---
    param_layer = 1
    param_names = _sampling_parameter_names(config)
    param_node = _add_node(
        nodes_by_id,
        layer_nodes,
        "Parameters",
        kind="source",
        role="parameter_source",
        layer=param_layer,
        label="Parameters",
        metadata={"module_type": "Parameter", "required_modules": []},
    )
    for port_name in param_names:
        _add_port(param_node, "out", _port(port_name, "parameter", label=port_name))
        variable_node = get_variable_node(
            port_name,
            param_layer,
            role="parameter",
            metadata={"origin": "Parameters"},
        )
        _add_edge(
            edges,
            edge_keys,
            "Parameters",
            port_name,
            variable_node,
            "in",
            "parameterflow",
        )

    # --- Modules by layer ---
    for layer_id in sorted(module_layers.keys()):
        for module_name in module_layers[layer_id]:
            specs = module_specs[module_name]
            module_node = _add_node(
                nodes_by_id,
                layer_nodes,
                module_name,
                kind=specs["kind"],
                role=specs["role"],
                layer=layer_id,
                label=str(specs.get("display_label") or module_name),
                metadata={
                    "module_type": specs["type"],
                    "required_modules": list(specs.get("required_modules") or []),
                },
            )
            if specs.get("operator"):
                module_node.setdefault("metadata", {})["operator"] = specs["operator"]
            if specs.get("call_mode"):
                module_node.setdefault("metadata", {})["call_mode"] = specs["call_mode"]
            if specs.get("shared_runtime"):
                module_node.setdefault("metadata", {})["shared_runtime"] = specs["shared_runtime"]

            selection = specs.get("selection")
            if selection:
                module_node["selection"] = selection
                module_node.setdefault("metadata", {})["selection"] = selection
                for selection_var in selection.get("variables") or []:
                    selection_port = f"selection::{selection_var}"
                    _add_port(
                        module_node,
                        "in",
                        _port(
                            selection_port,
                            "selection",
                            label=selection_var,
                            metadata={"selection": True},
                        ),
                    )
                    source_node = get_variable_node(
                        selection_var,
                        get_producer_layer(selection_var, layer_id),
                        role="observable",
                        metadata={"selection_consumer": module_name},
                    )
                    _bridge_chain(
                        nodes_by_id,
                        layer_nodes,
                        edges,
                        edge_keys,
                        source_node,
                        "out",
                        nodes_by_id[source_node]["layer_index"],
                        module_name,
                        selection_port,
                        layer_id,
                        "selectionflow",
                        selection_var,
                        metadata={
                            "module": module_name,
                            "selection": True,
                            "variable": selection_var,
                        },
                    )

            for spec in specs.get("inputs") or []:
                _add_port(
                    module_node,
                    "in",
                    _port(
                        spec["id"],
                        spec["kind"],
                        label=spec.get("label"),
                        metadata=spec.get("metadata"),
                    ),
                )
                if spec["kind"] == "file":
                    file_id = f"file::{module_name}::input::{spec['id']}"
                    file_node = _add_node(
                        nodes_by_id,
                        layer_nodes,
                        file_id,
                        kind="file",
                        role="input_file",
                        layer=layer_id,
                        label=str(spec["label"]),
                        metadata={
                            "direction": "input",
                            "file_type": (spec.get("metadata") or {}).get("file_type"),
                        },
                    )
                    _add_port(file_node, "in", _port("in", "file_input"))
                    _add_port(file_node, "out", _port("out", "file_output"))
                    for source_name in spec.get("source_names") or []:
                        source_node = get_variable_node(
                            source_name,
                            get_producer_layer(source_name, layer_id),
                            role="observable",
                            metadata={"consumer": module_name},
                        )
                        _add_edge(
                            edges,
                            edge_keys,
                            source_node,
                            "out",
                            file_id,
                            "in",
                            "fileflow",
                            metadata={
                                "variable": source_name,
                                "module": module_name,
                                "direction": "input",
                            },
                        )
                    _bridge_chain(
                        nodes_by_id,
                        layer_nodes,
                        edges,
                        edge_keys,
                        file_id,
                        "out",
                        layer_id,
                        module_name,
                        spec["id"],
                        layer_id,
                        "fileflow",
                        spec["id"],
                        metadata={
                            "module": module_name,
                            "input": spec["id"],
                            "direction": "input",
                        },
                    )
                else:
                    for source_name in spec.get("source_names") or []:
                        source_node = get_variable_node(
                            source_name,
                            get_producer_layer(source_name, layer_id),
                            role="observable",
                            metadata={"consumer": module_name},
                        )
                        _bridge_chain(
                            nodes_by_id,
                            layer_nodes,
                            edges,
                            edge_keys,
                            source_node,
                            "out",
                            nodes_by_id[source_node]["layer_index"],
                            module_name,
                            spec["id"],
                            layer_id,
                            "dataflow",
                            source_name,
                            metadata={"module": module_name, "input": spec["id"]},
                        )

            if specs["type"] == "Parameter":
                continue

            for spec in specs.get("outputs") or []:
                _add_port(
                    module_node,
                    "out",
                    _port(
                        spec["id"],
                        spec["kind"],
                        label=spec.get("label"),
                        metadata=spec.get("metadata"),
                    ),
                )
                if spec["kind"] == "file":
                    file_id = f"file::{module_name}::output::{spec['id']}"
                    file_node = _add_node(
                        nodes_by_id,
                        layer_nodes,
                        file_id,
                        kind="file",
                        role="output_file",
                        layer=layer_id,
                        label=str(spec["label"]),
                        metadata={
                            "direction": "output",
                            "file_type": (spec.get("metadata") or {}).get("file_type"),
                        },
                    )
                    _add_port(file_node, "in", _port("in", "file_input"))
                    _add_port(file_node, "out", _port("out", "file_output"))
                    _add_edge(
                        edges,
                        edge_keys,
                        module_name,
                        spec["id"],
                        file_id,
                        "in",
                        "fileflow",
                        metadata={
                            "module": module_name,
                            "output": spec["id"],
                            "direction": "output",
                        },
                    )
                    produced_names = list(spec.get("produced_names") or [])
                    if spec["id"] in file_input_names:
                        produced_names.append(spec["id"])
                    for target_name in _unique_names(produced_names):
                        target_node = get_variable_node(
                            target_name,
                            layer_id,
                            role="observable",
                            metadata={"producer": module_name},
                        )
                        _add_edge(
                            edges,
                            edge_keys,
                            file_id,
                            "out",
                            target_node,
                            "in",
                            "fileflow",
                            metadata={
                                "module": module_name,
                                "variable": target_name,
                                "direction": "output",
                            },
                        )
                else:
                    for target_name in spec.get("produced_names") or []:
                        target_node = get_variable_node(
                            target_name,
                            layer_id,
                            role="observable",
                            metadata={"producer": module_name},
                        )
                        _bridge_chain(
                            nodes_by_id,
                            layer_nodes,
                            edges,
                            edge_keys,
                            module_name,
                            spec["id"],
                            layer_id,
                            target_node,
                            "in",
                            layer_id,
                            "dataflow",
                            target_name,
                            metadata={"module": module_name, "output": spec["id"]},
                        )

    nodes: list[dict[str, Any]] = []
    for node in nodes_by_id.values():
        node.pop("_port_keys", None)
        layer_index = int(node.pop("layer_index"))
        node["layer"] = _layer_name(layer_index)
        nodes.append(node)

    layers: list[dict[str, Any]] = []
    for layer_id in sorted(layer_nodes.keys()):
        labels = [
            nodes_by_id[node_id]["label"]
            for node_id in layer_nodes[layer_id]
            if nodes_by_id[node_id]["kind"] in {"module", "source"}
        ]
        layer_name = _layer_name(layer_id)
        if layer_id == 1:
            label = "Parameters"
        elif labels:
            label = ", ".join(labels)
        else:
            label = f"Layer {layer_id}"
        layers.append(
            {
                "id": layer_name,
                "index": int(layer_id),
                "label": label,
                "nodes": [str(node_id) for node_id in layer_nodes[layer_id]],
            }
        )

    # A shared-pack calculator has several logical data-flow steps but only
    # one physical runtime pool.  Keep this topology explicit instead of
    # making a renderer infer it from dotted node names.
    shared_runtime_groups: dict[str, list[str]] = {}
    for module_name, specs in module_specs.items():
        runtime = specs.get("shared_runtime")
        if not isinstance(runtime, Mapping):
            continue
        parent = str(runtime.get("parent") or "").strip()
        if parent:
            shared_runtime_groups.setdefault(parent, []).append(module_name)
    groups = [
        {
            "id": f"shared_runtime::{parent}",
            "kind": "shared_runtime",
            "label": f"{parent} · shared PackID pool · modes serialize",
            "nodes": node_names,
            "metadata": {
                "parent": parent,
                "pool": parent,
                "execution_constraint": "serialized_with_sibling_modes",
            },
        }
        for parent, node_names in shared_runtime_groups.items()
    ]

    return {
        "schema": FLOWCHART_SCHEMA,
        "scene_type": FLOWCHART_SCENE_TYPE,
        "scene_id": FLOWCHART_SCENE_ID,
        "metadata": scene_metadata,
        "layers": layers,
        "nodes": nodes,
        "edges": edges,
        "groups": groups,
    }


def build_flowchart_scene(
    steps: Sequence[ExecutionStep] | Sequence[Mapping[str, Any]],
    *,
    scan_name: str = "scan",
    project_name: str = "",
    producer: str = "Jarvis-HEP",
    config: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a flowchart scene.

    Prefer ``build_flowchart_scene_from_config`` for full V1 IO-port graphs.
    When only an execution plan is available, emit a simplified layered graph
    (still ``jarvisplot.scene/v1``).
    """
    if config is not None:
        scene = build_flowchart_scene_from_config(config)
        if scan_name:
            scene.setdefault("metadata", {})["scan_name"] = str(scan_name)
        if project_name:
            scene.setdefault("metadata", {})["project_name"] = str(project_name)
        return scene

    # Plan-only fallback (no YAML IO): Parameters → layered modules.
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

    # Synthesize a minimal config-like structure for the full exporter when possible.
    calc = []
    opera = []
    for step in normalized:
        entry = {"name": step.name, "required_modules": []}
        if step.type == "calculator":
            calc.append(entry)
        elif step.type == "opera":
            opera.append(entry)
    if calc or opera:
        synthetic = {
            "Scan": {"name": scan_name},
            "project_name": project_name,
            "Calculators": {"Modules": calc} if calc else {},
            "Operas": {"Modules": opera} if opera else {},
        }
        # Plan-only fallback also omits Likelihood nodes by default.
        scene = build_flowchart_scene_from_config(
            synthetic, include_likelihood=False
        )
        scene.setdefault("metadata", {})["producer"] = producer
        return scene

    # Empty plan
    return {
        "schema": FLOWCHART_SCHEMA,
        "scene_type": FLOWCHART_SCENE_TYPE,
        "scene_id": FLOWCHART_SCENE_ID,
        "metadata": {
            "producer": producer,
            "producer_version": get_runtime_version(),
            "workflow_name": scan_name,
            "scan_name": scan_name,
            "project_name": project_name,
        },
        "layers": [
            {
                "id": "layer_1",
                "index": 1,
                "label": "Parameters",
                "nodes": ["Parameters"],
            }
        ],
        "nodes": [
            {
                "id": "Parameters",
                "kind": "source",
                "role": "parameter_source",
                "label": "Parameters",
                "layer": "layer_1",
                "in_ports": [],
                "out_ports": [],
                "metadata": {"module_type": "Parameter", "required_modules": []},
            }
        ],
        "edges": [],
    }


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
    """Render via stock JarvisPLOT (no V2 adapter). Returns PNG path or None."""
    try:
        from jarvisplot import render_flowchart
    except ImportError:
        return None
    import logging
    import warnings

    out = os.path.abspath(str(output_path))
    parent = os.path.dirname(out)
    if parent:
        os.makedirs(parent, exist_ok=True)
    # Defense in depth: mute matplotlib font-manager noise even if an older
    # JarvisPLOT is installed (cache rebuild / missing weight warnings).
    prev_mpl = logging.getLogger("matplotlib").level
    prev_fm = logging.getLogger("matplotlib.font_manager").level
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r".*[Ff]ont.*")
        warnings.filterwarnings("ignore", message=r".*findfont.*")
        try:
            rendered = render_flowchart(dict(scene), output_path=out)
        finally:
            logging.getLogger("matplotlib").setLevel(prev_mpl)
            logging.getLogger("matplotlib.font_manager").setLevel(prev_fm)
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
