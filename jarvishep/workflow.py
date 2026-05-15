#!/usr/bin/env python3
from __future__ import annotations

import json
import networkx as nx 
from jarvishep.base import Base
from jarvishep.Module.module import Module
from jarvishep.versioning import get_runtime_version
import numpy as np 
import contextlib
import os, io 
import asyncio
# from numpy.lib._type_check_impl import imag

class Workflow(Base):
    FLOWCHART_SCHEMA = "jarvisplot.scene/v1"
    FLOWCHART_SCENE_TYPE = "flowchart"
    FLOWCHART_SCENE_ID = "workflow_main"

    def __init__(self):
        super().__init__()
        self.modules = {}
        self.module_states = {}
        self.module_dependencies = {}
        self.parameter_module = None
        self.library_modules = {}
        self.graph = {}
        self.layer_info = {}
        self.layers = []
        self.workflow = {}

    @staticmethod
    def _flowchart_jsonable(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, dict):
            return {str(key): Workflow._flowchart_jsonable(val) for key, val in value.items()}
        if isinstance(value, (list, tuple)):
            return [Workflow._flowchart_jsonable(item) for item in value]
        if isinstance(value, set):
            return [Workflow._flowchart_jsonable(item) for item in sorted(value, key=str)]
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        return str(value)

    @staticmethod
    def _flowchart_layer_name(layer_index: int) -> str:
        return f"layer_{int(layer_index)}"

    @staticmethod
    def _flowchart_edge_role(role: str) -> str:
        normalized = str(role).strip().lower()
        allowed = {"parameterflow", "dataflow", "fileflow", "bridgeflow", "selectionflow"}
        return normalized if normalized in allowed else "dataflow"

    @staticmethod
    def _flowchart_module_role(module) -> tuple[str, str]:
        module_type = str(getattr(module, "type", "Module") or "Module")
        normalized = module_type.strip().lower().replace(" ", "_")
        if normalized == "parameter":
            return "source", "parameter_source"
        if normalized == "operas":
            return "module", "operas"
        if normalized == "calculator":
            return "module", "calculator"
        return "module", normalized or "module"

    @staticmethod
    def _flowchart_selection_spec(module):
        expression = getattr(module, "selection", None)
        if expression is None or str(expression).strip() == "":
            return None

        deps = getattr(module, "_selection_deps", None)
        if deps is None and hasattr(module, "_compile_selection"):
            module._compile_selection()
            deps = getattr(module, "_selection_deps", ())

        return {
            "expression": str(expression),
            "variables": Workflow._flowchart_unique_names([str(dep) for dep in (deps or [])]),
        }

    @staticmethod
    def _flowchart_port(port_id, role, label=None, metadata=None):
        port = {"id": str(port_id), "role": str(role)}
        if label is not None and str(label) != str(port_id):
            port["label"] = str(label)
        if metadata:
            port["metadata"] = Workflow._flowchart_jsonable(metadata)
        return port

    @staticmethod
    def _flowchart_unique_names(values):
        seen = set()
        ordered = []
        for value in values or []:
            text = str(value)
            if text in seen:
                continue
            seen.add(text)
            ordered.append(text)
        return ordered

    @staticmethod
    def _flowchart_input_names(spec):
        source_names = []
        if isinstance(spec, str):
            return [str(spec)]
        if not isinstance(spec, dict):
            return [str(spec)]

        variables = spec.get("variables")
        if isinstance(variables, dict):
            for item in variables.values():
                if isinstance(item, dict) and "inc" in item:
                    source_names.extend([str(val) for val in item.get("inc", [])])
                elif isinstance(item, dict) and item.get("name"):
                    source_names.append(str(item["name"]))
                elif isinstance(item, str):
                    source_names.append(str(item))
        elif isinstance(variables, list):
            for item in variables:
                if isinstance(item, dict) and item.get("name"):
                    source_names.append(str(item["name"]))
                elif isinstance(item, str):
                    source_names.append(str(item))

        if not source_names and spec.get("_inc"):
            source_names.extend([str(val) for val in spec.get("_inc", [])])
        if not source_names and spec.get("entry"):
            source_names.append(str(spec["entry"]))
        if not source_names and spec.get("name"):
            source_names.append(str(spec["name"]))
        return source_names

    @staticmethod
    def _flowchart_output_names(spec):
        output_names = []
        if isinstance(spec, str):
            return [str(spec)]
        if not isinstance(spec, dict):
            return [str(spec)]

        variables = spec.get("variables")
        if isinstance(variables, dict):
            output_names.extend([str(key) for key in variables.keys()])
        elif isinstance(variables, list):
            for item in variables:
                if isinstance(item, dict) and item.get("name"):
                    output_names.append(str(item["name"]))
                elif isinstance(item, str):
                    output_names.append(str(item))

        if not output_names and spec.get("name"):
            output_names.append(str(spec["name"]))
        return output_names

    @staticmethod
    def _flowchart_spec_is_file(spec, module_type, *, direction):
        if isinstance(spec, str):
            return False
        if not isinstance(spec, dict):
            return False

        spec_type = str(spec.get("type", "") or "").strip().lower()
        if spec.get("path"):
            return True
        if spec_type in {"file", "slha", "json", "xslha"}:
            return True

        if direction == "input":
            return module_type == "Calculator" and bool(spec.get("variables")) and not spec.get("expression")

        return bool(spec.get("variables"))

    def _flowchart_module_specs(self, module):
        raw_inputs = list(getattr(module, "input", []) or [])
        raw_outputs = list(getattr(module, "output", []) or [])

        if not raw_inputs and getattr(module, "inputs", None):
            raw_inputs = list(getattr(module, "inputs", {}).keys())
        if not raw_outputs and getattr(module, "outputs", None):
            raw_outputs = list(getattr(module, "outputs", {}).keys())

        module_type = str(getattr(module, "type", "Module") or "Module")
        inputs = []
        for index, item in enumerate(raw_inputs):
            name = item if isinstance(item, str) else str(item.get("name", f"input_{index}"))
            is_file = self._flowchart_spec_is_file(item, module_type, direction="input")
            source_names = self._flowchart_input_names(item)
            metadata = {}
            if isinstance(item, dict):
                file_type = str(item.get("type", "") or "").strip()
                if file_type:
                    metadata["file_type"] = file_type
                if item.get("entry"):
                    metadata["entry"] = str(item["entry"])
                if item.get("expression"):
                    metadata["expression"] = str(item["expression"])
            inputs.append(
                {
                    "id": str(name),
                    "label": str(name),
                    "kind": "file" if is_file else "variable",
                    "source_names": self._flowchart_unique_names(source_names),
                    "metadata": metadata,
                }
            )

        outputs = []
        for index, item in enumerate(raw_outputs):
            name = item if isinstance(item, str) else str(item.get("name", f"output_{index}"))
            is_file = self._flowchart_spec_is_file(item, module_type, direction="output")
            produced_names = self._flowchart_output_names(item)
            metadata = {}
            if isinstance(item, dict):
                file_type = str(item.get("type", "") or "").strip()
                if file_type:
                    metadata["file_type"] = file_type
                if item.get("entry"):
                    metadata["entry"] = str(item["entry"])
            outputs.append(
                {
                    "id": str(name),
                    "label": str(name),
                    "kind": "file" if is_file else "variable",
                    "produced_names": self._flowchart_unique_names(produced_names),
                    "metadata": metadata,
                }
            )

        return {
            "type": module_type,
            "role": self._flowchart_module_role(module)[1],
            "inputs": inputs,
            "outputs": outputs,
            "selection": self._flowchart_selection_spec(module),
            "required_modules": [str(val) for val in getattr(module, "required_modules", []) or []],
        }

    def _flowchart_resolve_layers(self):
        if not getattr(self, "calc_layer", None):
            self.resolve_dependencies()

        layer_map = {}
        for layer_id, layer in (self.calc_layer or {}).items():
            module_names = layer.get("module", []) if isinstance(layer, dict) else []
            layer_map[layer_id] = [str(name) for name in module_names]
        return layer_map

    def _flowchart_add_node(self, nodes_by_id, layer_nodes, node_id, *, kind, role, layer, label, metadata=None):
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
                node["metadata"] = self._flowchart_jsonable(metadata)
            nodes_by_id[node_id] = node
            layer_nodes.setdefault(int(layer), []).append(node_id)
            return node

        if int(layer) < int(node.get("layer_index", layer)):
            node["layer_index"] = int(layer)
        if metadata:
            existing_metadata = node.get("metadata", {})
            if not isinstance(existing_metadata, dict):
                existing_metadata = {}
            merged = dict(existing_metadata)
            merged.update(self._flowchart_jsonable(metadata))
            node["metadata"] = merged
        return node

    def _flowchart_add_port(self, node, port_direction, port):
        key = f"{port_direction}:{port['id']}"
        port_keys = node.setdefault("_port_keys", set())
        if key in port_keys:
            return
        node[f"{port_direction}_ports"].append(port)
        port_keys.add(key)

    def _flowchart_add_edge(self, edges, edge_keys, source_node, source_port, target_node, target_port, role, metadata=None):
        role = self._flowchart_edge_role(role)
        key = (source_node, source_port, target_node, target_port, role)
        if key in edge_keys:
            return
        edge = {
            "source": {"node": str(source_node), "port": str(source_port)},
            "target": {"node": str(target_node), "port": str(target_port)},
            "role": str(role),
        }
        if metadata:
            edge["metadata"] = self._flowchart_jsonable(metadata)
        edges.append(edge)
        edge_keys.add(key)

    def _flowchart_bridge_chain(
        self,
        nodes_by_id,
        layer_nodes,
        edges,
        edge_keys,
        source_node_id,
        source_port_id,
        source_layer,
        target_node_id,
        target_port_id,
        target_layer,
        role,
        bridge_label,
        metadata=None,
    ):
        if target_layer - source_layer <= 1:
            self._flowchart_add_edge(
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
            bridge_node = self._flowchart_add_node(
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
            self._flowchart_add_port(bridge_node, "in", self._flowchart_port("in", "bridge_in"))
            self._flowchart_add_port(bridge_node, "out", self._flowchart_port("out", "bridge_out"))
            self._flowchart_add_edge(
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

        self._flowchart_add_edge(
            edges,
            edge_keys,
            prev_node,
            prev_port,
            target_node_id,
            target_port_id,
            role,
            metadata=metadata,
        )

    def build_flowchart_semantics(self, workflow_name: str | None = None):
        if not getattr(self, "calc_layer", None):
            self.resolve_dependencies()

        module_layers = self._flowchart_resolve_layers()
        nodes_by_id = {}
        layer_nodes = {}
        edges = []
        edge_keys = set()
        variable_nodes = {}

        def ensure_variable_node(variable_name, layer, role="observable", metadata=None):
            node_id = f"var::{variable_name}"
            node = self._flowchart_add_node(
                nodes_by_id,
                layer_nodes,
                node_id,
                kind="variable",
                role=role,
                layer=layer,
                label=variable_name,
                metadata=metadata or {},
            )
            self._flowchart_add_port(node, "in", self._flowchart_port("in", role))
            self._flowchart_add_port(node, "out", self._flowchart_port("out", role))
            variable_nodes[variable_name] = node_id
            return node_id

        def get_variable_node(variable_name, layer, role="observable", metadata=None):
            node_id = variable_nodes.get(variable_name)
            if node_id is not None:
                node = nodes_by_id[node_id]
                if layer < node["layer_index"]:
                    node["layer_index"] = int(layer)
                return node_id
            return ensure_variable_node(variable_name, layer, role=role, metadata=metadata)

        def get_producer_layer(variable_name, fallback_layer):
            for node in nodes_by_id.values():
                if node["kind"] != "variable":
                    continue
                if node["id"] == f"var::{variable_name}":
                    return int(node["layer_index"])
            return int(fallback_layer)

        workflow_label = workflow_name or getattr(self, "workflow_name", None) or "Workflow"
        scene_metadata = {
            "producer": "Jarvis-HEP",
            "producer_version": get_runtime_version(),
            "workflow_name": workflow_label,
        }
        if isinstance(getattr(self, "info", None), dict):
            project_name = self.info.get("project_name")
            scan_name = self.info.get("scan_name")
            if project_name:
                scene_metadata["project_name"] = str(project_name)
            if scan_name:
                scene_metadata["scan_name"] = str(scan_name)
        parameter_module = getattr(self, "parameter_module", None) or self.modules.get("Parameters")
        if parameter_module is not None:
            param_layer = 1
            param_kind, param_role = self._flowchart_module_role(parameter_module)
            param_node = self._flowchart_add_node(
                nodes_by_id,
                layer_nodes,
                "Parameters",
                kind=param_kind,
                role=param_role,
                layer=param_layer,
                label="Parameters",
                metadata={
                    "module_type": getattr(parameter_module, "type", "Parameter"),
                    "required_modules": [str(val) for val in getattr(parameter_module, "required_modules", []) or []],
                },
            )
            for port_name, port_role in (getattr(parameter_module, "outputs", {}) or {}).items():
                normalized_port_role = str(port_role or "parameter").strip().lower()
                if normalized_port_role == "param":
                    normalized_port_role = "parameter"
                elif "nuisa" in normalized_port_role:
                    normalized_port_role = "nuisance"
                elif not normalized_port_role:
                    normalized_port_role = "observable"
                self._flowchart_add_port(
                    param_node,
                    "out",
                    self._flowchart_port(port_name, normalized_port_role, label=port_name),
                )
                variable_node = get_variable_node(
                    str(port_name),
                    param_layer,
                    role=(
                        "parameter"
                        if normalized_port_role == "parameter"
                        else "nuisance"
                        if normalized_port_role == "nuisance"
                        else "observable"
                    ),
                    metadata={"origin": "Parameters"},
                )
                self._flowchart_add_edge(
                    edges,
                    edge_keys,
                    "Parameters",
                    str(port_name),
                    variable_node,
                    "in",
                    "parameterflow",
                )

        for layer_id in sorted(module_layers.keys()):
            for module_name in module_layers[layer_id]:
                module = self.modules[module_name]
                module_kind, module_role = self._flowchart_module_role(module)
                module_specs = self._flowchart_module_specs(module)
                is_parameter_source = module_specs["type"] == "Parameter" or module_role == "parameter_source"
                module_node = self._flowchart_add_node(
                    nodes_by_id,
                    layer_nodes,
                    module_name,
                    kind=module_kind,
                    role=module_role,
                    layer=layer_id,
                    label=module_name,
                    metadata={
                        "module_type": module_specs["type"],
                        "required_modules": module_specs["required_modules"],
                    },
                )
                if module_specs["type"] == "Operas" and hasattr(module, "operator"):
                    module_node.setdefault("metadata", {})["operator"] = str(getattr(module, "operator", ""))
                if hasattr(module, "call_mode"):
                    module_node.setdefault("metadata", {})["call_mode"] = str(getattr(module, "call_mode", ""))
                if module_specs["selection"]:
                    module_node["selection"] = module_specs["selection"]
                    module_node.setdefault("metadata", {})["selection"] = module_specs["selection"]

                    for selection_var in module_specs["selection"]["variables"]:
                        selection_port = f"selection::{selection_var}"
                        self._flowchart_add_port(
                            module_node,
                            "in",
                            self._flowchart_port(
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
                        self._flowchart_bridge_chain(
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

                for spec in module_specs["inputs"]:
                    self._flowchart_add_port(
                        module_node,
                        "in",
                        self._flowchart_port(spec["id"], spec["kind"], label=spec["label"], metadata=spec.get("metadata")),
                    )
                    if spec["kind"] == "file":
                        file_id = f"file::{module_name}::input::{spec['id']}"
                        file_node = self._flowchart_add_node(
                            nodes_by_id,
                            layer_nodes,
                            file_id,
                            kind="file",
                            role="input_file",
                            layer=layer_id,
                            label=str(spec["label"]),
                            metadata={"direction": "input", "file_type": spec.get("metadata", {}).get("file_type")},
                        )
                        self._flowchart_add_port(file_node, "in", self._flowchart_port("in", "file_input"))
                        self._flowchart_add_port(file_node, "out", self._flowchart_port("out", "file_output"))
                        for source_name in spec["source_names"]:
                            source_node = get_variable_node(
                                source_name,
                                get_producer_layer(source_name, layer_id),
                                role="observable",
                                metadata={"consumer": module_name},
                            )
                            self._flowchart_add_edge(
                                edges,
                                edge_keys,
                                source_node,
                                "out",
                                file_id,
                                "in",
                                "fileflow",
                                metadata={"variable": source_name, "module": module_name, "direction": "input"},
                            )
                        self._flowchart_bridge_chain(
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
                            metadata={"module": module_name, "input": spec["id"], "direction": "input"},
                        )
                    else:
                        for source_name in spec["source_names"]:
                            source_node = get_variable_node(
                                source_name,
                                get_producer_layer(source_name, layer_id),
                                role="observable",
                                metadata={"consumer": module_name},
                            )
                            self._flowchart_bridge_chain(
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

                if is_parameter_source:
                    continue

                for spec in module_specs["outputs"]:
                    self._flowchart_add_port(
                        module_node,
                        "out",
                        self._flowchart_port(spec["id"], spec["kind"], label=spec["label"], metadata=spec.get("metadata")),
                    )
                    if spec["kind"] == "file":
                        file_id = f"file::{module_name}::output::{spec['id']}"
                        file_node = self._flowchart_add_node(
                            nodes_by_id,
                            layer_nodes,
                            file_id,
                            kind="file",
                            role="output_file",
                            layer=layer_id,
                            label=str(spec["label"]),
                            metadata={"direction": "output", "file_type": spec.get("metadata", {}).get("file_type")},
                        )
                        self._flowchart_add_port(file_node, "in", self._flowchart_port("in", "file_input"))
                        self._flowchart_add_port(file_node, "out", self._flowchart_port("out", "file_output"))
                        self._flowchart_add_edge(
                            edges,
                            edge_keys,
                            module_name,
                            spec["id"],
                            file_id,
                            "in",
                            "fileflow",
                            metadata={"module": module_name, "output": spec["id"], "direction": "output"},
                        )
                        for target_name in spec["produced_names"]:
                            target_node = get_variable_node(
                                target_name,
                                layer_id,
                                role="observable",
                                metadata={"producer": module_name},
                            )
                            self._flowchart_add_edge(
                                edges,
                                edge_keys,
                                file_id,
                                "out",
                                target_node,
                                "in",
                                "fileflow",
                                metadata={"module": module_name, "variable": target_name, "direction": "output"},
                            )
                    else:
                        for target_name in spec["produced_names"]:
                            target_node = get_variable_node(
                                target_name,
                                layer_id,
                                role="observable",
                                metadata={"producer": module_name},
                            )
                            self._flowchart_bridge_chain(
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

        nodes = []
        for node in nodes_by_id.values():
            node.pop("_port_keys", None)
            layer_index = int(node.pop("layer_index"))
            node["layer"] = self._flowchart_layer_name(layer_index)
            nodes.append(node)

        layers = []
        for layer_id in sorted(layer_nodes.keys()):
            labels = [nodes_by_id[node_id]["label"] for node_id in layer_nodes[layer_id] if nodes_by_id[node_id]["kind"] == "module"]
            layer_name = self._flowchart_layer_name(layer_id)
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

        return {
            "schema": self.FLOWCHART_SCHEMA,
            "scene_type": self.FLOWCHART_SCENE_TYPE,
            "scene_id": self.FLOWCHART_SCENE_ID,
            "metadata": scene_metadata,
            "layers": layers,
            "nodes": nodes,
            "edges": edges,
        }

    async def export_flowchart_semantics(self, save_path="flowchart.json", workflow_name: str | None = None):
        semantic_graph = self.build_flowchart_semantics(workflow_name=workflow_name)
        output_path = os.path.abspath(save_path)
        output_dir = os.path.dirname(output_path) or "."
        os.makedirs(output_dir, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as handle:
            json.dump(semantic_graph, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        return semantic_graph

    def add_module(self, module):
        self.modules[module.name] = module
        self.module_states[module.name] = 'waiting'

    def add_library_module(self, module):
        self.library_modules[module.name] = module
        self.modules[module.name] = module

    def set_modules(self, modules):
        from jarvishep.Module.parameters import Parameters
        from jarvishep.Module.library import LibraryModule
        from jarvishep.Module.calculator import CalculatorModule
        from jarvishep.Module.operas import OperasModule
        
        # print(modules.keys())
        parameter_module = Parameters("Parameters", modules['Parameter'])
        if modules.get("Nuisance", False):
            parameter_module.add_nuisance(modules['Nuisance'])
        parameter_module.analyze_ios()
        self.add_module(parameter_module)
        self.parameter_module = parameter_module
        
        if "Library" in modules:
            for lib in modules['Library']:
                module = LibraryModule(
                    name=lib['name'],
                    required_modules=lib.get("required_modules", []),
                    installed=lib.get("installed", False),
                    installation=lib.get("installation", {})
                )
                self.add_library_module(module)

        for calc in modules.get('Calculator', []) or []:
            module = CalculatorModule(name=calc['name'], config=calc)
            self.add_module(module)

        for oper in modules.get("Operas", []) or []:
            module = OperasModule(name=oper["name"], config=oper)
            self.add_module(module)

    def get_workflow_dict(self):
        for kk, layer in self.calc_layer.items():
            if kk > 1: 
                module_names = layer.get("module", []) if isinstance(layer, dict) else []
                self.workflow[kk] = list(module_names)
    
    def resolve_layers(self):
        self.calc_layer = {}  
        self.libr_layer = {}

        unresolved_all = []
        unresolved_lib = []
        unresolved = []
        for name in self.modules.keys():
            if name in self.library_modules:
                unresolved_lib.append(name)
            else:
                unresolved.append(name)
            unresolved_all.append(name)

        # resolve the library layers 
        current_layer = 1
        resolved = []
        while unresolved_lib:
            new_resolved = []
            for ii in range(len(unresolved_lib)):
                if current_layer == 1:
                    if not self.modules[unresolved_lib[ii]].required_modules:
                        new_resolved.append(self.modules[unresolved_lib[ii]].name)
                else:
                    if set(self.modules[unresolved_lib[ii]].required_modules).issubset(set(resolved)):
                        new_resolved.append(self.modules[unresolved_lib[ii]].name)

            resolved.extend(new_resolved)
            unresolved_lib = [item for item in unresolved_lib if item not in resolved]
            self.libr_layer[current_layer] = new_resolved
            current_layer += 1         
        
        current_layer = 1
        while unresolved:
            # achieve the first layer module
            new_resolved = []
            for ii in range(len(unresolved)):
                if current_layer == 1:
                    if not self.modules[unresolved[ii]].required_modules:
                        new_resolved.append(self.modules[unresolved[ii]].name)
                else:
                    if set(self.modules[unresolved[ii]].required_modules).issubset(set(resolved)):
                        new_resolved.append(self.modules[unresolved[ii]].name) 
            
            resolved.extend(new_resolved)
            unresolved = [item for item in unresolved if item not in resolved]
            self.calc_layer[current_layer] = {"module": new_resolved, "ipfs": {}, "opfs": {}, "ipvs": {}, "opvs": {}, "mdls": {}}
            current_layer += 1 

    def resolve_dependencies(self):
        # Resolve dependencies based on the IOs 
        output_to_module = {} 
        for module in self.modules.values():
            for output in module.outputs:
                output_to_module[output] = module.name

        # Update required_modules via IOs 
        for module in self.modules.values():
            module.required_modules = module.required_modules or []
            input_vars = [str(input_var) for input_var in module.inputs]
            selection = self._flowchart_selection_spec(module)
            if selection:
                input_vars.extend(selection["variables"])

            for input_var in self._flowchart_unique_names(input_vars):
                if input_var in output_to_module:
                    # if the output is the input of other module, then generate teh dependencies 
                    dep_module = output_to_module[input_var]
                    if dep_module != module.name and dep_module not in module.required_modules:
                        module.required_modules.append(dep_module)

        # Resolve the layers using the current dependencies  
        self.resolve_layers()

    def analyze(self):
        # First, add edges based on module dependencies
        for module in self.modules:
            for output in module.outputs:
                for dependent_module in self.modules:
                    if output in dependent_module.inputs:
                        self.graph.add_edge(module.name, dependent_module.name)
        
        # Check for cycles which would indicate a problem in the dependency graph
        try:
            cycle = nx.find_cycle(self.graph)
            raise Exception(f"Cyclic dependency found: {cycle}")
        except nx.NetworkXNoCycle:
            pass  # All good if there are no cycles

    def run(self, sample_point):
        try:
            for module in self.modules:
                self.module_states[module.name] = 'running'
                result = module.run(sample_point)
                
                if not module.validate_result(result):
                    self.module_states[module.name] = 'failed'
                    print(f"Module {module.name} failed. Stopping workflow for this sample.")
                    # Stop the workflow for this sample point
                    return False  
                
                self.module_states[module.name] = 'completed'
            # Workflow completed successfully for this sample point    
            return True  
        except Exception as e:
            print(f"An error occurred: {e}")
            # An error occurred, stopping workflow for this sample point
            return False  

    def reset_states(self):
        for module in self.modules:
            self.module_states[module.name] = 'not started'
    
    def run_all_samples(self, sample_points):
        for sample in sample_points:
            self.reset_states()
            if not self.run(sample):
                print(f"Workflow failed for sample point {sample}. Moving to next sample.")
            else:
                print(f"Workflow completed successfully for sample point {sample}.")

    async def draw_flowchart(self, save_path="flowchart.png"):
        # with io.StringIO() as buf, contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            import matplotlib 
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt 
            import logging
            logging.getLogger('matplotlib').setLevel(logging.CRITICAL)
            logging.getLogger('PIL.PngImagePlugin').setLevel(logging.CRITICAL)
            from matplotlib.collections import LineCollection
            from jarvishep.plot import create_round_square
            from PIL import Image
            icon_dir = os.path.join(os.path.dirname(__file__), "icons")

            layerW = 6.0 
            figL = layerW * len(self.calc_layer) - 3
            figH = 4

            def resolve_module(module, lid):
                module_type = getattr(self.modules[module], "type", "Calculator")
                if lid == 1:
                    res = {"name": module, "type": module_type, "ipf": {}, "opf": {}, "ipv": {}, "opv": {}, 'mhh':0.9, 'ohh': 0., "ihh": 0., "bp": np.array([0.85, 0.]), "width": 0.9}
                else:     
                    res = {"name": module, "type": module_type, "ipf": {}, "opf": {}, "ipv": {}, "opv": {}, 'mhh':0.9,'ohh': 0., "ihh": 0.0, "bp": np.array([0., 0.]), "width": 0.9}
                if module_type == "Operas" and self.modules[module].inputs:
                    ipfl = {}
                    deps = list(self.modules[module].inputs.keys())
                    nn = len(deps)
                    hh = 0.3 * nn - 0.2
                    for ii, dep in enumerate(deps):
                        ipfl[dep] = {
                            "nam": dep,
                            "pos": np.array([-0.35, ((nn - 1)/2 - ii) * 0.3]) + res["bp"],
                            "inc": [{"name": dep}],
                            "direct": True,
                        }
                    res["ipf"] = ipfl
                    if hh > 0.:
                        res["ihh"] = hh

                elif self.modules[module].input:
                    ipfl = {}
                    nn = len(self.modules[module].input)
                    hh = 0.7 * nn - 0.2
                    for ii in range(nn):
                        ipf = self.modules[module].input[ii]
                        if isinstance(ipf, dict):
                            ipf_name = ipf.get("name", f"input_{ii}")
                            ipf_path = ipf.get("path", f"{module}_input")
                            ipf_vars = ipf.get("variables", None)
                        else:
                            ipf_name = str(ipf)
                            ipf_path = f"{module}_input"
                            ipf_vars = None

                        ipfl[ipf_name] ={
                            "nam": ipf_name,
                            "pos": np.array([-1.1, ((nn - 1)/2 - ii) * 0.7 + 0.2]) + res["bp"],
                            "inc": [],
                            "fil": os.path.basename(ipf_path)
                        }
                        if isinstance(ipf_vars, dict):
                            for var in ipf_vars.values():
                                if "inc" in var.keys():
                                    for kk in var['inc']:
                                        if kk not in ipfl[ipf_name]['inc']:
                                            ipfl[ipf_name]['inc'].append(
                                                {
                                                    "name": str(kk)
                                                }
                                            )
                                else:
                                    ipfl[ipf_name]['inc'].append(var)
                        else:
                            deps = []
                            if isinstance(ipf, dict):
                                if ipf.get("_inc"):
                                    deps = [str(v) for v in ipf.get("_inc", [])]
                                elif ipf.get("entry"):
                                    deps = [str(ipf["entry"])]
                                elif ipf.get("name"):
                                    deps = [str(ipf["name"])]
                            else:
                                deps = [str(ipf)]
                            for dep in deps:
                                ipfl[ipf_name]['inc'].append(
                                            {
                                                "name": dep
                                            }
                                        )
                    res["ipf"] = ipfl
                    if hh > 0.:
                        res['ihh'] = hh 

                if self.modules[module].output and module_type != "Operas":
                    opfl = {}
                    nn = len(self.modules[module].output)
                    hh = 0.7 * nn - 0.2 
                    h0 = 0.
                    for ii in range(nn):
                        opf = self.modules[module].output[ii]
                        if isinstance(opf, dict):
                            opf_name = opf.get("name", f"output_{ii}")
                            opf_path = opf.get("path", f"{module}_output")
                            opf_vars = opf.get("variables", None)
                        else:
                            opf_name = str(opf)
                            opf_path = f"{module}_output"
                            opf_vars = None

                        if isinstance(opf_vars, list) and opf_vars:
                            inc_vars = [var['name'] for var in opf_vars if isinstance(var, dict) and 'name' in var]
                            if not inc_vars:
                                inc_vars = [opf_name]
                        else:
                            inc_vars = [opf_name]
                        nv = len(inc_vars)
                        if nv <= 2:
                            hl = 0.5 
                        else:
                            hl = 0.2 * nv
                        # print(opf['name'],)
                        op = {
                            "nam": opf_name,
                            "pos": np.array([1.1, -h0 - 0.5*hl -0.2 ]) + res["bp"],
                            "inc": inc_vars,
                            'fil': os.path.basename(opf_path)
                        }
                        h0 += hl 
                        h0 += 0.1
                        for jj, var_name in enumerate(inc_vars):
                            res['opv'][var_name] = {
                                "nam": var_name,
                                "pos": np.array([2.6, op["pos"][1] + ((nv-1)/2 - jj) * 0.2 ]) + res["bp"],
                                "wid": 0.,
                                "typ":  "Obser"
                            }
                        opfl[opf_name] = (op)
                    h0 -= 0.1
                    res["opf"] = opfl 
                    if h0 > 0.:
                        res['ohh'] = h0
                    for op in res['opf']:
                        res['opf'][op]['pos'] += np.array([0., + 0.5 * h0])
                    for op in res["opv"]:
                        res['opv'][op]["pos"] += np.array([0., + 0.5 * h0])

                elif self.modules[module].outputs:
                    nv = (len(self.modules[module].outputs) - 1)/2
                    h0 = nv * 0.2
                    ii = 0
                    for op, vv in self.modules[module].outputs.items():                  
                        res['opv'][op] = {
                            "nam": op, 
                            "pos": np.array([2.5, h0 - ii * 0.2]),
                            "wid": 0.,
                            "typ": "Obser",
                        }
                        if vv is not None: 
                            res['opv'][op]["typ"] = vv
                        ii += 1
                
                
                if res['ipf']:
                    res['width'] += 0.65 
                if res['opf']:
                    res['width'] += 0.65
                if res['opv']:
                    res['width'] += 1.5 
                # print(res)
                return res 

            def resolve_layer():
                for kk, item in self.calc_layer.items():
                    if kk == 1:
                        self.layer_info[1] = {
                            item["module"][0]: resolve_module(item["module"][0], kk)
                        }
                    else:
                        self.layer_info[kk] = {}
                        for module in item['module']:
                            self.layer_info[kk][module] = resolve_module(module, kk)

                layinfo = [None for ii in range(len(self.layer_info))]
                nlay = len(self.layer_info)
                for ii in range(nlay):
                    lid = len(self.layer_info) - ii
                    layer = self.layer_info[len(self.layer_info) - ii]
                    h0 = 0. 

                    y0 = {}
                    opvs = {}
                    for kk, item in layer.items():
                        mh = max([item['ihh'] + 0.4, item['ohh'] + 0.4, item['mhh']])
                        y00 = 0.5*mh + h0
                        item['bp'] += np.array([layerW * (lid - 1), y00])
                        h0 += mh 
                        y0[kk] = y00
                        for op, val in item['opv'].items():
                            opvs[op] = val
                        
                    layinfo[lid - 1] = {
                        "height": h0, 
                        "ipvars": {}, 
                        "y0":   y0, 
                        "opvars": opvs,
                        "bridge": {}
                    }
                for ii in range(nlay - 1):
                    clid = len(self.layer_info) - ii
                    # print(clid, "current layer")
                    layer = self.layer_info[clid]
                    layer_before = layinfo[clid - 2]
                    layer_current = layinfo[clid - 1]
                    # print(layer_before['opvars'])

                    bridge = {}
                    h0bf = layer_before['height']
                    for kk, item in layer.items():
                        # print(kk, item)
                        for iif, ipf in item['ipf'].items():
                            ipvs = {}
                            for ipvv in ipf['inc']:
                                # print(ipvv['name'], layer_before['opvars'], layer_before['bridge'])
                                # print("bg", ipvv, layer_before['opvars'].keys(), "end")
                                if ipvv["name"] not in layer_before['opvars'].keys() and ipvv['name'] not in layer_before['bridge'].keys():
                                    # print(ipvv, "Not Find in before in layer", clid - 1, len(layer_before["bridge"]))
                                    if not len(layer_before["bridge"]):
                                        h0bf += 1
                                    bridge[ipvv['name']] = {
                                        "nam": ipvv['name'],
                                        "pos": np.array([layerW*(clid - 2),  h0bf + 0.1])
                                    }
                                    layer_before['bridge'].update(bridge)
                                    ipvs[ipvv['name']] = {
                                        "nam": ipvv['name'],
                                        "pos": np.array([layerW*(clid - 2),  h0bf + 0.1])
                                    }
                                    h0bf += 0.2 
                                elif ipvv["name"] in layer_before['opvars'].keys(): 
                                    ipvs[ipvv['name']] = layer_before['opvars'][ipvv['name']]
                                elif ipvv['name'] in layer_before['bridge'].keys():
                                    # print(ipvv, "find in before", layer_before['bridge'])
                                    ipvs[ipvv['name']] = layer_before['bridge'][ipvv['name']]
                            
                            layer_current['ipvars'][iif] = [ipvs, ipf['pos']] 
                            layer_before['height'] = h0bf 
                
                for ii in range(nlay - 1):
                    clid = ii + 1 
                    layer = self.layer_info[clid]
                    layer_next = layinfo[clid]
                    if clid != 1:
                        for kk, mod in layer.items():
                            for off, opf in mod['opf'].items():
                                for op in opf['inc']:
                                    if op in layer_next['bridge'].keys():
                                        mod['opv'][op]['brg'] = layer_next['bridge'][op]['pos']
                    else: 
                        mod = layer['Parameters']
                        for op, opv in mod['opv'].items():
                            if op in layer_next['bridge'].keys():
                                opv['brg'] = layer_next['bridge'][op]['pos']
                
                maxh = max([layinfo[ii]['height'] for ii in range(nlay)]) + 0.4 
                figH = maxh
                for ii in range(nlay):
                    # print(layinfo[ii]['height'], layinfo[ii]['y0'], layinfo[ii]['bridge'])
                    lid =  ii + 1
                    layer = self.layer_info[lid]
                    for mm, mod in layer.items():
                        if ii != 0:
                            arrbp = mod['bp'] + np.array([0., 0.5 * (maxh - layinfo[ii]['height'])])
                        else: 
                            arrbp = np.array([0., 0.5 * (maxh - layinfo[ii]['height']) + mod['bp'][1]])
                        # update_module_pos(mod, arrbp)
                        mod['bp'] = mod['bp'] + np.array([0., 0.5 * (maxh - layinfo[ii]['height'])])
                        for oof, opf in mod['opf'].items():
                            opf['pos'] += arrbp
                        for oov, opv in mod['opv'].items():
                            opv['pos'] += arrbp 
                        for iif, ipf in mod['ipf'].items():
                            ipf['pos'] += arrbp
                for ii in range(nlay):
                    layer = self.layer_info[ii + 1]
                    layer_current = layinfo[ii]
                    for item in layer.values():
                        for kk, ipf in item['ipf'].items():
                            # pF = layer_current["ipvars"][kk][1]
                            # for op in layer_current["ipvars"][kk][0].values():
                            #     arraw_from_V2F(op['pos'], pF )
                            # print(kk, layer_current["ipvars"][kk])
                            ipf['link'] = layer_current["ipvars"][kk][0]
                            # break



                return layinfo, figH

            def update_module_pos(mod, arbp):
                from pprint import pprint
                pprint(mod)

            def draw_layer_module(klayer, mod):
                # print("Kleyaer", klayer)
                oh = np.array([0., + 0.5 * mod['ohh']])
                # ax.plot([(klayer-1)*layerW, (klayer-1)*layerW], [-100, 100], "-", c='grey', lw=0.6, alpha=0.2)
                if klayer != 1:
                    # print(mod)
                    module_at_pos(mod['bp'], mod.get("type", "Calculator"))
                    ax.text(mod['bp'][0], mod['bp'][1] - 0.4, mod['name'], ha="center", va='top', fontfamily="sans-serif", fontsize="medium", fontstyle="normal", fontweight="bold")
                    for kk, ipf in mod['ipf'].items():
                        if ipf.get("direct", False):
                            for link in ipf['link'].values():
                                arraw_from_V2F(link['pos'], mod['bp'])
                        else:
                            # print(kk, ipf.keys())
                            input_file_at_pos(ipf['pos'] )
                            ax.text(
                                ipf['pos'][0] - 0.25, ipf['pos'][1] - 0.27, ipf['fil'], 
                                ha="center", va='top', 
                                fontfamily="sans-serif", fontsize="x-small", fontstyle="normal", fontweight="light"
                            )
                            line_from_F2M(ipf['pos'] , mod['bp'])
                            for link in ipf['link'].values():
                                # print(kk, link)
                                arraw_from_V2F(link['pos'], ipf['pos'])
                    for op in mod['opv']:
                        mod['opv'][op]['pos'] = mod['opv'][op]['pos'] 
                        mod['opv'][op] = outvar_at_OP(mod['opv'][op])
                    # ax.plot([bs[0] + 2.6, bs[0] + 2.6], [-oh[1] - 0.2, oh[1] -0.2], "-", c='grey', lw=0.6, alpha=0.2)
                    for kk, opf in mod['opf'].items():
                        # print(kk, opf.keys())
                        output_file_at_pos(opf['pos'])
                        ax.text(
                            opf['pos'][0] + 0.25, opf['pos'][1] - 0.27, opf['fil'], 
                            ha="center", va='top', 
                            fontfamily="sans-serif", fontsize="x-small", fontstyle="normal", fontweight="light"
                        )
                        
                        line_from_M2F(opf['pos'], mod['bp'])

                        for op in opf['inc']:
                            opp = mod['opv'][op]
                            line_from_F2V(opf['pos'], opp)
                    if not mod['opf']:
                        for op in mod['opv'].values():
                            line_from_F2V(mod['bp'], op)
                else:
                    sampler_at_pos(mod['bp'])
                    ax.text(mod['bp'][0], mod['bp'][1] - 0.4, mod['name'], ha="center", va='top', fontfamily="sans-serif", fontsize="medium", fontstyle="normal", fontweight="bold")

                    # print(mod['opv'])
                    for op in mod['opv']:
                        mod['opv'][op] = outvar_at_OP(mod['opv'][op])
                        line_from_F2V(mod['bp'], mod['opv'][op])
                        if "brg" in mod['opv'][op].keys():
                            line_from_V2B(mod['opv'][op]['pos'], mod['opv'][op]['brg'])

            def arraw_from_V2F(pA, pB):
                tt = np.linspace(-0.5 * np.pi, 0.4 * np.pi, 100)
                xx = np.linspace(pA[0] + 0.14, pB[0] - 0.4, 100)
                yy = pA[1] + (np.sin(np.sin(tt)) + 0.8414709848078965)/1.6555005931675704 * (pB[1]-pA[1])
                points = np.array([xx, yy]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                cmap = plt.get_cmap('coolwarm')  
                norm = plt.Normalize(-1.5, 1.5)
                colors = cmap(norm(tt[:-1]))             
                lc = LineCollection(segments, colors=colors, linewidth=2)
                ax.add_collection(lc)
                ax.plot(pA[0] + 0.1, pA[1], 's', c="gray", markersize=4)
                ax.plot(pA[0] + 0.14, pA[1], 's', c="darkgray", markersize=6)

            def line_from_V2B(pV, pB):
                tt = np.linspace(-0.15 * np.pi, 0.5 * np.pi, 100)
                xx = np.linspace(pV[0] + 0.16, pB[0], 100)
                yy = pV[1] + (np.sin(np.sin(np.sin(np.sin(np.sin(tt))))) + 0.4004294116620898)/1.028001243711249 * (pB[1]-pV[1])
                ax.plot(xx, yy, "-", c="#3b4dc0", lw=2)

            def line_from_V2M(pA, pB):
                tt = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 100)
                xx = np.linspace(pA[0], pB[0], 100)
                yy = pA[1] + ((np.sin(np.sin(tt)) / 2 / 0.8414709848078965) + 0.5)  * (pB[1]-pA[1])
                ax.plot(xx, yy, "-", c="#3b4dc0", lw=2)

            def line_from_M2F(pF, pM):
                tt = np.linspace(0., 0.5*np.pi, 100)
                xx = np.linspace(0.3 + pM[0], 0.1 + pF[0], 100 )
                yy = pM[1] - 0.2 + (np.sin(np.sin(tt))) / 0.8414709848078965 * (pF[1] + 0.2 - pM[1])
                ax.plot(xx, yy, "-", c="#3b4dc0", lw=3)

            def line_from_F2V(pF, op):
                tt = np.linspace(-0.5*np.pi, 0.5*np.pi, 100)
                xx = np.linspace(pF[0] + 0.3, op['pos'][0] - op['wid'] - 0.1 , 100)
                yy = pF[1] + (np.sin(np.sin(tt)) / 0.8414709848078965  + 1)/2 * (op['pos'][1] - pF[1])  
                ax.plot(xx, yy, "-", c="#3b4dc0", lw=0.8, alpha=0.7)
                s = "o"
                if op['typ'] == "Param":
                    s = "s"
                elif op['typ'] == "Nuisa":
                    s = ">"
                ax.plot([op['pos'][0] - op['wid'] -0.1], [op['pos'][1]], s, markersize=4, c="#3b4dc0", alpha=1)

            def line_from_F2M(pF, pM):
                tt = np.linspace(0., 0.5*np.pi, 100)
                xx = np.linspace(-0.3 + pM[0], -0.1 + pF[0] , 100 )
                yy = pM[1] + 0.2 + (np.sin(np.sin(tt))) / 0.8414709848078965 * (pF[1] - 0.2 - pM[1])
                ax.plot(xx, yy, "-", c="#d45040", lw=3)

            def input_file_at_pos(pos):
                image_path = os.path.join(icon_dir, "inputfile.png")
                with Image.open(image_path) as image:
                    image = np.array(image)
                    ax.imshow(image, extent=[pos[0]-0.5, pos[0], pos[1]-0.25, pos[1]+0.25], zorder=100)

            def output_file_at_pos(pos):
                image_path = os.path.join(icon_dir, "outputfile.png")
                with Image.open(image_path) as image: 
                    image = np.array(image)
                    ax.imshow(image, extent=[pos[0], pos[0]+0.5, pos[1]-0.25, pos[1]+0.25], zorder=100)

            def outvar_at_OP(op):
                text = ax.text(op["pos"][0], op["pos"][1], op['nam'], ha="right", va="center", fontfamily="monospace", variant="small-caps", fontsize="x-small", fontstyle="normal", fontweight="bold")
                bbox = text.get_window_extent(renderer=fig.canvas.get_renderer())
                bbox_data = inv.transform([(bbox.x0, bbox.y0), (bbox.x1, bbox.y1)])
                op['wid'] = bbox_data[1][0] - bbox_data[0][0]
                return op 
                # ax.plot([op["pos"][0] - op['wid'], op["pos"][0]], [op["pos"][1], op["pos"][1]], "-")

            def module_at_pos(pos, module_type="Calculator"):
                icon_name = "calculator.png"
                if module_type == "Operas":
                    icon_name = "opera.png"
                image_path = os.path.join(icon_dir, icon_name)
                with Image.open(image_path) as image:
                    image = np.array(image.convert("RGBA"))  # Explicitly preserve transparency
                    ax.imshow(image, extent=[pos[0]-0.45, pos[0]+0.45, pos[1]-0.45, pos[1]+0.45], zorder=100)

            def lib_at_pos(pos):
                image_path = os.path.join(icon_dir, "library.png")
                with Image.open(image_path) as image:
                    image = np.array(image)
                    ax.imshow(image, extent=[pos[0]-0.45, pos[0]+0.45, pos[1]-0.45, pos[1]+0.45], zorder=100)

            def sampler_at_pos(pos):
                image_path = os.path.join(icon_dir, "sampler.png")
                with Image.open(image_path) as image:
                    image = np.array(image)
                    ax.imshow(image, extent=[pos[0]-0.45, pos[0]+0.45, pos[1]-0.45, pos[1]+0.45], zorder=100)

            def logo_at_pos(pos):
                image_path = os.path.join(icon_dir, "JarvisHEP.png")
                with Image.open(image_path) as image:
                    image = np.array(image.convert("RGBA"))
                    ax.imshow(image, extent=[pos[0]-0.25, pos[0]+0.25, pos[1]-0.25, pos[1]+0.25], zorder=100)
                    ax.text(pos[0]+0.27, pos[1]+0.15, "Jarvis-HEP", ha="left", va='top', color="#0F66C3", fontfamily="sans-serif", fontsize="small", fontstyle="normal", fontweight="bold")


            layerInfo, figH = resolve_layer()
            if figH < 1.8: 
                figH = 1.8
            # figH += 
            fig = plt.figure(figsize=(figL, figH))

            ax  = fig.add_axes([0., 0., 1., 1.])
            ax.axis("off")
            logo_at_pos([0.7, figH - 0.35])
            ax.set_xlim([0, figL])
            ax.set_ylim([0, figH])
            inv = ax.transData.inverted()

            # ax.plot([-100, 100], [0., 0.], "-", c='grey', lw=0.4, alpha=0.2)
            # ax.plot([-100, 100], [-0.2, -0.2], "-", c='grey', lw=0.4, alpha=0.2)
            # ax.plot([-100, 100], [0.2, 0.2], "-", c='grey', lw=0.4, alpha=0.2)    
            outframe = create_round_square([2, 4,])
            x_coords = [point[0] for point in outframe.exterior.coords]
            y_coords = [point[1] for point in outframe.exterior.coords]  

            for kk, item in self.layer_info.items():
                for mod, info in item.items():
                    draw_layer_module(kk, info)
            
            fig.savefig(save_path, dpi=300)
            plt.close(fig)
            self._release_matplotlib_font_handles()

    @staticmethod
    def _release_matplotlib_font_handles() -> None:
        """Best-effort release of matplotlib font file handles.

        Matplotlib may keep a persistent handle on selected TTF files
        (for example DejaVuSans-Bold.ttf) via its internal font cache.
        This method clears that cache when available.
        """
        try:
            from matplotlib import font_manager as _fm
        except Exception:
            return
        try:
            get_font = getattr(_fm, "_get_font", None)
            if hasattr(get_font, "cache_clear"):
                get_font.cache_clear()
        except Exception:
            # Best-effort cleanup; never fail workflow rendering on cache release.
            pass
