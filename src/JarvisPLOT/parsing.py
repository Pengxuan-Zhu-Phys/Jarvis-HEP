#!/usr/bin/env python3 
# parsing.py
from dataclasses import dataclass
from typing import List, Dict, Any
from .figures.figure import FigureSpec, AxesSpec, LayerSpec


@dataclass
class OutputSpec:
    dir: str
    formats: List[str]
    dpi: int

@dataclass
class DataSourceSpec:
    name: str
    type: str
    path: str
    options: Dict[str, Any] = None

@dataclass
class VariableSpec:
    # supports either dataset+expr OR name+expr
    dataset: str | None
    name: str | None
    expr: str
    lim: str | list | None
    scale: str | None
    label: str | None
    cmap: str | None
    norm: str | None

@dataclass
class LayerSpec:
    name: str | None
    type: str
    axes: str
    cfg: Dict[str, Any]   # layer-specific fields, includes x/y/c (objects or strings)

@dataclass
class AxesSpec:
    name: str
    preset: str | None    # optional explicit preset
    cfg: Dict[str, Any]   # rect/projection/grid... (overrides)

@dataclass
class FigureSpec:
    name: str
    size: list | None
    axes: List[AxesSpec]
    layers: List[LayerSpec]
    datasources: List[DataSourceSpec] | None = None

@dataclass
class ConfigSpec:
    yaml_dir: str
    style_name: str | None
    output: OutputSpec
    datasources: List[DataSourceSpec]
    variables: 'VariableSet'     # wrapper with resolve_all()
    figures: List[FigureSpec]


# --- YAML config loader and helpers ---
import os, yaml
from .variables import VariableSet

def _as_output(o: Dict[str, Any]) -> OutputSpec:
    if not o:
        o = {}
    return OutputSpec(
        dir=o.get("dir", "."),
        formats=o.get("formats", ["png"]),
        dpi=int(o.get("dpi", 150)),
    )

def _as_datasources(lst: List[Dict[str, Any]] | None) -> List[DataSourceSpec]:
    out: List[DataSourceSpec] = []
    for d in (lst or []):
        out.append(DataSourceSpec(
            name=d["name"],
            type=d.get("type", "csv"),
            path=d["path"],
            options={k:v for k,v in d.items() if k not in {"name","type","path"}},
        ))
    return out

def _as_variables(lst: List[Dict[str, Any]] | None) -> VariableSet:
    from dataclasses import asdict
    items: List[VariableSpec] = []
    for v in (lst or []):
        items.append(VariableSpec(
            dataset=v.get("dataset"),
            name=v.get("name"),
            expr=v.get("expr"),
            lim=v.get("lim"),
            scale=v.get("scale"),
            label=v.get("label"),
            cmap=v.get("cmap"),
            norm=v.get("norm"),
        ))
    return VariableSet(items)

def _as_axes(lst: List[Dict[str, Any]] | None) -> List[AxesSpec]:
    out: List[AxesSpec] = []
    for a in (lst or []):
        name = a["name"]
        preset = a.get("preset")
        cfg = {k:v for k,v in a.items() if k not in {"name","preset"}}
        out.append(AxesSpec(name=name, preset=preset, cfg=cfg))
    return out

def _as_layers(lst: List[Dict[str, Any]] | None) -> List[LayerSpec]:
    out: List[LayerSpec] = []
    for l in (lst or []):
        ltype = l["type"]
        axes = l.get("axes","main")
        name = l.get("name")
        cfg = {k:v for k,v in l.items() if k not in {"type","axes","name"}}
        out.append(LayerSpec(name=name, type=ltype, axes=axes, cfg=cfg))
    return out

def _as_figures(lst: List[Dict[str, Any]] | None) -> List[FigureSpec]:
    out: List[FigureSpec] = []
    for f in (lst or []):
        name = f["name"]
        size = f.get("size")
        axes = _as_axes(f.get("axes"))
        layers = _as_layers(f.get("layers"))
        out.append(FigureSpec(name=name, size=size, axes=axes, layers=layers))
    return out

def load_config_spec(yaml_path: str) -> ConfigSpec:
    yaml_abs = os.path.abspath(yaml_path)
    yaml_dir = os.path.dirname(yaml_abs)
    with open(yaml_abs,"r",encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    style_name = data.get("style") or None
    output = _as_output(data.get("output", {}))
    datasources = _as_datasources(data.get("datasources"))
    variables = _as_variables(data.get("variables"))
    figures = _as_figures(data.get("figures"))
    return ConfigSpec(
        yaml_dir=yaml_dir,
        style_name=style_name or "paper_1x1",
        output=output,
        datasources=datasources,
        variables=variables,
        figures=figures,
    )