# JarvisPLOT variables: resolve dataset+expr and name+expr
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class VariableSpec:
    dataset: Optional[str]
    name: Optional[str]
    expr: str
    lim: object = None
    scale: str | None = None
    label: str | None = None
    cmap: str | None = None
    norm: str | None = None

@dataclass
class VariableSet:
    items: List[VariableSpec]

    def resolve_all(self, data_env):
        ns = data_env.namespace()
        # 1) dataset+expr without name → evaluate and expose by label/expr (optional)
        for v in self.items:
            if v.dataset and not v.name:
                try:
                    val = eval(f"{v.dataset}.{v.expr}", ns, {})
                    key = v.label or v.expr
                    ns[key] = val
                except Exception:
                    pass
        # 2) named variables (computed or named dataset expr)
        for v in self.items:
            if v.name:
                try:
                    expr = f"{v.dataset}.{v.expr}" if v.dataset else v.expr
                    val = eval(expr, ns, {})
                    ns[v.name] = val
                except Exception:
                    pass
        return ns