# layers/scatter.py (minimal, using unified lim/scale/label if objects are passed)
import numpy as np
from .base import LayerBase
from matplotlib.colors import Normalize, LogNorm

class ScatterLayer(LayerBase):
    TYPE = "scatter"
    def render(self, ax, ctx):
        cfg = self.spec.cfg
        # x/y/c may be str or dict variable objects
        x = _resolve_field(cfg["x"], ctx)
        y = _resolve_field(cfg["y"], ctx)
        c = cfg.get("c")
        mappable = None

        # apply axis meta if dict variables
        _apply_axis(ax, 'x', cfg["x"], x)
        _apply_axis(ax, 'y', cfg["y"], y)

        style = {k:v for k,v in cfg.items() if k not in {"x","y","c","axes","type","name"}}
        if c is not None:
            c_arr, c_meta = _resolve_color(c, ctx)
            norm = _build_norm(c_meta, c_arr)
            sc = ax.scatter(x[0], y[0], c=c_arr, norm=norm, cmap=c_meta.get("cmap") or style.pop("cmap", None), **style)
            mappable = sc
        else:
            sc = ax.scatter(x[0], y[0], **style)
        return mappable

def _resolve_field(val, ctx):
    # returns (arr, meta) — here meta for axis; keep simple for now
    ns = ctx.namespace
    if isinstance(val, str): return (eval(val, ns, {}), {})
    if isinstance(val, dict):
        expr = f"{val['dataset']}.{val['expr']}" if val.get('dataset') else val['expr']
        arr = eval(expr, ns, {})
        return (arr, val)  # carry lim/scale/label inside
    return (np.asarray(val), {})

def _apply_axis(ax, which, var_obj, resolved):
    arr, meta = resolved
    if isinstance(var_obj, dict):
        if lab := meta.get("label"): getattr(ax, f"set_{which}label")(lab)
        if sc := meta.get("scale"): getattr(ax, f"set_{which}scale")(sc)
        lim = meta.get("lim")
        if isinstance(lim, list) and len(lim)==2:
            getattr(ax, f"set_{which}lim")(lim)

def _resolve_color(val, ctx):
    ns = ctx.namespace
    if isinstance(val, str): return eval(val, ns, {}), {}
    expr = f"{val['dataset']}.{val['expr']}" if val.get('dataset') else val['expr']
    arr = eval(expr, ns, {})
    return arr, val

def _build_norm(meta, data):
    nm = (meta or {}).get("norm", "linear")
    lim = (meta or {}).get("lim")
    vmin = min(data) if lim=="auto" or lim is None else lim[0]
    vmax = max(data) if lim=="auto" or lim is None else lim[1]
    if nm == "log": return LogNorm(vmin=max(vmin,1e-12), vmax=vmax)
    return Normalize(vmin=vmin, vmax=vmax)