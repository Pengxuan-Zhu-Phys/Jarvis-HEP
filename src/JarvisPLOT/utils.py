#!/usr/bin/env python3 
# JarvisPLOT utils: minimal helpers
from __future__ import annotations
from typing import Any, Mapping

__all__ = [
    "deep_merge",
    "rc_update_safely",
]

def deep_merge(a: Any, b: Any) -> Any:
    """Deep-merge two mapping-like objects (dicts). Right-hand side wins.
    - If either side is not a dict, return `b`.
    - For keys present in both and both values are dicts → merge recursively.
    - Otherwise, overwrite with the right-hand value.
    """
    if not isinstance(a, dict) or not isinstance(b, dict):
        return b
    out: dict[str, Any] = dict(a)
    for k, v in b.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def rc_update_safely(rc_params_like: Mapping[str, Any] | None) -> None:
    """Update matplotlib rcParams with a mapping, ignoring invalid keys.
    This keeps style bundles resilient across matplotlib versions.
    """
    if not rc_params_like:
        return
    try:
        import matplotlib as mpl
    except Exception:
        return
    for k, v in rc_params_like.items():
        try:
            mpl.rcParams[k] = v
        except Exception:
            # Silently ignore unknown/invalid rc keys to avoid breaking runs
            pass