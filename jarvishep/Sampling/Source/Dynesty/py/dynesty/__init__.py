#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
dynesty is nested sampling package.
The main functionality of dynesty is performed by the
dynesty.NestedSampler and dynesty.DynamicNestedSampler
classes
"""
from importlib import import_module

from ._version import __version__

__all__ = [
    "NestedSampler",
    "DynamicNestedSampler",
    "bounding",
    "utils",
    "pool",
    "__version__",
]


def __getattr__(name):
    if name in {"NestedSampler", "DynamicNestedSampler"}:
        from .dynesty import DynamicNestedSampler, NestedSampler

        return {
            "NestedSampler": NestedSampler,
            "DynamicNestedSampler": DynamicNestedSampler,
        }[name]
    if name == "bounding":
        return import_module(f"{__name__}.bounding")
    if name == "utils":
        return import_module(f"{__name__}.utils")
    if name == "pool":
        return import_module(f"{__name__}.pool")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
