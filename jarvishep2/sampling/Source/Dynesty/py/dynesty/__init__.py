#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
dynesty is nested sampling package.
The main functionality of dynesty is performed by the
dynesty.NestedSampler and dynesty.DynamicNestedSampler
classes
"""
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("dynesty")
except PackageNotFoundError:
    # Vendored into jarvishep2 — not installed as top-level 'dynesty'.
    __version__ = "3.1.0+jarvishep2"


from .dynesty import NestedSampler, DynamicNestedSampler
from . import bounding
from . import utils
from . import pool

# Route dynesty warnings + default progress to Jarvis logger when available.
try:
    from .jarvis_logging import install_warnings_bridge

    install_warnings_bridge()
except Exception:
    pass
