#!/usr/bin/env python3
from __future__ import annotations

import importlib


def test_bridson_module_imports_without_lib2to3_dependency():
    module = importlib.import_module("jarvishep.Sampling.bridson")
    assert hasattr(module, "Bridson")
