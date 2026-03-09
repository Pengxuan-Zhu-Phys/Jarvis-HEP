#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import unittest
import warnings


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.config import ConfigLoader  # noqa: E402
from jarvishep.distributor import Distributor  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class MultiNestExampleConfigTests(unittest.TestCase):
    def _validate_yaml(self, yaml_path: str):
        sampler = Distributor.set_method("MultiNest")
        loader = ConfigLoader()
        loader.logger = _NoopLogger()
        loader.load_config(yaml_path)
        loader.set_schema(sampler.schema)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", DeprecationWarning)
            loader.validate_config()
        remote_ref_warnings = [
            item
            for item in caught
            if isinstance(item.message, DeprecationWarning)
            and "Automatically retrieving remote references" in str(item.message)
        ]
        self.assertEqual(
            len(remote_ref_warnings),
            0,
            msg="jsonschema validation should not rely on implicit remote-reference retrieval.",
        )
        self.assertTrue(loader.config.get("Sampling", {}).get("Method") == "MultiNest")

    def test_multinest_calculator_example_schema_valid(self):
        self._validate_yaml(os.path.join(PROJECT_ROOT, "bin", "EggBox", "Example_MultiNest.yaml"))

    def test_multinest_calculator_quick_example_schema_valid(self):
        self._validate_yaml(os.path.join(PROJECT_ROOT, "bin", "EggBox", "Example_MultiNest_Quick.yaml"))

    def test_multinest_operas_example_schema_valid(self):
        self._validate_yaml(os.path.join(PROJECT_ROOT, "bin", "EggBox", "Example_MultiNest_Operas.yaml"))

    def test_multinest_operas_quick_example_schema_valid(self):
        self._validate_yaml(os.path.join(PROJECT_ROOT, "bin", "EggBox", "Example_MultiNest_Operas_Quick.yaml"))


if __name__ == "__main__":
    unittest.main()
