#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.config import ConfigValidator  # noqa: E402


class _NoopLogger:
    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None

    def info(self, *_args, **_kwargs):
        return None


class ConfigValidatorSchemaBlockTests(unittest.TestCase):
    def test_validate_yaml_skips_missing_schema_blocks_without_keyerror(self):
        validator = ConfigValidator()
        validator.logger = _NoopLogger()
        validator.set_config({})
        validator.schema = {
            "type": "object",
            "schemaBlock": {
                "input": {"$ref": "PATH-TO-input-JSON-SCHEMA"},
                "output": {"$ref": "PATH-TO-output-JSON-SCHEMA"},
            },
        }

        validator.validate_yaml()

        self.assertTrue(validator.passcheck)
        self.assertIn("input", validator.schema["schemaBlock"])
        self.assertIn("output", validator.schema["schemaBlock"])
        self.assertNotIn("Nuisance", validator.schema["schemaBlock"])


if __name__ == "__main__":
    unittest.main()
