#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import textwrap
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
    def _write_example_yaml(self, *, quick: bool, use_operas: bool) -> str:
        maxcall = 25 if quick else 250
        operas_block = ""
        if use_operas:
            operas_block = """
Operas:
  make_paraller: 1
  Modules:
    - name: EggBox
      operator: "helper.eggbox2d"
      call_mode: call
      required_modules: []
      input:
        - {name: x, expression: "x"}
        - {name: y, expression: "y"}
      output:
        - {name: z, entry: z}
"""

        yaml_text = textwrap.dedent(
            f"""
            Scan:
              name: "test_multinest"
              save_dir: "&J/outputs"
            Sampling:
              Method: "MultiNest"
              Variables:
                - name: x
                  description: "x"
                  distribution:
                    type: Flat
                    parameters:
                      min: 0
                      max: 1
                - name: y
                  description: "y"
                  distribution:
                    type: Flat
                    parameters:
                      min: 0
                      max: 1
              Bounds:
                nlive: 50
                rseed: 1
                run_nested:
                  maxcall: {maxcall}
                  print_progress: false
              LogLikelihood:
                - name: LogL_Z
                  expression: "z"
            EnvReqs:
              OS:
                - name: Darwin
                  version: ">=10.14"
              CERN_ROOT:
                required: false
                version: ">=0"
                get_path_command: "echo"
                Dependencies: []
              Python:
                version: ">=3.10"
                Dependencies: []
            Calculators:
              make_paraller: 1
              Modules:
                - name: EggBox
                  required_modules: []
                  clone_shadow: false
                  installation: []
                  initialization: []
                  execution:
                    path: "&J/calculators/runtime/program"
                    commands:
                      - "echo run"
                    input:
                      - name: params
                        path: "input.json"
                        type: JSON
                        save: true
                        actions:
                          - type: Dump
                            variables:
                              - name: x
                              - name: y
                    output:
                      - name: observables
                        path: "output.json"
                        type: JSON
                        save: false
                        variables:
                          - name: z
            """
        ).strip()
        if operas_block:
            yaml_text = f"{yaml_text}\n{operas_block.strip()}\n"
        else:
            yaml_text = f"{yaml_text}\n"

        tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(tempdir.cleanup)
        yaml_path = os.path.join(tempdir.name, "example_multinest.yaml")
        with open(yaml_path, "w", encoding="utf-8") as f1:
            f1.write(yaml_text)
        return yaml_path

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
        self._validate_yaml(self._write_example_yaml(quick=False, use_operas=False))

    def test_multinest_calculator_quick_example_schema_valid(self):
        self._validate_yaml(self._write_example_yaml(quick=True, use_operas=False))

    def test_multinest_operas_example_schema_valid(self):
        self._validate_yaml(self._write_example_yaml(quick=False, use_operas=True))

    def test_multinest_operas_quick_example_schema_valid(self):
        self._validate_yaml(self._write_example_yaml(quick=True, use_operas=True))


if __name__ == "__main__":
    unittest.main()
