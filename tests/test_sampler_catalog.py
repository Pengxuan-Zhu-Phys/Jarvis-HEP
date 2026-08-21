#!/usr/bin/env python3
"""D25.1/D25.2: builtin sampler catalog is the method-name source of truth."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

from jarvishep2.contracts.methods import (
    _MCMC_METHODS,
    _METHODS_NEED_VARIABLES,
    _NESTED_METHODS,
)
from jarvishep2.distributor import Distributor
from jarvishep2.mapper import STATISTICAL_METHODS
from jarvishep2.sampler_catalog import (
    available_methods,
    builtin_specs,
    mcmc_names,
    methods_need_variables,
    names,
    nested_names,
    schema_ids,
    statistical_names,
)
from jarvishep2.task_schema import MANIFEST_PATH


class SamplerCatalogConsistencyTests(unittest.TestCase):
    def test_catalog_matches_manifest_and_distributor(self) -> None:
        with MANIFEST_PATH.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
        method_schemas = manifest["sampling_methods"]
        self.assertEqual(set(method_schemas), set(names()))
        self.assertEqual(dict(method_schemas), schema_ids())
        self.assertEqual(set(method_schemas), set(Distributor.available_methods()))
        self.assertEqual(available_methods(), sorted(names()))

    def test_every_spec_has_a_factory(self) -> None:
        registered = set(Distributor.available_methods())
        self.assertEqual({spec.name for spec in builtin_specs()}, registered)

    def test_family_views_match_contracts_and_mapper(self) -> None:
        self.assertEqual(_NESTED_METHODS, nested_names())
        self.assertEqual(_MCMC_METHODS, mcmc_names())
        self.assertEqual(_METHODS_NEED_VARIABLES, methods_need_variables())
        self.assertEqual(STATISTICAL_METHODS, statistical_names())
        self.assertEqual(statistical_names(), nested_names() | mcmc_names())
        self.assertNotIn("CSV", methods_need_variables())
        self.assertIn("CSV", names())


class SamplerCatalogImportIsolationTests(unittest.TestCase):
    _BANNED = (
        "jarvishep2.worker",
        "jarvishep2.Sampling.mcmc_sampler",
        "jarvishep2.Sampling.adaptive_bridson",
        "jarvishep2.distributor",
        "jarvishep2.core",
    )

    def _assert_clean_import(self, module: str) -> None:
        banned = ", ".join(repr(name) for name in self._BANNED)
        script = f"""
import sys
import {module}
loaded = [name for name in ({banned},) if name in sys.modules]
if loaded:
    raise SystemExit("unexpected modules: " + ",".join(loaded))
"""
        completed = subprocess.run(
            [sys.executable, "-c", script],
            cwd=str(Path(__file__).resolve().parents[1]),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_catalog_import_is_light(self) -> None:
        self._assert_clean_import("jarvishep2.sampler_catalog")

    def test_task_validation_import_is_light(self) -> None:
        self._assert_clean_import("jarvishep2.task_validation")

    def test_client_import_is_light(self) -> None:
        self._assert_clean_import("jarvishep2.client")

    def test_package_import_is_light(self) -> None:
        self._assert_clean_import("jarvishep2")

    def test_distributor_register_does_not_import_samplers(self) -> None:
        """Factories stay lazy; registering builtins must not import MCMC/ALS."""
        script = """
import sys
import jarvishep2.distributor
loaded = [
    name for name in (
        'jarvishep2.Sampling.mcmc_sampler',
        'jarvishep2.Sampling.adaptive_bridson',
        'jarvishep2.Sampling.dynesty_sampler',
    )
    if name in sys.modules
]
if loaded:
    raise SystemExit('unexpected modules: ' + ','.join(loaded))
"""
        completed = subprocess.run(
            [sys.executable, "-c", script],
            cwd=str(Path(__file__).resolve().parents[1]),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
