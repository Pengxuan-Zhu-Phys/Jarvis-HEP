#!/usr/bin/env python3
"""D25.10: layer directories keep legacy import paths working."""

from __future__ import annotations

import unittest
from importlib import import_module, metadata


class PackageLayoutCompatibilityTests(unittest.TestCase):
    def test_sampling_legacy_and_canonical_bridson_are_the_same(self) -> None:
        from jarvishep2.Sampling.bridson import Bridson as legacy
        from jarvishep2.sampling.bridson import Bridson as canonical

        self.assertIs(legacy, canonical)
        self.assertEqual(canonical.__module__, "jarvishep2.sampling.bridson")

    def test_module_legacy_and_canonical_calculator_are_the_same(self) -> None:
        from jarvishep2.Module.calculator import CalculatorModule as legacy
        from jarvishep2.calculators.calculator import CalculatorModule as canonical

        self.assertIs(legacy, canonical)

    def test_contracts_legacy_and_canonical_are_the_same(self) -> None:
        from jarvishep2.contracts.methods import validate_method_sampling as legacy
        from jarvishep2.cardload.contracts.methods import (
            validate_method_sampling as canonical,
        )

        self.assertIs(legacy, canonical)

    def test_sampling_virtial_alias_is_sampling_virtual(self) -> None:
        from jarvishep2.Sampling.sampler import SamplingVirtial
        from jarvishep2.sampling.sampler import SamplingVirtual

        self.assertIs(SamplingVirtial, SamplingVirtual)
        self.assertEqual(SamplingVirtual.__name__, "SamplingVirtual")

    def test_client_entry_point_still_resolves_main(self) -> None:
        from jarvishep2.cli.client import main as canonical
        from jarvishep2.client import main as legacy

        self.assertIs(legacy, canonical)
        entry = metadata.entry_points().select(group="console_scripts", name="Jarvis")
        matches = list(entry)
        self.assertTrue(matches, "console script Jarvis is not installed")
        self.assertEqual(matches[0].value, "jarvishep2.client:main")

    def test_core_factory_worker_shims_are_identity(self) -> None:
        self.assertIs(import_module("jarvishep2.core"), import_module("jarvishep2.runtime.core"))
        self.assertIs(
            import_module("jarvishep2.factory"),
            import_module("jarvishep2.runtime.factory"),
        )
        self.assertIs(
            import_module("jarvishep2.worker"),
            import_module("jarvishep2.runtime.worker"),
        )
        self.assertIs(
            import_module("jarvishep2.redis_queue"),
            import_module("jarvishep2.queue.redis_queue"),
        )


if __name__ == "__main__":
    unittest.main()
