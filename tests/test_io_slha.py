#!/usr/bin/env python3
"""D11.3: SLHA end-to-end via Portal V2 surface (no mocked adapters)."""

from __future__ import annotations

import asyncio
import os
import tempfile
import unittest

from jarvishep2.io_portal import (
    available_io_formats,
    build_io_context,
    get_io_registry,
    read_io_output,
    reset_io_registry_for_tests,
    write_io_input,
)

MINIMAL_SLHA = """\
Block MASS
   25  125.0  # higgs
Block MINPAR
   1   0.1
"""


class PortalV2SlhaTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_io_registry_for_tests()

    def tearDown(self) -> None:
        reset_io_registry_for_tests()

    def test_v2_registry_exposes_slha(self) -> None:
        formats = available_io_formats("input")
        self.assertIn("SLHA", formats)
        self.assertIn("JSON", formats)
        out_formats = available_io_formats("output")
        self.assertIn("SLHA", out_formats)
        self.assertIn("xSLHA", out_formats)

    def test_slha_write_and_read_end_to_end(self) -> None:
        registry = get_io_registry()
        adapter = registry.get("SLHA", "input")
        self.assertEqual(adapter.format_name, "SLHA")

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "input.slha")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(MINIMAL_SLHA)

            context = build_io_context(sample_info={"uuid": "slha-test"}, pack_id="001")

            async def _run():
                await write_io_input(
                    {
                        "name": "card",
                        "type": "SLHA",
                        "path": path,
                        "operations": [
                            {
                                "op": "set",
                                "name": "mH",
                                "block": "MASS",
                                "entry": 25,
                                "value": 130.0,
                            }
                        ],
                    },
                    {"mH": 130.0},
                    context=context,
                    path=path,
                )
                return await read_io_output(
                    {
                        "name": "card_out",
                        "type": "SLHA",
                        "path": path,
                        "variables": [
                            {"name": "mH", "block": "MASS", "entry": 25},
                        ],
                    },
                    context=context,
                    path=path,
                )

            observables = asyncio.run(_run())
            self.assertAlmostEqual(float(observables["mH"]), 130.0, places=6)


if __name__ == "__main__":
    unittest.main()
