#!/usr/bin/env python3
"""IO save:true is HEP FileOperation-owned (Portal only does format R/W)."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.async_subprocess import AsyncSubprocessScheduler
from jarvishep2.file_operation_service import FileOperationService
from jarvishep2.io_portal import reset_io_registry_for_tests
from jarvishep2.sample import Sample


class IOSaveTrueTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_io_registry_for_tests()

    def _assert_eggbox_save(
        self,
        *,
        file_ops: FileOperationService | None,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proj = Path(tmp)
            runtime = proj / "calc" / "001"
            runtime.mkdir(parents=True)
            (runtime / "input.json").write_text('{"xx": 0, "yy": 0}\n', encoding="utf-8")
            (runtime / "eggbox.py").write_text(
                "import json\n"
                "d = json.load(open('input.json'))\n"
                "json.dump({'z': float(d.get('xx', 0)) + float(d.get('yy', 0))}, "
                "open('output.json', 'w'))\n",
                encoding="utf-8",
            )
            sample_root = proj / "SAMPLE"
            sample_root.mkdir()
            cfg = {
                "name": "EggBox",
                "path": str(proj / "calc" / "@PackID"),
                "execution": {
                    "path": str(proj / "calc" / "@PackID"),
                    "commands": [
                        {"cmd": "python3 eggbox.py", "cwd": str(runtime)}
                    ],
                    "input": [
                        {
                            "name": "inpjson",
                            "type": "JSON",
                            "path": str(runtime / "input.json"),
                            "save": True,
                            "actions": [
                                {
                                    "type": "Dump",
                                    "variables": [
                                        {"name": "xx", "expression": "x"},
                                        {"name": "yy", "expression": "y"},
                                    ],
                                }
                            ],
                        }
                    ],
                    "output": [
                        {
                            "name": "oupjson",
                            "type": "JSON",
                            "path": str(runtime / "output.json"),
                            "save": True,
                            "variables": [{"name": "z"}],
                        }
                    ],
                },
            }
            mod = CalculatorModule("EggBox", cfg)
            mod.attach_scheduler(AsyncSubprocessScheduler())
            if file_ops is not None:
                mod.attach_file_ops(file_ops)
            mod.acquire_pack_id("001")
            sample = Sample(uuid="save-uuid")
            sample.info = {
                "uuid": "save-uuid",
                "params": {"x": 1.5, "y": 2.5},
                "observables": {"x": 1.5, "y": 2.5, "uuid": "save-uuid"},
                "_bucket_parent": str(sample_root),
                "sample_artifacts": "always",
                "workflow_has_calculator": True,
            }
            sample.materialize()
            result = mod.execute(sample.info)
            self.assertAlmostEqual(float(result["z"]), 4.0)
            self.assertIn("SAMPLE", str(result["inpjson"]))
            self.assertIn("SAMPLE", str(result["oupjson"]))
            self.assertTrue(os.path.isfile(result["inpjson"]))
            self.assertTrue(os.path.isfile(result["oupjson"]))
            self.assertTrue(result["inpjson"].endswith("input.json@EggBox"))
            self.assertTrue(result["oupjson"].endswith("output.json@EggBox"))

    def test_json_save_true_copies_into_sample_dir(self) -> None:
        """Fallback path: Calculator apply_hep_io_save without a service client."""
        self._assert_eggbox_save(file_ops=None)

    def test_json_save_true_via_file_operation_service(self) -> None:
        service = FileOperationService.start(mode="inline")
        try:
            self._assert_eggbox_save(file_ops=service)
        finally:
            service.shutdown()


if __name__ == "__main__":
    unittest.main()
