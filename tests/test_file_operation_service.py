#!/usr/bin/env python3
"""FileOperationService: HEP-owned SAMPLE save/copy/delete process."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from jarvishep2.file_operation_service import FileOperationService
from jarvishep2.file_ops import apply_io_save_policy, save_io_copy
from jarvishep2.proc_title import file_operation_title


class SaveIoCopyTests(unittest.TestCase):
    def test_save_true_copies_into_sample_dir(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "runtime" / "input.json"
            src.parent.mkdir(parents=True)
            src.write_text('{"x": 1}\n', encoding="utf-8")
            sample = Path(tmp) / "SAMPLE" / "uuid-1"
            dest = save_io_copy(
                str(src),
                str(sample),
                module="EggBox",
                save=True,
            )
            self.assertIsNotNone(dest)
            assert dest is not None
            self.assertTrue(dest.endswith("input.json@EggBox"))
            self.assertTrue(os.path.isfile(dest))
            self.assertEqual(Path(dest).read_text(encoding="utf-8"), '{"x": 1}\n')

    def test_output_default_temp_when_save_false(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "runtime" / "output.json"
            src.parent.mkdir(parents=True)
            src.write_text('{"z": 2}\n', encoding="utf-8")
            sample = Path(tmp) / "SAMPLE" / "uuid-1"
            dest = apply_io_save_policy(
                source_path=str(src),
                sample_save_dir=str(sample),
                module="EggBox",
                spec={"save": False},
                direction="output",
            )
            self.assertIsNotNone(dest)
            assert dest is not None
            self.assertIn(os.path.join(".temp", "output.json@EggBox"), dest.replace("/", os.sep))
            self.assertTrue(os.path.isfile(dest))

    def test_input_no_copy_when_save_false(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "runtime" / "input.json"
            src.parent.mkdir(parents=True)
            src.write_text("{}", encoding="utf-8")
            sample = Path(tmp) / "SAMPLE" / "uuid-1"
            dest = apply_io_save_policy(
                source_path=str(src),
                sample_save_dir=str(sample),
                module="EggBox",
                spec={"save": False},
                direction="input",
            )
            self.assertIsNone(dest)


class FileOperationServiceTests(unittest.TestCase):
    def test_file_operation_title_includes_scan_name(self) -> None:
        self.assertEqual(file_operation_title(scan_name="Demo"), "Jarvis2-FileOperation:Demo")

    def test_inline_save_and_delete(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "out.json"
            src.write_text('{"z": 1}\n', encoding="utf-8")
            sample = Path(tmp) / "SAMPLE"
            service = FileOperationService.start(mode="inline", delete_method="shutil")
            try:
                dest = service.save_io(
                    source_path=str(src),
                    sample_save_dir=str(sample),
                    module="Calc",
                    spec={"save": True, "name": "oup"},
                    direction="output",
                )
                self.assertIsNotNone(dest)
                assert dest is not None
                self.assertTrue(os.path.isfile(dest))
                service.delete([dest], missing_ok=False)
                self.assertFalse(os.path.lexists(dest))
            finally:
                service.shutdown()

    def test_process_mode_save_io(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "out.json"
            src.write_text('{"z": 3}\n', encoding="utf-8")
            sample = Path(tmp) / "SAMPLE"
            service = FileOperationService.start(mode="process", delete_method="shutil")
            try:
                dest = service.save_io(
                    source_path=str(src),
                    sample_save_dir=str(sample),
                    module="Calc",
                    spec={"save": True},
                    direction="output",
                )
                self.assertIsNotNone(dest)
                assert dest is not None
                self.assertTrue(os.path.isfile(dest))
                self.assertTrue(dest.endswith("out.json@Calc"))
            finally:
                service.shutdown()


if __name__ == "__main__":
    unittest.main()
