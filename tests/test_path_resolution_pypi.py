#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from unittest import mock


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.base import Base  # noqa: E402
from jarvishep.IOs.IOs import IOfile  # noqa: E402


class PathResolutionPyPITests(unittest.TestCase):
    def test_src_marker_points_to_package_root(self):
        base = Base()
        self.assertTrue(base.path["src_root"].endswith("jarvishep"))
        self.assertTrue(os.path.exists(base.path["args_info"]))
        self.assertTrue(os.path.exists(base.path["logger_config_path"]))
        self.assertTrue(os.path.exists(base.path["preference"]))

    def test_legacy_src_marker_is_not_remapped(self):
        base = Base()
        legacy = "&SRC/" + "src/card/preference.json"
        with self.assertRaises(ValueError):
            base.decode_path(legacy)

    def test_runtime_legacy_src_marker_is_not_remapped(self):
        base = Base()
        with tempfile.TemporaryDirectory() as tmpdir:
            base.configure_runtime_context(cwd=tmpdir)
            legacy = "&J/" + "src/card/preference.json"
            with self.assertRaises(ValueError):
                base.decode_path(legacy)

    def test_project_marker_infers_task_root_from_configs_dir(self):
        base = Base()
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "Demo")
            cfg_dir = os.path.join(project_root, "configs")
            os.makedirs(cfg_dir, exist_ok=True)
            marker = os.path.join(project_root, ".jarvis-project.json")
            with open(marker, "w", encoding="utf-8") as f1:
                f1.write("{}")
            fake_cfg = os.path.join(cfg_dir, "task.yaml")
            with open(fake_cfg, "w", encoding="utf-8") as f1:
                f1.write("Scan: {name: test, save_dir: '&J/outputs'}\n")

            inherited_root = os.path.join(tmpdir, "Inherited")
            os.makedirs(inherited_root, exist_ok=True)
            with mock.patch.dict(
                os.environ,
                {
                    "JARVIS_HEP_TASK_ROOT": inherited_root,
                    "JHEP_TASK_ROOT": inherited_root,
                },
                clear=False,
            ):
                base.configure_runtime_context(config_path=fake_cfg)
                self.assertEqual(
                    os.path.realpath(base.path["task_root"]),
                    os.path.realpath(project_root),
                )
                self.assertEqual(
                    os.path.realpath(os.getenv("JARVIS_HEP_TASK_ROOT", "")),
                    os.path.realpath(project_root),
                )

    def test_iofile_decode_path_prefers_project_marker_when_env_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "Standalone")
            nested = os.path.join(project_root, "configs", "sub")
            os.makedirs(nested, exist_ok=True)
            with open(os.path.join(project_root, ".jarvis-project.json"), "w", encoding="utf-8") as f1:
                f1.write("{}")

            old_cwd = os.getcwd()
            try:
                os.chdir(nested)
                with mock.patch.dict(
                    os.environ,
                    {"JARVIS_HEP_TASK_ROOT": "", "JHEP_TASK_ROOT": ""},
                    clear=False,
                ):
                    iofile = IOfile(
                        name="demo",
                        path="",
                        file_type="Json",
                        variables=[],
                        save=False,
                        logger=None,
                        PackID=None,
                        sample_save_dir="",
                        module="Demo",
                        funcs={},
                    )
                    self.assertEqual(
                        os.path.realpath(iofile.decode_path("&J/outputs")),
                        os.path.realpath(os.path.join(project_root, "outputs")),
                    )
            finally:
                os.chdir(old_cwd)


if __name__ == "__main__":
    unittest.main()
