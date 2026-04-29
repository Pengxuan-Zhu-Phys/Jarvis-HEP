#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from types import SimpleNamespace


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.core import Core  # noqa: E402
from jarvishep.project_scaffold import create_project_scaffold  # noqa: E402


class ProjectRuntimeLayoutTests(unittest.TestCase):
    @staticmethod
    def _real(path: str) -> str:
        return os.path.realpath(path)

    def test_init_project_places_logs_and_flowchart_under_project_scope(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = create_project_scaffold("StandaloneA", cwd=tmpdir)
            yaml_path = os.path.join(
                project_root, "bin", "quickstart_mcmc_operas.yaml"
            )

            core = Core()
            core.configure_runtime_context(config_path=yaml_path)
            core.args = SimpleNamespace(file=yaml_path, plot=False)
            core.init_project()

            scan_name = "quickstart_mcmc_operas"
            outputs_root = os.path.join(project_root, "outputs", scan_name)
            logs_root = os.path.join(project_root, "logs", scan_name)
            images_root = os.path.join(project_root, "images", scan_name)

            self.assertEqual(
                self._real(core.info["sample"]["task_result_dir"]),
                self._real(outputs_root),
            )
            self.assertEqual(self._real(core.info["logs_dir"]), self._real(logs_root))
            self.assertEqual(self._real(core.info["images_dir"]), self._real(images_root))
            self.assertEqual(
                self._real(core.info["jarvis_log"]),
                self._real(os.path.join(logs_root, f"{scan_name}.log")),
            )
            self.assertEqual(
                self._real(core.info["sampler_log"]),
                self._real(os.path.join(logs_root, "MCMC.log")),
            )
            self.assertEqual(
                self._real(core.info["factory_log"]),
                self._real(os.path.join(logs_root, "Factory.log")),
            )
            self.assertEqual(
                self._real(core.info["flowchart_path"]),
                self._real(os.path.join(images_root, "flowchart.png")),
            )
            self.assertEqual(
                self._real(core.info["flowchart_semantic_path"]),
                self._real(os.path.join(images_root, "flowchart.json")),
            )
            self.assertFalse(
                self._real(core.info["jarvis_log"]).startswith(
                    self._real(outputs_root) + os.sep
                )
            )
            self.assertFalse(
                self._real(core.info["sampler_log"]).startswith(
                    self._real(outputs_root) + os.sep
                )
            )
            self.assertFalse(
                self._real(core.info["factory_log"]).startswith(
                    self._real(outputs_root) + os.sep
                )
            )
            self.assertFalse(
                self._real(core.info["flowchart_path"]).startswith(
                    self._real(outputs_root) + os.sep
                )
            )
            self.assertFalse(
                self._real(core.info["flowchart_semantic_path"]).startswith(
                    self._real(outputs_root) + os.sep
                )
            )

            self.assertTrue(os.path.isdir(outputs_root))
            self.assertTrue(os.path.isdir(os.path.join(outputs_root, "SAMPLE")))
            self.assertTrue(os.path.isdir(os.path.join(outputs_root, "DATABASE")))
            self.assertTrue(os.path.isdir(logs_root))
            self.assertTrue(os.path.isdir(images_root))

    def test_plot_mode_uses_scan_scoped_images_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = create_project_scaffold("StandaloneB", cwd=tmpdir)
            yaml_path = os.path.join(
                project_root, "bin", "quickstart_mcmc_operas.yaml"
            )

            core = Core()
            core.configure_runtime_context(config_path=yaml_path)
            core.args = SimpleNamespace(file=yaml_path, plot=True)
            core.init_project()

            scan_name = "quickstart_mcmc_operas"
            images_root = os.path.join(project_root, "images", scan_name)
            self.assertEqual(
                self._real(core.info["plot"]["save_path"]),
                self._real(images_root),
            )
            self.assertEqual(
                self._real(core.info["plot"]["config"]),
                self._real(os.path.join(images_root, f"{scan_name}.yaml")),
            )

    def test_check_modules_mode_uses_tests_subdirectory_for_sample_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = create_project_scaffold("StandaloneC", cwd=tmpdir)
            yaml_path = os.path.join(
                project_root, "bin", "quickstart_mcmc_operas.yaml"
            )

            core = Core()
            core.configure_runtime_context(config_path=yaml_path)
            core.mode = "1PC"
            core.args = SimpleNamespace(file=yaml_path, plot=False)
            core.init_project()

            scan_name = "quickstart_mcmc_operas"
            expected = os.path.join(project_root, "outputs", scan_name, "SAMPLE", "tests")
            self.assertEqual(self._real(core.info["sample"]["sample_dirs"]), self._real(expected))
            self.assertTrue(os.path.isdir(expected))


if __name__ == "__main__":
    unittest.main()
