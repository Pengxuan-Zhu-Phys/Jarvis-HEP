#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.project_scaffold import (  # noqa: E402
    LEGACY_COMPAT_SUBDIRS,
    PROJECT_DESCRIPTOR_NAME,
    PROJECT_MARKER_NAME,
    PROJECT_NESTED_SUBDIRS,
    PROJECT_SUBDIRS,
    create_project_scaffold,
)


class ProjectScaffoldTemplateTests(unittest.TestCase):
    def test_create_project_scaffold_writes_marker_and_template_assets(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = create_project_scaffold("StandaloneA", cwd=tmpdir)
            self.assertTrue(os.path.isdir(project_root))
            for subdir in PROJECT_SUBDIRS:
                self.assertTrue(os.path.isdir(os.path.join(project_root, subdir)))
            for subdir in PROJECT_NESTED_SUBDIRS:
                self.assertTrue(os.path.isdir(os.path.join(project_root, subdir)))
            for subdir in LEGACY_COMPAT_SUBDIRS:
                self.assertTrue(os.path.isdir(os.path.join(project_root, subdir)))
            self.assertTrue(os.path.exists(os.path.join(project_root, PROJECT_MARKER_NAME)))
            self.assertTrue(os.path.exists(os.path.join(project_root, PROJECT_DESCRIPTOR_NAME)))
            self.assertTrue(
                os.path.exists(
                    os.path.join(project_root, "bin", "quickstart_mcmc_operas.yaml")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(project_root, "data", "points.csv")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(project_root, "deps", "environment_default.yaml")
                )
            )
            self.assertFalse(os.path.exists(os.path.join(project_root, "outputs")))
            self.assertFalse(os.path.exists(os.path.join(project_root, "logs")))
            self.assertFalse(os.path.exists(os.path.join(project_root, "images")))
            self.assertFalse(os.path.exists(os.path.join(project_root, "configs")))
            self.assertFalse(os.path.exists(os.path.join(project_root, "scripts")))
            self.assertFalse(os.path.exists(os.path.join(project_root, "docs")))
            self.assertFalse(os.path.exists(os.path.join(project_root, "assets")))
            self.assertFalse(os.path.exists(os.path.join(project_root, "cache")))
            self.assertFalse(os.path.exists(os.path.join(project_root, "calculators")))

            with open(
                os.path.join(project_root, "bin", "quickstart_mcmc_operas.yaml"),
                "r",
                encoding="utf-8",
            ) as f1:
                mcmc_yaml = f1.read()
            with open(
                os.path.join(project_root, "bin", "quickstart_csv_operas.yaml"),
                "r",
                encoding="utf-8",
            ) as f1:
                csv_yaml = f1.read()

            expected_env_path = 'default_yaml_path: "&J/deps/environment_default.yaml"'
            self.assertIn(expected_env_path, mcmc_yaml)
            self.assertIn(expected_env_path, csv_yaml)

            with open(
                os.path.join(project_root, PROJECT_DESCRIPTOR_NAME),
                "r",
                encoding="utf-8",
            ) as f1:
                descriptor = f1.read()
            self.assertIn("artifact_roots:", descriptor)
            self.assertIn("logs: logs", descriptor)
            self.assertIn("images: images", descriptor)

            with open(
                os.path.join(project_root, "README.md"),
                "r",
                encoding="utf-8",
            ) as f1:
                scaffold_readme = f1.read()
            self.assertIn("Jarvis project create", scaffold_readme)
            self.assertNotIn("Jarvis --mkproject", scaffold_readme)


if __name__ == "__main__":
    unittest.main()
