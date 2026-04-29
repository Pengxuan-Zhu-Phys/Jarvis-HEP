#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tarfile
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.project_packager import (  # noqa: E402
    ProjectNotFoundError,
    create_project_package,
)


class ProjectPackagerTests(unittest.TestCase):
    def _make_demo_project(self, root: str) -> None:
        os.makedirs(root, exist_ok=True)
        for dname in (
            "bin",
            "data",
            "outputs",
            "images",
            "logs",
            "checkpoints",
            "calculators",
            "deps",
        ):
            os.makedirs(os.path.join(root, dname), exist_ok=True)

        with open(os.path.join(root, ".jarvis-project.json"), "w", encoding="utf-8") as f1:
            f1.write("{}\n")
        with open(os.path.join(root, "jarvis.project.yaml"), "w", encoding="utf-8") as f1:
            f1.write("project: {}\n")
        with open(os.path.join(root, "README.md"), "w", encoding="utf-8") as f1:
            f1.write("# Demo\n")
        with open(os.path.join(root, ".hidden-note"), "w", encoding="utf-8") as f1:
            f1.write("secret\n")

        with open(os.path.join(root, "bin", "task.yaml"), "w", encoding="utf-8") as f1:
            f1.write("Scan: {name: demo, save_dir: '&J/outputs'}\n")
        with open(os.path.join(root, "data", "points.csv"), "w", encoding="utf-8") as f1:
            f1.write("uuid,x,y\n0,1,2\n")
        with open(os.path.join(root, "outputs", "result.txt"), "w", encoding="utf-8") as f1:
            f1.write("done\n")
        with open(os.path.join(root, "images", "plot.txt"), "w", encoding="utf-8") as f1:
            f1.write("plot\n")
        with open(os.path.join(root, "logs", "run.log"), "w", encoding="utf-8") as f1:
            f1.write("log\n")
        with open(os.path.join(root, "checkpoints", "state.ckpt"), "w", encoding="utf-8") as f1:
            f1.write("checkpoint\n")
        with open(os.path.join(root, "calculators", "runtime.txt"), "w", encoding="utf-8") as f1:
            f1.write("calc\n")
        with open(os.path.join(root, "deps", "pkg.txt"), "w", encoding="utf-8") as f1:
            f1.write("dep\n")

    def _tar_names(self, archive_path: str) -> set[str]:
        with tarfile.open(archive_path, "r:gz") as tf:
            return set(tf.getnames())

    def test_share_repro_full_profile_selection(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)

            share_report = create_project_package(project_root, profile="share")
            repro_report = create_project_package(project_root, profile="repro")
            full_report = create_project_package(project_root, profile="full")

            self.assertEqual(os.path.dirname(share_report.archive_path), tmpdir)
            self.assertEqual(os.path.dirname(repro_report.archive_path), tmpdir)
            self.assertEqual(os.path.dirname(full_report.archive_path), tmpdir)

            share_names = self._tar_names(share_report.archive_path)
            repro_names = self._tar_names(repro_report.archive_path)
            full_names = self._tar_names(full_report.archive_path)

            self.assertIn("DemoProject/bin/task.yaml", share_names)
            self.assertIn("DemoProject/data/points.csv", share_names)
            self.assertNotIn("DemoProject/deps/pkg.txt", share_names)
            self.assertNotIn("DemoProject/calculators/runtime.txt", share_names)
            self.assertNotIn("DemoProject/logs/run.log", share_names)
            self.assertNotIn("DemoProject/checkpoints/state.ckpt", share_names)

            self.assertIn("DemoProject/deps/pkg.txt", repro_names)
            self.assertIn("DemoProject/calculators/runtime.txt", repro_names)
            self.assertIn("DemoProject/logs/run.log", repro_names)
            self.assertIn("DemoProject/checkpoints/state.ckpt", repro_names)
            self.assertNotIn("DemoProject/.hidden-note", repro_names)

            self.assertIn("DemoProject/.hidden-note", full_names)
            self.assertIn("DemoProject/.jarvis-pack/manifest.json", full_names)
            self.assertIn("DemoProject/.jarvis-pack/checksums.sha256", full_names)

    def test_default_target_uses_cwd(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)
            old_cwd = os.getcwd()
            try:
                os.chdir(project_root)
                report = create_project_package(None, profile="repro")
            finally:
                os.chdir(old_cwd)
            self.assertTrue(os.path.exists(report.archive_path))
            self.assertEqual(os.path.realpath(report.project_root), os.path.realpath(project_root))

    def test_non_project_path_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ProjectNotFoundError):
                create_project_package(tmpdir, profile="repro")


if __name__ == "__main__":
    unittest.main()
