#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import io
import os
import sys
import tarfile
import tempfile
import unittest

import yaml


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.project_packager import (  # noqa: E402
    ProjectPackError,
    ProjectNotFoundError,
    create_project_pack_manifest,
    create_project_package,
    create_project_package_from_manifest,
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

    def test_manifest_generation_creates_incrementing_yaml_without_archive(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            manifest_dir = os.path.join(tmpdir, "manifests")
            self._make_demo_project(project_root)

            first = create_project_pack_manifest(
                project_root,
                profile="repro",
                manifest_dir=manifest_dir,
            )
            second = create_project_pack_manifest(
                project_root,
                profile="repro",
                manifest_dir=manifest_dir,
            )

            self.assertRegex(
                os.path.basename(first.manifest_path),
                r"^pack_\d{8}_001\.yaml$",
            )
            self.assertRegex(
                os.path.basename(second.manifest_path),
                r"^pack_\d{8}_002\.yaml$",
            )
            self.assertTrue(os.path.exists(first.manifest_path))
            self.assertTrue(os.path.exists(second.manifest_path))
            self.assertEqual(
                [],
                [
                    name for name in os.listdir(manifest_dir)
                    if name.endswith(".tar.gz")
                ],
            )

    def test_cli_man_creates_manifest_without_archive(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)

            old_cwd = os.getcwd()
            out = io.StringIO()
            try:
                os.chdir(tmpdir)
                with contextlib.redirect_stdout(out):
                    rc = main([
                        "Jarvis",
                        "project",
                        "pack",
                        project_root,
                        "--repro",
                        "--man",
                    ])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 0, msg=out.getvalue())
            manifests = [
                name for name in os.listdir(tmpdir)
                if name.startswith("pack_") and name.endswith(".yaml")
            ]
            archives = [
                name for name in os.listdir(tmpdir)
                if name.endswith(".tar.gz")
            ]
            self.assertEqual(1, len(manifests), msg=out.getvalue())
            self.assertEqual([], archives, msg=out.getvalue())

    def test_generated_manifest_has_required_structure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)

            report = create_project_pack_manifest(
                project_root,
                profile="share",
                manifest_dir=tmpdir,
            )
            with open(report.manifest_path, "r", encoding="utf-8") as f1:
                payload = yaml.safe_load(f1)

            self.assertEqual(report.pack_id, payload["pack_id"])
            self.assertEqual("share", payload["mode"])
            self.assertEqual(
                os.path.realpath(project_root),
                os.path.realpath(payload["project_root"]),
            )
            self.assertTrue(payload["output"].endswith(".tar.gz"))
            self.assertIsInstance(payload["include"], list)
            self.assertIsInstance(payload["exclude"], list)
            self.assertIn("jarvis.project.yaml", payload["include"])

    def test_manifest_input_creates_archive_and_preserves_relative_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)
            archive_path = os.path.join(tmpdir, "from_manifest.tar.gz")
            manifest_path = os.path.join(tmpdir, "pack_20260508_001.yaml")
            payload = {
                "pack_id": "pack_20260508_001",
                "mode": "repro",
                "project_root": project_root,
                "output": archive_path,
                "include": [
                    "jarvis.project.yaml",
                    "README.md",
                    "bin/task.yaml",
                    "logs/run.log",
                    "checkpoints/state.ckpt",
                ],
                "exclude": [
                    "logs/",
                    "checkpoints/state.ckpt",
                ],
            }
            with open(manifest_path, "w", encoding="utf-8") as f1:
                yaml.safe_dump(payload, f1, sort_keys=False)

            report = create_project_package_from_manifest(manifest_path)
            self.assertEqual(archive_path, report.archive_path)
            self.assertTrue(os.path.exists(archive_path))

            names = self._tar_names(archive_path)
            self.assertIn("jarvis.project.yaml", names)
            self.assertIn("README.md", names)
            self.assertIn("bin/task.yaml", names)
            self.assertNotIn("logs/run.log", names)
            self.assertNotIn("checkpoints/state.ckpt", names)
            self.assertNotIn("DemoProject/bin/task.yaml", names)

    def test_cli_manifest_input_creates_archive(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)
            manifest_path = os.path.join(tmpdir, "pack_20260508_001.yaml")
            archive_path = os.path.join(tmpdir, "cli_manifest.tar.gz")
            payload = {
                "pack_id": "pack_20260508_001",
                "mode": "share",
                "project_root": project_root,
                "output": archive_path,
                "include": ["README.md", "bin/task.yaml"],
                "exclude": [],
            }
            with open(manifest_path, "w", encoding="utf-8") as f1:
                yaml.safe_dump(payload, f1, sort_keys=False)

            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                rc = main(["Jarvis", "project", "pack", manifest_path])

            self.assertEqual(rc, 0, msg=out.getvalue())
            self.assertTrue(os.path.exists(archive_path), msg=out.getvalue())

    def test_manifest_output_relative_to_current_working_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)
            manifest_path = os.path.join(tmpdir, "pack_20260508_001.yaml")
            payload = {
                "pack_id": "pack_20260508_001",
                "mode": "share",
                "project_root": project_root,
                "output": "relative_manifest.tar.gz",
                "include": ["README.md"],
                "exclude": [],
            }
            with open(manifest_path, "w", encoding="utf-8") as f1:
                yaml.safe_dump(payload, f1, sort_keys=False)

            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                report = create_project_package_from_manifest(manifest_path)
            finally:
                os.chdir(old_cwd)

            self.assertEqual(
                os.path.realpath(os.path.join(tmpdir, "relative_manifest.tar.gz")),
                os.path.realpath(report.archive_path),
            )
            self.assertTrue(os.path.exists(report.archive_path))

    def test_manifest_rejects_absolute_and_parent_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)
            manifest_path = os.path.join(tmpdir, "pack_20260508_001.yaml")
            base_payload = {
                "pack_id": "pack_20260508_001",
                "mode": "share",
                "project_root": project_root,
                "output": os.path.join(tmpdir, "bad.tar.gz"),
                "include": ["README.md"],
                "exclude": [],
            }

            for key, bad_path in (
                ("include", "/tmp/secret.txt"),
                ("include", "../secret.txt"),
                ("exclude", "/tmp/secret.txt"),
                ("exclude", "../secret.txt"),
            ):
                with self.subTest(key=key, bad_path=bad_path):
                    payload = dict(base_payload)
                    payload["include"] = list(base_payload["include"])
                    payload["exclude"] = list(base_payload["exclude"])
                    payload[key] = [bad_path]
                    with open(manifest_path, "w", encoding="utf-8") as f1:
                        yaml.safe_dump(payload, f1, sort_keys=False)
                    with self.assertRaises(ProjectPackError):
                        create_project_package_from_manifest(manifest_path)

    def test_manifest_fails_loudly_for_missing_included_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)
            manifest_path = os.path.join(tmpdir, "pack_20260508_001.yaml")
            payload = {
                "pack_id": "pack_20260508_001",
                "mode": "share",
                "project_root": project_root,
                "output": os.path.join(tmpdir, "missing.tar.gz"),
                "include": ["missing.txt"],
                "exclude": [],
            }
            with open(manifest_path, "w", encoding="utf-8") as f1:
                yaml.safe_dump(payload, f1, sort_keys=False)

            with self.assertRaises(ProjectPackError):
                create_project_package_from_manifest(manifest_path)


if __name__ == "__main__":
    unittest.main()
