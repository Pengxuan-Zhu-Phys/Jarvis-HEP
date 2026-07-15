#!/usr/bin/env python3
"""D12.5: Jarvis2 project create/pack/browse/info (scaffold + library)."""

from __future__ import annotations

import json
import os
import tarfile
import tempfile
import unittest
from unittest import mock

from jarvishep2.client import dispatch_project, main
from jarvishep2.official_project_library import (
    OfficialCatalogError,
    OfficialProjectNotFoundError,
    get_official_project,
    list_official_projects,
    load_official_project_catalog,
)
from jarvishep2.project_packager import (
    create_project_package,
    create_project_pack_manifest,
)
from jarvishep2.project_scaffold import (
    PROJECT_DESCRIPTOR_NAME,
    PROJECT_MARKER_NAME,
    PROJECT_SUBDIRS,
    create_project_scaffold,
)
from jarvishep2.task_config import load_task_yaml


class ProjectScaffoldTests(unittest.TestCase):
    def test_create_scaffold_layout_and_v2_yaml(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = create_project_scaffold("DemoProject", cwd=tmp)
            self.assertTrue(os.path.isdir(root))
            for sub in PROJECT_SUBDIRS:
                self.assertTrue(os.path.isdir(os.path.join(root, sub)))
            self.assertTrue(os.path.isfile(os.path.join(root, PROJECT_MARKER_NAME)))
            self.assertTrue(os.path.isfile(os.path.join(root, PROJECT_DESCRIPTOR_NAME)))
            quick = os.path.join(root, "bin", "quickstart_bridson_operas.yaml")
            self.assertTrue(os.path.isfile(quick))
            defaults = os.path.join(root, "deps", "environment_default.yaml")
            self.assertTrue(os.path.isfile(defaults))

            config = load_task_yaml(quick)
            self.assertEqual(config["Runtime"]["mode"], "redis")
            self.assertIn("workers", config["Runtime"])
            # No forbidden top-level Runtime in source YAML (loader synthesizes it).
            with open(quick, encoding="utf-8") as handle:
                text = handle.read()
            self.assertNotIn("\nRuntime:", text)

    def test_create_refuses_existing_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            create_project_scaffold("Exists", cwd=tmp)
            with self.assertRaises(FileExistsError):
                create_project_scaffold("Exists", cwd=tmp)

    def test_cli_project_create(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cwd = os.getcwd()
            try:
                os.chdir(tmp)
                code = dispatch_project(["create", "CliDemo"])
            finally:
                os.chdir(cwd)
            self.assertEqual(code, 0)
            self.assertTrue(os.path.isdir(os.path.join(tmp, "CliDemo", "bin")))


class ProjectPackTests(unittest.TestCase):
    def test_pack_share_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = create_project_scaffold("PackMe", cwd=tmp)
            report = create_project_package(root, profile="share")
            self.assertTrue(os.path.isfile(report.archive_path))
            self.assertGreater(report.included_files, 0)
            with tarfile.open(report.archive_path, "r:gz") as tar:
                names = tar.getnames()
            self.assertTrue(any(n.endswith("jarvis.project.yaml") for n in names))
            self.assertTrue(any("quickstart_bridson_operas.yaml" in n for n in names))

            man = create_project_pack_manifest(root, profile="repro", manifest_dir=tmp)
            self.assertTrue(os.path.isfile(man.manifest_path))
            self.assertEqual(man.profile, "repro")


class OfficialLibraryTests(unittest.TestCase):
    def test_packaged_catalog_browse_and_info(self) -> None:
        # Force packaged fallback (no network).
        with mock.patch(
            "jarvishep2.official_project_library._read_catalog_from_url",
            side_effect=OfficialCatalogError("offline"),
        ):
            projects = list_official_projects()
        self.assertTrue(projects)
        names = {p["name"] for p in projects}
        self.assertIn("Eggbox", names)

        with mock.patch(
            "jarvishep2.official_project_library._read_catalog_from_url",
            side_effect=OfficialCatalogError("offline"),
        ):
            egg = get_official_project("Eggbox")
        self.assertEqual(egg["name"], "Eggbox")
        self.assertTrue(egg.get("entrypoint"))

        with self.assertRaises(OfficialProjectNotFoundError):
            with mock.patch(
                "jarvishep2.official_project_library._read_catalog_from_url",
                side_effect=OfficialCatalogError("offline"),
            ):
                get_official_project("DefinitelyMissingProjectXYZ")

    def test_schema_version_reject_future_major(self) -> None:
        payload = {
            "schema_version": 99,
            "library_name": "test",
            "projects": [{"name": "X", "entrypoint": "bin/a.yaml"}],
        }
        with mock.patch(
            "jarvishep2.official_project_library._load_catalog_payload",
            return_value=payload,
        ):
            with self.assertRaises(OfficialCatalogError):
                load_official_project_catalog()

    def test_cli_browse(self) -> None:
        with mock.patch(
            "jarvishep2.official_project_library._read_catalog_from_url",
            side_effect=OfficialCatalogError("offline"),
        ):
            code = dispatch_project(["browse"])
        self.assertEqual(code, 0)

    def test_main_project_passthrough(self) -> None:
        code = main(["project", "-h"])
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
