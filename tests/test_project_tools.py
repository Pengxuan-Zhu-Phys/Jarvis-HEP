#!/usr/bin/env python3
"""D12.5 / D12.6: project tools + Examples catalog + restricted key fetch."""

from __future__ import annotations

import io
import json
import os
import tarfile
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
from unittest import mock

from jarvishep2.client import dispatch_project, main
from jarvishep2.official_project_library import (
    DEFAULT_OFFICIAL_LIBRARY_INDEX_URL,
    OfficialCatalogError,
    OfficialProjectFetchError,
    OfficialProjectNotFoundError,
    fetch_official_project,
    format_project_list_table,
    get_official_project,
    list_official_projects,
    load_official_project_catalog,
)
from jarvishep2.project_crypto import (
    encrypt_file,
    is_openssl_encrypted_file,
    maybe_decrypt_archive,
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
            for name in (
                "quickstart_random.yaml",
                "quickstart_random_operas.yaml",
                "quickstart_random_calculator.yaml",
            ):
                self.assertTrue(os.path.isfile(os.path.join(root, "bin", name)))
            defaults = os.path.join(root, "deps", "environment_default.yaml")
            self.assertTrue(os.path.isfile(defaults))

            calculator = os.path.join(root, "bin", "quickstart_calculator.yaml")
            with open(calculator, encoding="utf-8") as handle:
                self.assertNotIn("calculators/echo", handle.read())

            config = load_task_yaml(quick)
            self.assertEqual(config["Runtime"]["mode"], "redis")
            self.assertIn("workers", config["Runtime"])
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
    def test_archive_names_use_hour_stamp_and_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = create_project_scaffold("NamedProject", cwd=tmp)
            stamp = datetime.now(timezone.utc).strftime("%m%d%H")

            first = create_project_package(root, profile="share")
            second = create_project_package(root, profile="share")

            self.assertEqual(
                os.path.basename(first.archive_path),
                f"NamedProject_share_{stamp}_01.tar.gz",
            )
            self.assertEqual(
                os.path.basename(second.archive_path),
                f"NamedProject_share_{stamp}_02.tar.gz",
            )

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

    def test_share_pack_includes_deps_unless_project_excludes_them(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = create_project_scaffold("DependencyPack", cwd=tmp)
            included = os.path.join(root, "deps", "runtime", "source.txt")
            excluded = os.path.join(root, "deps", "private", "weights.bin")
            calculator = os.path.join(root, "calculators", "cached", "main")
            log = os.path.join(root, "logs", "run.log")
            for path, content in (
                (included, "required dependency"),
                (excluded, "private dependency"),
                (calculator, "local build"),
                (log, "local log"),
            ):
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write(content)
            with open(os.path.join(root, "jarvis.project.yaml"), "a", encoding="utf-8") as handle:
                handle.write("\npack:\n  exclude:\n    - deps/private/\n")

            report = create_project_package(root, profile="share")
            with tarfile.open(report.archive_path, "r:gz") as tar:
                names = tar.getnames()

            self.assertIn("DependencyPack/deps/runtime/source.txt", names)
            self.assertNotIn("DependencyPack/deps/private/weights.bin", names)
            self.assertNotIn("DependencyPack/calculators/cached/main", names)
            self.assertNotIn("DependencyPack/logs/run.log", names)

    def test_project_pack_excludes_are_local_and_apply_to_all_profiles(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = create_project_scaffold("ExcludePrivate", cwd=tmp)
            private_dir = os.path.join(root, "deps", "private")
            os.makedirs(private_dir, exist_ok=True)
            with open(os.path.join(private_dir, "weights.bin"), "wb") as handle:
                handle.write(b"private weights")
            with open(os.path.join(root, "jarvis.project.yaml"), "a", encoding="utf-8") as handle:
                handle.write("\npack:\n  exclude:\n    - deps/private/\n")

            for profile in ("share", "repro", "full"):
                report = create_project_package(root, profile=profile)
                with tarfile.open(report.archive_path, "r:gz") as tar:
                    names = tar.getnames()
                self.assertFalse(
                    any("deps/private/" in name for name in names),
                    f"deps/private was included in {profile} package",
                )

            other_root = create_project_scaffold("KeepPrivate", cwd=tmp)
            other_private = os.path.join(other_root, "deps", "private", "weights.bin")
            os.makedirs(os.path.dirname(other_private), exist_ok=True)
            with open(other_private, "wb") as handle:
                handle.write(b"keep this in an unrelated project")
            report = create_project_package(other_root, profile="full")
            with tarfile.open(report.archive_path, "r:gz") as tar:
                self.assertIn("KeepPrivate/deps/private/weights.bin", tar.getnames())

    def test_share_pack_can_allowlist_files_under_an_excluded_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = create_project_scaffold("BinFilter", cwd=tmp)
            for name in (
                "adaptive_alpha.yaml",
                "adaptive_beta.yaml",
                "other.yaml",
            ):
                with open(os.path.join(root, "bin", name), "w", encoding="utf-8") as handle:
                    handle.write("Scan: {}\n")
            with open(os.path.join(root, "jarvis.project.yaml"), "a", encoding="utf-8") as handle:
                handle.write(
                    "\npack:\n"
                    "  share:\n"
                    "    exclude:\n"
                    "      - bin/\n"
                    "    include:\n"
                    "      - bin/adaptive_alpha.yaml\n"
                    "      - bin/adaptive_beta.yaml\n"
                )

            report = create_project_package(root, profile="share")
            with tarfile.open(report.archive_path, "r:gz") as tar:
                names = tar.getnames()

            bin_names = {
                name.split("/", 1)[1]
                for name in names
                if "/bin/" in name
            }
            self.assertEqual(
                bin_names,
                {
                    "bin/adaptive_alpha.yaml",
                    "bin/adaptive_beta.yaml",
                },
            )


class OfficialLibraryTests(unittest.TestCase):
    def test_default_index_points_at_examples_github(self) -> None:
        self.assertIn("Jarvis-Examples", DEFAULT_OFFICIAL_LIBRARY_INDEX_URL)
        self.assertIn("catalog/official_project_library.json", DEFAULT_OFFICIAL_LIBRARY_INDEX_URL)

    def test_packaged_catalog_browse_and_info(self) -> None:
        with mock.patch(
            "jarvishep2.official_project_library._read_catalog_from_url",
            side_effect=OfficialCatalogError("offline"),
        ):
            projects = list_official_projects()
        self.assertTrue(projects)
        names = {p["name"] for p in projects}
        self.assertIn("Eggbox", names)
        egg = next(p for p in projects if p["name"] == "Eggbox")
        self.assertEqual(egg["access"], "public")
        self.assertFalse(egg["requires_key"])

        table = format_project_list_table(projects)
        self.assertIn("Eggbox", table)
        self.assertIn("public", table)
        self.assertIn("no", table)  # Key column

        with mock.patch(
            "jarvishep2.official_project_library._read_catalog_from_url",
            side_effect=OfficialCatalogError("offline"),
        ):
            egg = get_official_project("Eggbox")
        self.assertEqual(egg["name"], "Eggbox")

        with self.assertRaises(OfficialProjectNotFoundError):
            with mock.patch(
                "jarvishep2.official_project_library._read_catalog_from_url",
                side_effect=OfficialCatalogError("offline"),
            ):
                get_official_project("DefinitelyMissingProjectXYZ")

    def test_restricted_entry_normalized(self) -> None:
        payload = {
            "schema_version": 1,
            "projects": [
                {
                    "name": "SecretToy",
                    "access": "restricted",
                    "requires_key": True,
                    "encryption": {
                        "scheme": "openssl-aes-256-cbc",
                        "hint": "ask Alice",
                    },
                    "archive_url": "https://example.test/secret.jenc",
                    "summary": "private",
                },
                {
                    "name": "OpenToy",
                    "access": "public",
                    "archive_url": "https://example.test/open.tar.gz",
                },
            ],
        }
        with mock.patch(
            "jarvishep2.official_project_library._load_catalog_payload",
            return_value=payload,
        ):
            projects = list_official_projects()
        by_name = {p["name"]: p for p in projects}
        self.assertTrue(by_name["SecretToy"]["requires_key"])
        self.assertEqual(by_name["SecretToy"]["access"], "restricted")
        self.assertEqual(by_name["SecretToy"]["encryption_hint"], "ask Alice")
        self.assertFalse(by_name["OpenToy"]["requires_key"])
        table = format_project_list_table(projects)
        self.assertIn("required", table)
        self.assertIn("SecretToy", table)

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

    def test_cli_browse_and_removed_list_command(self) -> None:
        with mock.patch(
            "jarvishep2.official_project_library._read_catalog_from_url",
            side_effect=OfficialCatalogError("offline"),
        ):
            self.assertEqual(dispatch_project(["browse"]), 0)
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                self.assertNotEqual(dispatch_project(["list"]), 0)

    def test_cli_browse_uses_rich_project_panel(self) -> None:
        output = io.StringIO()
        with mock.patch(
            "jarvishep2.official_project_library._read_catalog_from_url",
            side_effect=OfficialCatalogError("offline"),
        ), redirect_stdout(output):
            self.assertEqual(dispatch_project(["browse"]), 0)

        rendered = output.getvalue()
        self.assertIn("Official library", rendered)
        self.assertIn("NAME", rendered)
        self.assertIn("CATEGORY", rendered)
        self.assertIn("Key: required", rendered)
        self.assertIn("Jarvis project fetch Eggbox", rendered)
        self.assertIn("Jarvis project fetch iDM --key 'YOUR_KEY'", rendered)

    def test_cli_info_uses_rich_project_card(self) -> None:
        output = io.StringIO()
        with mock.patch(
            "jarvishep2.official_project_library._read_catalog_from_url",
            side_effect=OfficialCatalogError("offline"),
        ), redirect_stdout(output):
            self.assertEqual(dispatch_project(["info", "Eggbox"]), 0)

        rendered = output.getvalue()
        self.assertIn("Official project", rendered)
        self.assertIn("Eggbox", rendered)
        self.assertIn("KEY REQUIRED", rendered)
        self.assertIn("COMPATIBILITY", rendered)
        self.assertIn("bin/Example_Bridson_Operas.yaml", rendered)

    def test_main_project_passthrough(self) -> None:
        code = main(["project", "-h"])
        self.assertEqual(code, 0)

    def test_fetch_restricted_without_key_fails(self) -> None:
        payload = {
            "schema_version": 1,
            "projects": [
                {
                    "name": "Locked",
                    "access": "restricted",
                    "requires_key": True,
                    "archive_url": "https://example.test/x.jenc",
                    "encryption": {"scheme": "openssl-aes-256-cbc"},
                }
            ],
        }
        with mock.patch(
            "jarvishep2.official_project_library._load_catalog_payload",
            return_value=payload,
        ):
            with mock.patch.dict(os.environ, {}, clear=False):
                os.environ.pop("JARVIS_PROJECT_FETCH_KEY", None)
                with self.assertRaises(OfficialProjectFetchError):
                    fetch_official_project("Locked")


class ProjectCryptoTests(unittest.TestCase):
    def test_cli_encrypt_and_pack_encrypt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = create_project_scaffold("EncMe", cwd=tmp)
            cwd = os.getcwd()
            try:
                os.chdir(tmp)
                code = dispatch_project(
                    ["pack", "EncMe", "--share", "--encrypt", "--key", "cli-secret"]
                )
            finally:
                os.chdir(cwd)
            self.assertEqual(code, 0)
            jencs = [
                os.path.join(tmp, name)
                for name in os.listdir(tmp)
                if name.endswith(".jenc")
            ]
            self.assertTrue(jencs, "expected encrypted pack next to project parent")
            self.assertTrue(is_openssl_encrypted_file(jencs[0]))

            # encrypt subcommand on a plain file
            plain = os.path.join(tmp, "manual.bin")
            with open(plain, "wb") as handle:
                handle.write(b"payload-for-encrypt-cmd")
            code = dispatch_project(["encrypt", plain, "--key", "cli-secret"])
            self.assertEqual(code, 0)
            self.assertTrue(os.path.isfile(plain + ".jenc"))

    def test_openssl_encrypt_decrypt_roundtrip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            plain = os.path.join(tmp, "payload.tar.gz")
            enc = os.path.join(tmp, "payload.tar.gz.jenc")
            # Minimal gzip-ish content is fine; openssl only needs bytes.
            with open(plain, "wb") as handle:
                handle.write(b"hello-jarvis-secret-payload-0123456789")
            encrypt_file(plain, enc, key="unit-test-key")
            self.assertTrue(is_openssl_encrypted_file(enc))
            out = maybe_decrypt_archive(
                enc,
                requires_key=True,
                encryption_scheme="openssl-aes-256-cbc",
                key="unit-test-key",
            )
            self.assertNotEqual(out, enc)
            with open(out, "rb") as handle:
                self.assertEqual(handle.read(), b"hello-jarvis-secret-payload-0123456789")
            os.unlink(out)

    def test_fetch_encrypted_archive_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            # Build a tiny project tarball then encrypt it.
            proj = os.path.join(tmp, "Mini")
            os.makedirs(os.path.join(proj, "bin"))
            with open(os.path.join(proj, ".jarvis-project.json"), "w", encoding="utf-8") as handle:
                json.dump({"format": "jarvis-hep-standalone-project", "version": 1}, handle)
            with open(os.path.join(proj, "bin", "run.yaml"), "w", encoding="utf-8") as handle:
                handle.write("Scan:\n  name: mini\n")
            tar_path = os.path.join(tmp, "mini.tar.gz")
            with tarfile.open(tar_path, "w:gz") as tar:
                tar.add(proj, arcname="Mini")
            enc_path = os.path.join(tmp, "mini.tar.gz.jenc")
            encrypt_file(tar_path, enc_path, key="secret-mini")

            payload = {
                "schema_version": 1,
                "projects": [
                    {
                        "name": "Mini",
                        "access": "restricted",
                        "requires_key": True,
                        "archive_url": f"file://{enc_path}",
                        "archive_root": ".",
                        "entrypoint": "bin/run.yaml",
                        "encryption": {"scheme": "openssl-aes-256-cbc"},
                    }
                ],
            }
            dest = os.path.join(tmp, "out", "Mini")
            with mock.patch(
                "jarvishep2.official_project_library._load_catalog_payload",
                return_value=payload,
            ):
                report = fetch_official_project(
                    "Mini",
                    target_dir=dest,
                    key="secret-mini",
                )
            self.assertEqual(report.project_name, "Mini")
            self.assertTrue(os.path.isfile(os.path.join(dest, "bin", "run.yaml")))
            self.assertTrue(report.required_key)


if __name__ == "__main__":
    unittest.main()
