#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import glob
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest import mock


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


class CliEntrypointTests(unittest.TestCase):
    @staticmethod
    def _build_argparser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="Jarvis Program Help Center",
            formatter_class=argparse.RawTextHelpFormatter,
        )
        args_cfg = os.path.join(PROJECT_ROOT, "jarvishep", "card", "argparser.json")
        with open(args_cfg, "r", encoding="utf-8") as f1:
            config = json.load(f1)

        for pos_arg in config.get("positionals", []):
            kwargs = {"help": pos_arg["help"]}
            if pos_arg["name"] == "file":
                kwargs["nargs"] = "?"
            parser.add_argument(pos_arg["name"], **kwargs)

        for opt in config.get("options", []):
            kwargs = {
                "help": opt["help"],
                "action": opt.get("action", "store"),
                "dest": opt.get("dest"),
            }
            if "default" in opt:
                kwargs["default"] = opt["default"]
            if "nargs" in opt:
                kwargs["nargs"] = opt["nargs"]
            if "const" in opt:
                kwargs["const"] = opt["const"]
            if "choices" in opt:
                kwargs["choices"] = opt["choices"]
            if "type" in opt:
                if opt["type"] == "int":
                    kwargs["type"] = int
                elif opt["type"] == "float":
                    kwargs["type"] = float
                else:
                    kwargs["type"] = str
            if "short" in opt and "long" in opt:
                parser.add_argument(opt["short"], opt["long"], **kwargs)
            elif "long" in opt:
                parser.add_argument(opt["long"], **kwargs)
        return parser

    @staticmethod
    def _make_demo_project(project_root: str) -> None:
        os.makedirs(os.path.join(project_root, "bin"), exist_ok=True)
        os.makedirs(os.path.join(project_root, "data"), exist_ok=True)
        os.makedirs(os.path.join(project_root, "deps"), exist_ok=True)
        with open(os.path.join(project_root, ".jarvis-project.json"), "w", encoding="utf-8") as f1:
            f1.write("{}\n")
        with open(os.path.join(project_root, "jarvis.project.yaml"), "w", encoding="utf-8") as f1:
            f1.write("project: {}\n")
        with open(os.path.join(project_root, "README.md"), "w", encoding="utf-8") as f1:
            f1.write("# Demo\n")
        with open(os.path.join(project_root, "bin", "task.yaml"), "w", encoding="utf-8") as f1:
            f1.write("Scan: {name: demo, save_dir: '&J/outputs'}\n")
        with open(os.path.join(project_root, "deps", "pkg.txt"), "w", encoding="utf-8") as f1:
            f1.write("dep\n")

    @staticmethod
    def _write_mock_official_library(tmpdir: str) -> tuple[str, str]:
        project_name = "Example_Bridson"
        source_root = os.path.join(tmpdir, "source", project_name)
        os.makedirs(os.path.join(source_root, "bin"), exist_ok=True)
        with open(os.path.join(source_root, ".jarvis-project.json"), "w", encoding="utf-8") as f1:
            f1.write("{}\n")
        with open(os.path.join(source_root, "jarvis.project.yaml"), "w", encoding="utf-8") as f1:
            f1.write("project: {}\n")
        with open(
            os.path.join(source_root, "bin", "Example_Bridson_Operas.yaml"),
            "w",
            encoding="utf-8",
        ) as f1:
            f1.write("Scan: {name: demo, save_dir: '&J/outputs'}\n")

        archive_path = os.path.join(tmpdir, "Example_Bridson.tar.gz")
        with tarfile.open(archive_path, "w:gz") as tf:
            tf.add(source_root, arcname=f"mock-official-library/{project_name}")

        library_path = os.path.join(tmpdir, "official_project_library.json")
        payload = {
            "library_name": "official Jarvis library",
            "projects": [
                {
                    "name": project_name,
                    "category": "sampling",
                    "summary": "Mock Bridson project.",
                    "entrypoint": "bin/Example_Bridson_Operas.yaml",
                    "archive_url": Path(archive_path).as_uri(),
                    "archive_root": f"mock-official-library/{project_name}",
                    "compatibility_notes": "Mock official Jarvis library entry for tests.",
                }
            ],
        }
        with open(library_path, "w", encoding="utf-8") as f1:
            json.dump(payload, f1, indent=2, ensure_ascii=False)
            f1.write("\n")

        return Path(library_path).as_uri(), project_name

    def test_import_jarvishep_client_main(self):
        from jarvishep.client import main

        self.assertTrue(callable(main))

    def test_pyproject_jarvis_entrypoint(self):
        pyproject = os.path.join(PROJECT_ROOT, "pyproject.toml")
        with open(pyproject, "r", encoding="utf-8") as f1:
            content = f1.read()
        self.assertIn('[project.scripts]', content)
        self.assertIn('Jarvis = "jarvishep.client:main"', content)

    def test_help_contract_parser_excludes_removed_flags(self):
        parser = self._build_argparser()
        help_text = parser.format_help()
        self.assertIn("Jarvis Program Help Center", help_text)
        self.assertIn("--plot", help_text)
        self.assertIn("--convert", help_text)
        self.assertIn("--monitor", help_text)
        self.assertIn("--resume", help_text)
        self.assertIn("--check-modules", help_text)
        self.assertNotIn("--mkproject", help_text)
        self.assertNotIn("--packproject", help_text)
        self.assertNotIn("--profile", help_text)
        self.assertNotIn("--max-concurrency", help_text)
        self.assertNotIn("--per-task-timeout-sec", help_text)
        self.assertNotIn("--progress-interval-sec", help_text)
        self.assertNotIn("--log-policy", help_text)

    def test_removed_flags_are_rejected_by_argparser(self):
        parser = self._build_argparser()
        removed = [
            ["demo.yaml", "--mkproject", "Demo"],
            ["demo.yaml", "--packproject"],
            ["demo.yaml", "--profile", "repro"],
            ["demo.yaml", "--max-concurrency", "4"],
            ["demo.yaml", "--per-task-timeout-sec", "30"],
            ["demo.yaml", "--progress-interval-sec", "2"],
            ["demo.yaml", "--log-policy", "quiet"],
        ]
        for argv in removed:
            with self.subTest(argv=argv):
                err = io.StringIO()
                with contextlib.redirect_stderr(err):
                    with self.assertRaises(SystemExit) as cm:
                        parser.parse_args(argv)
                self.assertEqual(cm.exception.code, 2)
                self.assertIn("unrecognized arguments", err.getvalue())

    def test_main_top_level_help_is_clean_and_entry_oriented(self):
        from jarvishep.client import main

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["Jarvis", "-h"])

        self.assertEqual(rc, 0)
        text = out.getvalue()
        self.assertIn("Usage:", text)
        self.assertIn("Jarvis [file] [options]", text)
        self.assertIn("Jarvis project <command> [arguments]", text)
        self.assertIn("Jarvis portal formats", text)
        self.assertIn("Main entry points:", text)
        self.assertIn("General options:", text)
        self.assertIn("Workflow options:", text)
        self.assertIn("Run `Jarvis project -h`", text)
        self.assertIn("Run `Jarvis portal -h`", text)
        self.assertNotIn("--mkproject", text)
        self.assertNotIn("--packproject", text)
        self.assertNotIn("--profile", text)
        self.assertNotIn("--max-concurrency", text)

    def test_main_removed_legacy_project_flags_fail(self):
        from jarvishep.client import main

        removed = [
            ["Jarvis", "--mkproject", "DemoProject"],
            ["Jarvis", "--packproject"],
            ["Jarvis", "--profile", "repro"],
        ]
        for argv in removed:
            with self.subTest(argv=argv):
                err = io.StringIO()
                with contextlib.redirect_stderr(err):
                    with self.assertRaises(SystemExit) as cm:
                        main(argv)
                self.assertEqual(cm.exception.code, 2)
                self.assertIn("unrecognized arguments", err.getvalue())

    def test_main_major_workflow_mode_flags_still_parse(self):
        parser = self._build_argparser()
        cases = [
            ("--plot", "plot"),
            ("--convert", "cvtDB"),
            ("--monitor", "monitor"),
            ("--resume", "resume"),
            ("--check-modules", "OPC"),
        ]
        for flag, dest in cases:
            with self.subTest(flag=flag):
                args = parser.parse_args(["demo.yaml", flag])
                self.assertTrue(getattr(args, dest))

    def test_main_project_help_pages_match(self):
        from jarvishep.client import main

        outputs = {}
        for argv in (
            ["Jarvis", "project"],
            ["Jarvis", "project", "-h"],
            ["Jarvis", "project", "--help"],
        ):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                rc = main(argv)
            self.assertEqual(rc, 0)
            outputs[tuple(argv)] = out.getvalue()

        self.assertEqual(outputs[("Jarvis", "project")], outputs[("Jarvis", "project", "-h")])
        self.assertEqual(outputs[("Jarvis", "project")], outputs[("Jarvis", "project", "--help")])
        text = outputs[("Jarvis", "project")]
        self.assertIn("Jarvis project <command> [arguments]", text)
        self.assertIn("Pack a local project for sharing, reproduction, or full export", text)

    def test_main_portal_help_pages_match(self):
        from jarvishep.client import main

        outputs = {}
        for argv in (
            ["Jarvis", "portal"],
            ["Jarvis", "portal", "-h"],
            ["Jarvis", "portal", "--help"],
        ):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                rc = main(argv)
            self.assertEqual(rc, 0)
            outputs[tuple(argv)] = out.getvalue()

        self.assertEqual(outputs[("Jarvis", "portal")], outputs[("Jarvis", "portal", "-h")])
        self.assertEqual(outputs[("Jarvis", "portal")], outputs[("Jarvis", "portal", "--help")])
        text = outputs[("Jarvis", "portal")]
        self.assertIn("Jarvis portal formats", text)
        self.assertIn("Show calculator IO formats", text)

    def test_main_portal_formats_lists_jarvis_portal_registry(self):
        from jarvishep.client import main

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["Jarvis", "portal", "formats"])

        self.assertEqual(rc, 0)
        text = out.getvalue()
        self.assertIn("Jarvis-Portal calculator IO formats", text)
        self.assertIn("Format", text)
        self.assertIn("All", text)
        self.assertIn("Input", text)
        self.assertIn("Output", text)
        self.assertIn("✓", text)
        self.assertNotIn("Formats", text)
        self.assertNotIn("| Scope", text)
        self.assertIn("JSON", text)

    def test_main_portal_formats_rejects_extra_args(self):
        from jarvishep.client import main

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["Jarvis", "portal", "formats", "extra"])

        self.assertEqual(rc, 2)
        self.assertIn("Usage error: Jarvis portal formats", out.getvalue())

    def test_main_project_pack_help(self):
        from jarvishep.client import main

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["Jarvis", "project", "pack", "-h"])

        self.assertEqual(rc, 0)
        text = out.getvalue()
        self.assertIn("Jarvis project pack [path] [--share | --repro | --full] [--man]", text)
        self.assertIn("Jarvis project pack <pack_manifest.yaml>", text)
        self.assertIn("--share", text)
        self.assertIn("--repro", text)
        self.assertIn("--full", text)
        self.assertIn("--man", text)
        self.assertIn("--share` is used by default", text)

    def test_main_project_subcommand_help_pages(self):
        from jarvishep.client import main

        commands = ("create", "browse", "fetch", "info")
        for command in commands:
            with self.subTest(command=command):
                out = io.StringIO()
                with contextlib.redirect_stdout(out):
                    rc = main(["Jarvis", "project", command, "--help"])
                self.assertEqual(rc, 0)
                text = out.getvalue()
                self.assertIn("Usage:", text)
                self.assertIn(f"Jarvis project {command}", text)

    def test_main_project_create_success(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                rc = main(["Jarvis", "project", "create", "NewProject"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 0)
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "NewProject", "bin")))
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "NewProject", "data")))
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "NewProject", "deps")))
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "NewProject", "outputs")))

    def test_main_project_pack_modes(self):
        from jarvishep.client import main

        for mode_flag, profile in (("--share", "share"), ("--repro", "repro"), ("--full", "full")):
            with self.subTest(mode_flag=mode_flag):
                with tempfile.TemporaryDirectory() as tmpdir:
                    project_root = os.path.join(tmpdir, "DemoProject")
                    self._make_demo_project(project_root)

                    old_cwd = os.getcwd()
                    out = io.StringIO()
                    try:
                        os.chdir(project_root)
                        with contextlib.redirect_stdout(out):
                            rc = main(["Jarvis", "project", "pack", ".", mode_flag])
                    finally:
                        os.chdir(old_cwd)

                    self.assertEqual(rc, 0, msg=out.getvalue())
                    archives = glob.glob(os.path.join(tmpdir, f"DemoProject_{profile}_*.tar.gz"))
                    self.assertEqual(len(archives), 1, msg=out.getvalue())

    def test_main_project_pack_defaults_to_share(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)

            old_cwd = os.getcwd()
            out = io.StringIO()
            try:
                os.chdir(project_root)
                with contextlib.redirect_stdout(out):
                    rc = main(["Jarvis", "project", "pack"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 0, msg=out.getvalue())
            archives = glob.glob(os.path.join(tmpdir, "DemoProject_share_*.tar.gz"))
            self.assertEqual(len(archives), 1, msg=out.getvalue())

    def test_main_project_pack_rejects_multiple_modes(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            self._make_demo_project(project_root)

            old_cwd = os.getcwd()
            out = io.StringIO()
            try:
                os.chdir(project_root)
                with contextlib.redirect_stdout(out):
                    rc = main(["Jarvis", "project", "pack", ".", "--share", "--repro"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 2)
            self.assertIn("mutually exclusive", out.getvalue())

    def test_main_project_browse_info_and_fetch_with_mock_library(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            library_url, project_name = self._write_mock_official_library(tmpdir)
            env_patch = {
                "JARVIS_OFFICIAL_LIBRARY_INDEX_URL": library_url,
            }

            browse_out = io.StringIO()
            with mock.patch.dict(os.environ, env_patch, clear=False):
                with contextlib.redirect_stdout(browse_out):
                    rc_browse = main(["Jarvis", "project", "browse"])
            self.assertEqual(rc_browse, 0)
            self.assertIn(project_name, browse_out.getvalue())

            info_out = io.StringIO()
            with mock.patch.dict(os.environ, env_patch, clear=False):
                with contextlib.redirect_stdout(info_out):
                    rc_info = main(["Jarvis", "project", "info", project_name])
            self.assertEqual(rc_info, 0)
            self.assertIn("Official project", info_out.getvalue())
            self.assertIn("Entrypoint: bin/Example_Bridson_Operas.yaml", info_out.getvalue())

            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                fetch_out = io.StringIO()
                with mock.patch.dict(os.environ, env_patch, clear=False):
                    with contextlib.redirect_stdout(fetch_out):
                        rc_fetch = main(["Jarvis", "project", "fetch", project_name])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc_fetch, 0, msg=fetch_out.getvalue())
            fetched_root = os.path.join(tmpdir, project_name)
            self.assertTrue(os.path.isdir(fetched_root))
            self.assertTrue(
                os.path.isfile(
                    os.path.join(fetched_root, "bin", "Example_Bridson_Operas.yaml")
                )
            )

    def test_main_project_fetch_cleans_partial_target_on_failure(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            broken_archive = os.path.join(tmpdir, "broken.tar.gz")
            with open(broken_archive, "wb") as f1:
                f1.write(b"not-an-archive")

            library_path = os.path.join(tmpdir, "official_project_library.json")
            payload = {
                "library_name": "official Jarvis library",
                "projects": [
                    {
                        "name": "BrokenProject",
                        "entrypoint": "bin/task.yaml",
                        "archive_url": Path(broken_archive).as_uri(),
                        "archive_root": "BrokenProject",
                    }
                ],
            }
            with open(library_path, "w", encoding="utf-8") as f1:
                json.dump(payload, f1, indent=2, ensure_ascii=False)
                f1.write("\n")

            env_patch = {"JARVIS_OFFICIAL_LIBRARY_INDEX_URL": Path(library_path).as_uri()}
            old_cwd = os.getcwd()
            out = io.StringIO()
            try:
                os.chdir(tmpdir)
                with mock.patch.dict(os.environ, env_patch, clear=False):
                    with contextlib.redirect_stdout(out):
                        rc = main(["Jarvis", "project", "fetch", "BrokenProject"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 1)
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "BrokenProject")))

    def test_main_missing_yaml_file_exits_cleanly(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            err = io.StringIO()
            try:
                os.chdir(tmpdir)
                with contextlib.redirect_stderr(err):
                    with self.assertRaises(SystemExit) as cm:
                        main(["Jarvis", "./bin/Example_Bridson_Operas.yaml"])
            finally:
                os.chdir(old_cwd)

        self.assertEqual(cm.exception.code, 2)
        text = err.getvalue()
        self.assertIn("YAML file not found:", text)
        self.assertNotIn("Traceback", text)

    def test_python_module_project_create_success(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            env = dict(os.environ)
            pythonpath = env.get("PYTHONPATH", "")
            env["PYTHONPATH"] = PROJECT_ROOT if not pythonpath else f"{PROJECT_ROOT}{os.pathsep}{pythonpath}"
            proc = subprocess.run(
                [sys.executable, "-m", "jarvishep", "project", "create", "SubprocProject"],
                cwd=tmpdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stdout)
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "SubprocProject", "bin")))

    def test_main_version_fast_path(self):
        from jarvishep.client import main

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["Jarvis", "--version"])

        self.assertEqual(rc, 0)
        text = out.getvalue()
        self.assertIn("Author:", text)
        self.assertIn("Version:", text)
        self.assertIn("Resources:", text)
        self.assertIn("References:", text)


if __name__ == "__main__":
    unittest.main()
