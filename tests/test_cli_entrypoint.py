#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import glob
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
        with open(args_cfg, "r", encoding="utf-8") as f:
            config = json.load(f)

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

    def test_import_jarvishep_client_main(self):
        from jarvishep.client import main

        self.assertTrue(callable(main))

    def test_pyproject_jarvis_entrypoint(self):
        pyproject = os.path.join(PROJECT_ROOT, "pyproject.toml")
        with open(pyproject, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn('[project.scripts]', content)
        self.assertIn('Jarvis = "jarvishep.client:main"', content)

    def test_pyproject_includes_jarvis_operas_runtime_dependency(self):
        pyproject = os.path.join(PROJECT_ROOT, "pyproject.toml")
        with open(pyproject, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn('"Jarvis-Operas>=1.3.2"', content)

    def test_help_contract_does_not_expose_install_dependencies(self):
        parser = self._build_argparser()
        help_text = parser.format_help()
        self.assertIn("Jarvis Program Help Center", help_text)
        self.assertNotIn("--install-dependencies", help_text)
        self.assertIn("--version", help_text)
        self.assertIn("--plot", help_text)
        self.assertIn("--convert", help_text)
        self.assertIn("--monitor", help_text)
        self.assertIn("--mkproject", help_text)
        self.assertIn("--packproject", help_text)
        self.assertIn("--profile", help_text)

    def test_removed_flag_is_rejected_by_argparser(self):
        parser = self._build_argparser()
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            with self.assertRaises(SystemExit) as cm:
                parser.parse_args(["demo.yaml", "--install-dependencies"])
        self.assertEqual(cm.exception.code, 2)
        self.assertIn("unrecognized arguments: --install-dependencies", err.getvalue())

    def test_major_mode_flags_still_parse(self):
        parser = self._build_argparser()
        cases = [
            ("--plot", "plot"),
            ("--convert", "cvtDB"),
            ("--monitor", "monitor"),
            ("--check-modules", "OPC"),
        ]
        for flag, dest in cases:
            with self.subTest(flag=flag):
                args = parser.parse_args(["demo.yaml", flag])
                self.assertTrue(getattr(args, dest))

    def test_main_mkproject_success(self):
        from jarvishep.client import main
        from jarvishep.project_scaffold import PROJECT_DESCRIPTOR_NAME, PROJECT_SUBDIRS

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                rc = main(["Jarvis", "--mkproject", "DemoProject"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 0)
            project_root = os.path.join(tmpdir, "DemoProject")
            self.assertTrue(os.path.isdir(project_root))
            for subdir in PROJECT_SUBDIRS:
                self.assertTrue(os.path.isdir(os.path.join(project_root, subdir)))
            self.assertTrue(os.path.isfile(os.path.join(project_root, ".jarvis-project.json")))
            self.assertTrue(os.path.isfile(os.path.join(project_root, PROJECT_DESCRIPTOR_NAME)))
            self.assertTrue(
                os.path.isfile(
                    os.path.join(project_root, "bin", "quickstart_mcmc_operas.yaml")
                )
            )
            self.assertTrue(
                os.path.isfile(
                    os.path.join(project_root, "bin", "quickstart_csv_operas.yaml")
                )
            )
            self.assertTrue(
                os.path.isfile(
                    os.path.join(project_root, "deps", "environment_default.yaml")
                )
            )

    def test_main_mkproject_conflict_flag(self):
        from jarvishep.client import main

        rc = main(["Jarvis", "--mkproject", "DemoProject", "--plot"])
        self.assertEqual(rc, 2)

    def test_main_plot_emits_plot_yaml_under_project_images_root(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                with mock.patch.dict(
                    os.environ,
                    {"JARVIS_HEP_TASK_ROOT": "", "JHEP_TASK_ROOT": ""},
                    clear=False,
                ):
                    rc_mk = main(["Jarvis", "--mkproject", "DemoProject"])
                    self.assertEqual(rc_mk, 0)

                    project_root = os.path.join(tmpdir, "DemoProject")
                    os.chdir(project_root)
                    with mock.patch(
                        "jarvishep.core.Core.init_operas_functions",
                        side_effect=AssertionError("plot mode should not initialize Jarvis-Operas"),
                    ):
                        rc_plot = main(["Jarvis", "bin/quickstart_mcmc_operas.yaml", "--plot"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc_plot, 0)
            project_root = os.path.join(tmpdir, "DemoProject")
            scan_name = "quickstart_mcmc_operas"
            self.assertTrue(
                os.path.isfile(os.path.join(project_root, "images", f"{scan_name}.yaml"))
            )
            self.assertFalse(
                os.path.exists(
                    os.path.join(project_root, "outputs", scan_name, "IMAGE", f"{scan_name}.yaml")
                )
            )

    def test_main_mkproject_rejects_removed_flag(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                rc = main(["Jarvis", "--mkproject", "DemoProject", "--install-dependencies"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 2)
            self.assertFalse(os.path.isdir(os.path.join(tmpdir, "DemoProject")))

    def test_main_mkproject_missing_name(self):
        from jarvishep.client import main

        rc = main(["Jarvis", "--mkproject"])
        self.assertEqual(rc, 2)

    def test_main_mkproject_existing_dir(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            os.makedirs(os.path.join(tmpdir, "DemoProject"), exist_ok=True)
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                rc = main(["Jarvis", "--mkproject", "DemoProject"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 1)

    def test_main_packproject_default_current_dir_success(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = os.path.join(tmpdir, "DemoProject")
            os.makedirs(os.path.join(project_root, "bin"), exist_ok=True)
            os.makedirs(os.path.join(project_root, "data"), exist_ok=True)
            with open(os.path.join(project_root, ".jarvis-project.json"), "w", encoding="utf-8") as f1:
                f1.write("{}\n")
            with open(os.path.join(project_root, "jarvis.project.yaml"), "w", encoding="utf-8") as f1:
                f1.write("project: {}\n")
            with open(os.path.join(project_root, "README.md"), "w", encoding="utf-8") as f1:
                f1.write("# Demo\n")
            with open(os.path.join(project_root, "bin", "task.yaml"), "w", encoding="utf-8") as f1:
                f1.write("Scan: {name: demo, save_dir: '&J/outputs'}\n")

            old_cwd = os.getcwd()
            out = io.StringIO()
            try:
                os.chdir(project_root)
                with contextlib.redirect_stdout(out):
                    rc = main(["Jarvis", "--packproject"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(rc, 0, msg=out.getvalue())
            archives = glob.glob(os.path.join(tmpdir, "DemoProject_repro_*.tar.gz"))
            self.assertEqual(len(archives), 1, msg=out.getvalue())
            with tarfile.open(archives[0], "r:gz") as tf:
                names = set(tf.getnames())
            self.assertIn("DemoProject/.jarvis-project.json", names)
            self.assertIn("DemoProject/.jarvis-pack/manifest.json", names)

    def test_main_packproject_non_project_path_fails(self):
        from jarvishep.client import main

        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            out = io.StringIO()
            try:
                os.chdir(tmpdir)
                with contextlib.redirect_stdout(out):
                    with self.assertRaises(SystemExit) as cm:
                        main(["Jarvis", "--packproject"])
            finally:
                os.chdir(old_cwd)

            self.assertEqual(cm.exception.code, 2)
            self.assertIn("Jarvis project not found", out.getvalue())

    def test_python_module_mkproject_success(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            env = dict(os.environ)
            pythonpath = env.get("PYTHONPATH", "")
            env["PYTHONPATH"] = PROJECT_ROOT if not pythonpath else f"{PROJECT_ROOT}{os.pathsep}{pythonpath}"
            proc = subprocess.run(
                [sys.executable, "-m", "jarvishep", "--mkproject", "SubprocProject"],
                cwd=tmpdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stdout)
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "SubprocProject", "bin")))

    def test_main_version_fast_path_long_flag(self):
        from jarvishep.client import main

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["Jarvis", "--version"])

        self.assertEqual(rc, 0)
        text = out.getvalue()
        self.assertIn("=== Jarvis-HEP ===", text)
        self.assertIn("Author:", text)
        self.assertIn("Version:", text)
        self.assertIn("Resources:", text)
        self.assertIn("\tOnline docs:\t", text)
        self.assertIn("\tHomepage:\t", text)
        self.assertIn("References:", text)
        self.assertIn("\tJarvis-HEP:", text)
        self.assertIn("\tBuilt-in Scanners:", text)
        self.assertNotIn("\tBuilt-in Tools:", text)
        self.assertIn("[1] Core Jarvis-HEP framework paper", text)
        self.assertIn("Title:\t", text)
        self.assertIn("arXiv:\t", text)
        self.assertIn("DOI:\t", text)
        self.assertNotIn("edit:", text)
        self.assertNotIn("N/A", text)
        self.assertNotIn("TODO", text)

    def test_main_version_fast_path_short_flag(self):
        from jarvishep.client import main

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["Jarvis", "-v"])

        self.assertEqual(rc, 0)
        text = out.getvalue()
        self.assertIn("=== Jarvis-HEP ===", text)
        self.assertIn("Author:", text)
        self.assertIn("Version:", text)
        self.assertIn("Resources:", text)
        self.assertIn("\tOnline docs:\t", text)
        self.assertIn("\tHomepage:\t", text)
        self.assertIn("References:", text)
        self.assertIn("\tJarvis-HEP:", text)
        self.assertIn("\tBuilt-in Scanners:", text)
        self.assertNotIn("\tBuilt-in Tools:", text)
        self.assertIn("[1] Core Jarvis-HEP framework paper", text)
        self.assertIn("Title:\t", text)
        self.assertIn("arXiv:\t", text)
        self.assertIn("DOI:\t", text)
        self.assertNotIn("edit:", text)
        self.assertNotIn("N/A", text)
        self.assertNotIn("TODO", text)


if __name__ == "__main__":
    unittest.main()
