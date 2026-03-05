#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest


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
        from jarvishep.project_scaffold import PROJECT_SUBDIRS

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

    def test_main_mkproject_conflict_flag(self):
        from jarvishep.client import main

        rc = main(["Jarvis", "--mkproject", "DemoProject", "--plot"])
        self.assertEqual(rc, 2)

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
