#!/usr/bin/env python3
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


class CliEntrypointTests(unittest.TestCase):
    def test_import_jarvishep_client_main(self):
        from jarvishep.client import main

        self.assertTrue(callable(main))

    def test_pyproject_jarvis_entrypoint(self):
        pyproject = os.path.join(PROJECT_ROOT, "pyproject.toml")
        with open(pyproject, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn('[project.scripts]', content)
        self.assertIn('Jarvis = "jarvishep.client:main"', content)

    def test_python_module_help(self):
        proc = subprocess.run(
            [sys.executable, "-m", "jarvishep", "--help"],
            cwd=PROJECT_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        self.assertEqual(proc.returncode, 0, msg=proc.stdout)
        self.assertIn("Jarvis Program Help Center", proc.stdout)

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


if __name__ == "__main__":
    unittest.main()
