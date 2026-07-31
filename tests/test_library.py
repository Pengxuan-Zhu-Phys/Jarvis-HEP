#!/usr/bin/env python3
"""LibDeps / LibraryManager tests (WP-D2.3 design-doc parity)."""

from __future__ import annotations

import json
import os
import tempfile
import threading
import unittest
from unittest import mock

from jarvishep2.command_parser import CommandParser
from jarvishep2.library import (
    LIBRARY_INSTALL_CONTROL_BASENAME,
    LIBRARY_INSTALL_STAMP_BASENAME,
    LibraryInstallError,
    LibraryInstaller,
    LibraryManager,
)
from jarvishep2.core import Jarvis2Core

TESTS_ROOT = os.path.dirname(__file__)
SAFE_SOURCE = os.path.join(TESTS_ROOT, "fixtures", "safe_calc", "source")
SHADOW_SOURCE = os.path.join(TESTS_ROOT, "fixtures", "shadow_calc", "source")


class LibraryManagerTests(unittest.TestCase):
    def test_requires_shadow_flag(self) -> None:
        self.assertTrue(LibraryManager.requires_shadow(True))
        self.assertFalse(LibraryManager.requires_shadow(False))

    def test_link_into_sample_creates_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, "SAMPLE", "000001")
            link_path = LibraryManager().link_into_sample(SAFE_SOURCE, sample_dir, "SafeCalc")
            self.assertTrue(os.path.islink(link_path))
            self.assertEqual(os.path.realpath(link_path), os.path.realpath(SAFE_SOURCE))

    def test_safe_vs_shadow_isolation_policy(self) -> None:
        manager = LibraryManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, "SAMPLE", "safe-sample")
            link_path = manager.link_into_sample(SAFE_SOURCE, sample_dir, "SafeCalc")
            self.assertTrue(os.path.islink(link_path))
            self.assertFalse(manager.requires_shadow(False))
            self.assertTrue(manager.requires_shadow(True))

    def test_concurrent_link_into_sample_is_race_safe(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            errors: list[BaseException] = []
            lock = threading.Lock()

            def _link(sample_index: int) -> None:
                try:
                    sample_dir = os.path.join(tmpdir, "SAMPLE", f"{sample_index:06d}")
                    link_path = LibraryManager().link_into_sample(
                        SAFE_SOURCE,
                        sample_dir,
                        "SafeCalc",
                    )
                    if not os.path.islink(link_path):
                        raise AssertionError(f"expected symlink at {link_path}")
                except BaseException as exc:
                    with lock:
                        errors.append(exc)

            threads = [threading.Thread(target=_link, args=(index,)) for index in range(8)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=5.0)
            self.assertEqual(errors, [])

    def test_missing_source_raises_clear_boot_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, "SAMPLE", "missing")
            missing = os.path.join(tmpdir, "does-not-exist", "tool")
            with self.assertRaises(FileNotFoundError) as ctx:
                LibraryManager().link_into_sample(missing, sample_dir, "MissingTool")
            self.assertIn("LibDeps source does not exist", str(ctx.exception))


class LibraryInstallerTests(unittest.TestCase):
    def _installer(self, root: str, *, modules: list[dict], skip: bool = False) -> LibraryInstaller:
        libdeps_path = os.path.join(root, "libraries")
        config = {
            "task_result_dir": root,
            "LibDeps": {
                "path": libdeps_path,
                "make_paraller": 2,
                "Modules": modules,
            },
        }
        parser = CommandParser.from_config(config, project_root=root)
        return LibraryInstaller(
            config,
            parser=parser,
            logs_dir=os.path.join(root, "logs", "library-test"),
            skip_installation=skip,
        )

    def test_installs_dependency_layers_and_reuses_matching_stamps(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            core_path = os.path.join(root, "libraries", "Core")
            addon_path = os.path.join(root, "libraries", "Addon")
            modules = [
                {
                    "name": "Core",
                    "installation": {
                        "path": core_path,
                        "commands": [
                            "mkdir -p ${path}",
                            "cd ${path}",
                            "printf core > installed.txt",
                        ],
                    },
                },
                {
                    "name": "Addon",
                    "required_modules": ["Core"],
                    "installation": {
                        "path": addon_path,
                        "commands": [
                            "test -f ${LibDeps:Core}/installed.txt",
                            "mkdir -p ${path}",
                            "cd ${path}",
                            "printf addon > installed.txt",
                        ],
                    },
                },
            ]
            installer = self._installer(root, modules=modules)
            first = installer.prepare()
            self.assertEqual({name: item["status"] for name, item in first.items()}, {
                "Core": "installed", "Addon": "installed"
            })
            self.assertTrue(os.path.isfile(os.path.join(core_path, LIBRARY_INSTALL_STAMP_BASENAME)))
            self.assertTrue(os.path.isfile(os.path.join(addon_path, LIBRARY_INSTALL_STAMP_BASENAME)))
            self.assertTrue(os.path.isfile(os.path.join(root, "libraries", LIBRARY_INSTALL_CONTROL_BASENAME)))
            self.assertTrue(os.path.isfile(os.path.join(root, "logs", "library-test", "library-Core.log")))

            second = self._installer(root, modules=modules).prepare()
            self.assertEqual({name: item["status"] for name, item in second.items()}, {
                "Core": "reused", "Addon": "reused"
            })

    def test_control_reinstall_requests_next_build_without_interaction(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "libraries", "Demo")
            modules = [{"name": "Demo", "installation": {"path": path, "commands": ["true"]}}]
            installer = self._installer(root, modules=modules)
            installer.prepare()
            control_path = os.path.join(root, "libraries", LIBRARY_INSTALL_CONTROL_BASENAME)
            with open(control_path, encoding="utf-8") as handle:
                control = json.load(handle)
            control["reinstall"] = True
            with open(control_path, "w", encoding="utf-8") as handle:
                json.dump(control, handle)

            rebuilt = self._installer(root, modules=modules).prepare()
            self.assertEqual(rebuilt["Demo"]["status"], "installed")
            with open(control_path, encoding="utf-8") as handle:
                refreshed = json.load(handle)
            self.assertFalse(refreshed["reinstall"])
            self.assertEqual(refreshed["reinstall_epoch"], 1)

    def test_v1_command_can_create_its_own_installation_path(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "libraries", "Pythia8")
            result = self._installer(
                root,
                modules=[
                    {
                        "name": "Pythia8",
                        "installation": {
                            "path": path,
                            "commands": ["cd ${LibDeps:path}", "mkdir Pythia8"],
                        },
                    }
                ],
            ).prepare()
            self.assertEqual(result["Pythia8"]["status"], "installed")
            self.assertTrue(os.path.isdir(path))

    def test_failure_is_clear_and_leaves_no_success_stamp(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "libraries", "Broken")
            installer = self._installer(
                root,
                modules=[{"name": "Broken", "installation": {"path": path, "commands": ["false"]}}],
            )
            with self.assertRaisesRegex(LibraryInstallError, r"library-Broken\.log"):
                installer.prepare()
            self.assertFalse(os.path.exists(os.path.join(path, LIBRARY_INSTALL_STAMP_BASENAME)))

    def test_skip_installation_checks_the_declared_library_path(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "libraries", "External")
            modules = [{"name": "External", "installation": {"path": path, "commands": ["false"]}}]
            with self.assertRaisesRegex(LibraryInstallError, r"skip-library-installation"):
                self._installer(root, modules=modules, skip=True).prepare()
            os.makedirs(path)
            result = self._installer(root, modules=modules, skip=True).prepare()
            self.assertEqual(result["External"]["status"], "skipped")

    def test_invalid_dependency_graph_is_rejected_before_commands_run(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            modules = [
                {"name": "A", "required_modules": ["B"], "installation": {"commands": ["false"]}},
                {"name": "B", "required_modules": ["A"], "installation": {"commands": ["false"]}},
            ]
            with self.assertRaisesRegex(LibraryInstallError, "dependency cycle"):
                self._installer(root, modules=modules).prepare()

    def test_prebuild_parser_defers_registered_executable_until_library_exists(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            library_path = os.path.join(root, "libraries", "Tool")
            executable = os.path.join(library_path, "tool")
            config = {
                "task_result_dir": root,
                "LibDeps": {
                    "path": os.path.join(root, "libraries"),
                    "Modules": [
                        {
                            "name": "Tool",
                            "installation": {
                                "path": library_path,
                                "commands": [
                                    "mkdir -p ${path}",
                                    "printf '#!/bin/sh\\nexit 0\\n' > ${path}/tool",
                                ],
                            },
                        }
                    ],
                    "registered_executables": [
                        {
                            "name": "tool",
                            "source": "${LibDeps:Tool}/tool",
                            "resolution": "direct_path",
                        }
                    ],
                },
            }
            with self.assertRaises(FileNotFoundError):
                CommandParser.from_config(config, project_root=root)
            prebuild = CommandParser.from_config(
                config, project_root=root, register_executables=False
            )
            LibraryInstaller(
                config,
                parser=prebuild,
                logs_dir=os.path.join(root, "logs", "library-test"),
            ).prepare()
            parser = CommandParser.from_config(config, project_root=root)
            self.assertEqual(parser.registered["tool"].path, executable)

    def test_build_failure_stops_core_before_redis_initialization(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            core = Jarvis2Core(
                {
                    "task_root": root,
                    "scan_name": "broken-library",
                    "LibDeps": {
                        "path": os.path.join(root, "libraries"),
                        "Modules": [
                            {"name": "Broken", "installation": {"commands": ["false"]}}
                        ],
                    },
                }
            )
            core.info.update({"task_root": root, "scan_name": "broken-library"})
            core.init_logger = mock.Mock()
            core.check_environment_requirements = mock.Mock()
            core.init_redis = mock.Mock()
            with mock.patch.object(core, "is_redis_runtime", return_value=True), \
                 mock.patch("jarvishep2.process_cleanup.ensure_scan_name_available"):
                with self.assertRaises(LibraryInstallError):
                    core.bootstrap_distributed_runtime()
            core.init_redis.assert_not_called()


if __name__ == "__main__":
    unittest.main()
