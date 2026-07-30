#!/usr/bin/env python3
"""Jarvis2 CLI parse + dispatch routing tests (D11.1 / D11.2)."""

from __future__ import annotations

import os
import sys
import unittest
from unittest import mock

from jarvishep2.client import (
    build_parser,
    dispatch,
    main,
    normalize_argv,
    resolve_intent,
)
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED, EXIT_USAGE, RunOutcome


TESTS_ROOT = os.path.dirname(__file__)
CHECK_MODULES_YAML = os.path.join(TESTS_ROOT, "parity_project", "check_modules.yaml")


class CliNormalizeArgvTests(unittest.TestCase):
    def test_bare_yaml_becomes_run(self) -> None:
        self.assertEqual(
            normalize_argv([CHECK_MODULES_YAML, "--resume"]),
            ["run", CHECK_MODULES_YAML, "--resume"],
        )

    def test_legacy_check_modules_flag(self) -> None:
        self.assertEqual(
            normalize_argv([CHECK_MODULES_YAML, "--check-modules"]),
            ["check", CHECK_MODULES_YAML],
        )

    def test_legacy_plot_flag(self) -> None:
        self.assertEqual(
            normalize_argv(["scene.yaml", "--plot"]),
            ["plot", "scene.yaml"],
        )

    def test_explicit_subcommand_unchanged(self) -> None:
        self.assertEqual(
            normalize_argv(["monitor"]),
            ["monitor"],
        )

    def test_legacy_top_level_monitor_flag(self) -> None:
        """``Jarvis2 --monitor`` keeps the top-level flag (resolve_intent → monitor)."""
        self.assertEqual(normalize_argv(["--monitor"]), ["--monitor"])


class CliParseTests(unittest.TestCase):
    def test_help_uses_jarvis_box_sections(self) -> None:
        help_text = build_parser().format_help()
        self.assertIn("Usage: Jarvis2 [OPTIONS] COMMAND [ARGS]...", help_text)
        self.assertIn("╭─ Commands", help_text)
        self.assertIn("╭─ Options", help_text)
        self.assertIn("Run a distributed scan task YAML", help_text)
        self.assertIn("--help", help_text)
        self.assertIn("-h", help_text)
        self.assertIn("╰", help_text)

    def test_parse_run_subcommand(self) -> None:
        args = build_parser().parse_args(["run", CHECK_MODULES_YAML, "--resume"])
        self.assertEqual(args.command, "run")
        self.assertEqual(args.task_yaml, CHECK_MODULES_YAML)
        self.assertTrue(args.resume)

    def test_parse_check_subcommand(self) -> None:
        args = build_parser().parse_args(["check", CHECK_MODULES_YAML])
        self.assertEqual(args.command, "check")

    def test_parse_monitor_subcommand(self) -> None:
        args = build_parser().parse_args(["monitor"])
        self.assertEqual(args.command, "monitor")

    def test_parse_plot_scene_subcommand(self) -> None:
        args = build_parser().parse_args(["plot-scene", CHECK_MODULES_YAML])
        self.assertEqual(args.command, "plot-scene")
        self.assertEqual(args.task_yaml, CHECK_MODULES_YAML)


class CliResolveIntentTests(unittest.TestCase):
    def test_conflict_monitor_and_plot(self) -> None:
        args = build_parser().parse_args(["--monitor", "--plot"])
        with self.assertRaises(ValueError):
            resolve_intent(args)

    def test_pid_hard_fails(self) -> None:
        args = build_parser().parse_args(["run", CHECK_MODULES_YAML])
        args.pid = 123  # retired flag still rejected when present
        with self.assertRaisesRegex(ValueError, "--pid"):
            resolve_intent(args)


class CliVersionTests(unittest.TestCase):
    def test_parse_version_flags(self) -> None:
        for flag in ("-v", "--version"):
            args = build_parser().parse_args([flag])
            self.assertTrue(args.version)

    def test_dispatch_version_prints_logo(self) -> None:
        from jarvishep2.client import dispatch_version

        with mock.patch("builtins.print") as printer:
            code = dispatch_version()
        self.assertEqual(code, EXIT_OK)
        printer.assert_called_once()
        banner = str(printer.call_args[0][0])
        self.assertIn("Version:", banner)
        self.assertIn("Author:", banner)

    def test_main_version_short_circuit(self) -> None:
        with mock.patch("jarvishep2.client.dispatch_version", return_value=0) as ver:
            code = main(["--version"])
        self.assertEqual(code, 0)
        ver.assert_called_once()


class CliDispatchTests(unittest.TestCase):
    def test_dispatch_routes_plot_scene(self) -> None:
        args = build_parser().parse_args(["plot-scene", CHECK_MODULES_YAML])
        with mock.patch("jarvishep2.client.dispatch_plot_scene", return_value=0) as scene:
            self.assertEqual(dispatch(args), 0)
        scene.assert_called_once_with(CHECK_MODULES_YAML)

    def test_dispatch_routes_check_modules_to_core(self) -> None:
        args = build_parser().parse_args(
            normalize_argv([CHECK_MODULES_YAML, "--check-modules"])
        )
        core = mock.Mock()
        core.check_modules.return_value = RunOutcome(
            submitted=10, completed=10, failed=0, archived=10
        )
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, EXIT_OK)
        core.load_task_yaml.assert_called_once_with(
            CHECK_MODULES_YAML,
            validate=True,
            strict=False,
            check_modules=True,
        )
        core.check_modules.assert_called_once_with()
        core.run.assert_not_called()

    def test_dispatch_routes_resume_run_to_core(self) -> None:
        args = build_parser().parse_args(
            normalize_argv([CHECK_MODULES_YAML, "--resume"])
        )
        core = mock.Mock()
        core.run.return_value = RunOutcome(
            submitted=10, completed=10, failed=0, archived=10
        )
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, EXIT_OK)
        core.run.assert_called_once_with(resume=True)

    def test_dispatch_all_failed_exits_nonzero(self) -> None:
        args = build_parser().parse_args(["run", CHECK_MODULES_YAML])
        core = mock.Mock()
        # Archived rows exist but every sample failed — must not exit 0.
        core.run.return_value = RunOutcome(
            submitted=5, completed=0, failed=5, archived=5
        )
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, EXIT_RUN_FAILED)

    def test_dispatch_partial_failure_exits_nonzero(self) -> None:
        args = build_parser().parse_args(["run", CHECK_MODULES_YAML])
        core = mock.Mock()
        core.run.return_value = RunOutcome(
            submitted=5, completed=3, failed=2, archived=5
        )
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, EXIT_RUN_FAILED)

    def test_dispatch_missing_task_yaml_returns_usage_exit(self) -> None:
        # No subcommand / no path → usage.
        args = build_parser().parse_args([])
        code = dispatch(args)
        self.assertEqual(code, EXIT_USAGE)

    def test_dispatch_missing_file_is_clear_error(self) -> None:
        args = build_parser().parse_args(["check", "/no/such/task.yaml"])
        code = dispatch(args)
        self.assertEqual(code, EXIT_RUN_FAILED)

    def test_dispatch_not_implemented_returns_exit_2(self) -> None:
        args = build_parser().parse_args(["run", CHECK_MODULES_YAML])
        core = mock.Mock()
        core.run.side_effect = NotImplementedError("unsupported sampler")
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, EXIT_USAGE)

    def test_dispatch_mode_conflict_exits_usage(self) -> None:
        args = build_parser().parse_args(["--monitor", "--plot"])
        code = dispatch(args)
        self.assertEqual(code, EXIT_USAGE)

    def test_main_does_not_import_v1_client(self) -> None:
        self.assertNotIn("jarvishep.client", sys.modules)
        self.assertNotIn(
            "jarvishep",
            {name.split(".")[0] for name in sys.modules if name.startswith("jarvishep")},
        )

    def test_portal_passthrough_lists_v2_formats(self) -> None:
        from jarvishep2.client import dispatch_portal

        # Same surface as jportal man; V2 registry includes SLHA.
        code = dispatch_portal(["man"])
        self.assertEqual(code, EXIT_OK)

    def test_portal_formats_alias(self) -> None:
        from jarvishep2.client import dispatch_portal

        self.assertEqual(dispatch_portal(["formats"]), EXIT_OK)


if __name__ == "__main__":
    unittest.main()
