#!/usr/bin/env python3
"""Jarvis CLI parse + dispatch routing tests (D11.1 / D11.2)."""

from __future__ import annotations

import contextlib
import io
import os
import subprocess
import sys
import unittest
from unittest import mock

from jarvishep2.client import (
    _log_outcome,
    build_parser,
    build_project_help_parser,
    dispatch,
    main,
    normalize_argv,
    resolve_intent,
)
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED, EXIT_USAGE, RunOutcome


TESTS_ROOT = os.path.dirname(__file__)
CHECK_MODULES_YAML = os.path.join(TESTS_ROOT, "parity_project", "check_modules.yaml")


class CliNormalizeArgvTests(unittest.TestCase):
    def test_normalize_argv_rewrites_only_legacy_shortcuts(self) -> None:
        cases = (
            # A bare path is not rewritten into ``run``.
            ([CHECK_MODULES_YAML, "--resume"], [CHECK_MODULES_YAML, "--resume"]),
            (
                [CHECK_MODULES_YAML, "--check-modules"],
                [CHECK_MODULES_YAML, "--check-modules"],
            ),
            (["scene.yaml", "--plot"], ["plot", "scene.yaml"]),
            (["task.yaml", "--convert"], ["convert", "task.yaml"]),
            (["task.yaml", "--validate", "--strict"], ["validate", "task.yaml", "--strict"]),
            (["monitor"], ["monitor"]),
            # ``Jarvis --monitor`` keeps the top-level flag (resolve_intent → monitor).
            (["--monitor"], ["--monitor"]),
        )
        for argv, expected in cases:
            with self.subTest(argv=argv):
                self.assertEqual(normalize_argv(argv), expected)


class CliHelpFallbackTests(unittest.TestCase):
    def _root_help(self) -> str:
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            build_parser().print_help(stdout)
        return stdout.getvalue()

    def test_unknown_help_prints_root_help(self) -> None:
        expected = self._root_help()
        cases = (
            ["opera", "-h"],
            ["opera", "--help"],
            ["opera", "--validate", "-h"],
            ["task.yaml", "--plot", "--help"],
            ["--validate", "-h"],
        )
        for argv in cases:
            with self.subTest(argv=argv):
                stdout = io.StringIO()
                stderr = io.StringIO()
                with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                    code = main(argv)
                self.assertEqual(code, EXIT_OK)
                self.assertEqual(stdout.getvalue(), expected)
                self.assertNotIn("Usage: Jarvis run [OPTIONS]", stdout.getvalue())
                self.assertNotIn("Usage: Jarvis validate [OPTIONS]", stdout.getvalue())

    def test_legal_command_help_stays_on_that_command(self) -> None:
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as raised:
                main(["run", "-h"])
        self.assertEqual(raised.exception.code, 0)
        self.assertIn("Usage: Jarvis run [OPTIONS] <task_yaml>", stdout.getvalue())

    def test_bare_task_path_is_not_run(self) -> None:
        for argv in (["task.yaml"], ["task.yaml", "--resume"], ["opera"]):
            with self.subTest(argv=argv):
                stderr = io.StringIO()
                with contextlib.redirect_stderr(stderr), mock.patch(
                    "jarvishep2.core.Jarvis2Core"
                ) as core:
                    code = main(argv)
                self.assertEqual(code, EXIT_USAGE)
                self.assertIn("Unknown command", stderr.getvalue())
                self.assertIn("Jarvis -h", stderr.getvalue())
                core.assert_not_called()


class CliParseTests(unittest.TestCase):
    def test_run_debug_shortcut_enables_both_debug_levels(self) -> None:
        args = build_parser().parse_args(["run", "task.yaml", "-d"])
        self.assertTrue(args.debug)

        core = mock.Mock()
        core.run.return_value = RunOutcome(submitted=1, completed=1)
        with mock.patch("jarvishep2.core.Jarvis2Core", return_value=core):
            self.assertEqual(dispatch(args), EXIT_OK)
        core.set_logging_options.assert_called_once_with(
            console_level="DEBUG", silence=False
        )

    def test_removed_options_are_rejected(self) -> None:
        for argv in (
            ["run", "task.yaml", "--log-level", "DEBUG"],
            ["convert", "task.yaml", "--force"],
        ):
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    build_parser().parse_args(argv)

    def test_run_help_explains_console_level_values(self) -> None:
        help_text = build_parser()._subparsers._group_actions[0].choices["run"].format_help()
        self.assertIn("Usage: Jarvis run [OPTIONS] <task_yaml>", help_text)
        self.assertIn("--console-level <str>", help_text)
        for level in ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"):
            self.assertIn(level, help_text)
        self.assertIn("default:", help_text)
        self.assertIn("WARNING).", help_text)

    def test_parse_sets_command_arguments_and_defaults(self) -> None:
        cases = (
            (
                ["run", CHECK_MODULES_YAML, "--resume"],
                {"command": "run", "task_yaml": CHECK_MODULES_YAML, "resume": True},
            ),
            (["run", "task.yaml"], {"console_level": "WARNING"}),
            (["check", CHECK_MODULES_YAML], {"command": "check", "console_level": "DEBUG"}),
            (["run", "task.yaml", "--skip-library-installation"], {"skip_library_installation": True}),
            (
                ["check", "task.yaml", "--skip-library-installation"],
                {"skip_library_installation": True},
            ),
            (["monitor", "R1"], {"command": "monitor", "scan_ref": "R1"}),
            (
                ["gen-plot-yaml", CHECK_MODULES_YAML],
                {"command": "gen-plot-yaml", "task_yaml": CHECK_MODULES_YAML},
            ),
            (["convert", "task.yaml"], {"command": "convert", "task_yaml": "task.yaml"}),
            (["--refs"], {"refs": True}),
            (["-v"], {"version": True}),
            (["-V"], {"version": True}),
            (["--version"], {"version": True}),
        )
        for argv, expected in cases:
            with self.subTest(argv=argv):
                args = build_parser().parse_args(argv)
                for name, value in expected.items():
                    self.assertEqual(getattr(args, name), value, name)

    def test_help_uses_jarvis_box_sections(self) -> None:
        help_text = build_parser().format_help()
        self.assertIn("Usage: Jarvis [OPTIONS] COMMAND [ARGS]...", help_text)
        self.assertIn("Command help: Jarvis COMMAND -h", help_text)
        self.assertIn("See each command's arguments and options", help_text)
        self.assertIn("╭─ YAML authoring", help_text)
        self.assertIn("╭─ Scan workflow", help_text)
        self.assertIn("╭─ Data export", help_text)
        self.assertIn("╭─ Plots", help_text)
        self.assertIn("╭─ Runtime control", help_text)
        self.assertIn("╭─ Legacy options", help_text)
        self.assertGreater(help_text.index("╭─ Legacy options"), help_text.index("╭─ Extensions & projects"))
        self.assertIn("╭─ Options", help_text)
        self.assertIn("Run a distributed scan task YAML", help_text)
        self.assertIn("coding agents: use", help_text)
        self.assertIn("--json for structured keys", help_text)
        self.assertNotIn("│ man", help_text[help_text.index("╭─ Scan workflow"):help_text.index("╭─ Data export")])
        self.assertIn("--help", help_text)
        self.assertIn("-h", help_text)
        self.assertIn("╰", help_text)

    def test_help_uses_a_shared_fixed_description_column(self) -> None:
        help_text = build_parser().format_help()

        def column_of(fragment: str) -> int:
            line = next(line for line in help_text.splitlines() if fragment in line)
            return line.index(fragment)

        self.assertEqual(
            column_of("Print Jarvis-HEP logo"),
            column_of("Run a distributed scan task YAML"),
        )
        self.assertEqual(
            column_of("Run a distributed scan task YAML"),
            column_of("(legacy) List running scans"),
        )

    def test_help_orders_commands_by_workflow(self) -> None:
        help_text = build_parser().format_help()
        scan_steps = (
            "Validate task YAML",
            "Run fixed-point calculator smoke",
            "Run a distributed scan task YAML",
            "Refresh CSV snapshots",
        )
        runtime_steps = (
            "Open the live TUI, or print one snapshot",
            "List process groups, or show one scan / ZP",
            "List process groups, or terminate one",
        )
        self.assertEqual(
            sorted(scan_steps, key=help_text.index), list(scan_steps)
        )
        self.assertEqual(
            sorted(runtime_steps, key=help_text.index), list(runtime_steps)
        )

    def test_project_help_uses_the_shared_panel_layout(self) -> None:
        help_text = build_project_help_parser().format_help()
        self.assertIn("Create, package, browse, fetch", help_text)
        self.assertIn("╭─ Options", help_text)
        self.assertIn("╰", help_text)


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


class CliReferencesTests(unittest.TestCase):
    def test_old_refs_subcommand_is_rejected_with_migration_hint(self) -> None:
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            code = main(["refs"])
        self.assertEqual(code, EXIT_USAGE)
        self.assertIn("Jarvis --refs", stderr.getvalue())

    def test_refs_prints_framework_and_supported_sampler_citations(self) -> None:
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            code = main(["--refs"])
        self.assertEqual(code, EXIT_OK)
        rendered = stdout.getvalue()
        self.assertIn("Jarvis-HEP V2 references", rendered)
        self.assertIn("Bridson / AdaptiveBridson", rendered)
        self.assertIn("MultiNest", rendered)
        self.assertIn("Dynesty", rendered)
        self.assertIn("MCMC family", rendered)


class CliDispatchTests(unittest.TestCase):
    def test_run_outcome_is_written_to_process_logger(self) -> None:
        outcome = RunOutcome(
            submitted=10,
            completed=10,
            failed=0,
            archived=10,
            extras={"ncall": 42},
        )
        logger = mock.Mock()
        with mock.patch("jarvishep2.cli.run.get_jarvis_logger", return_value=logger):
            _log_outcome(outcome)

        logger.warning.assert_called_once_with(
            "RunOutcome status=success submitted=10 completed=10 "
            "failed=0 archived=10 ncall=42 exit=0"
        )
        logger.error.assert_not_called()

    def test_dispatch_routes_gen_plot_yaml(self) -> None:
        args = build_parser().parse_args(["gen-plot-yaml", CHECK_MODULES_YAML])
        with mock.patch("jarvishep2.client.dispatch_gen_plot_yaml", return_value=0) as scene:
            self.assertEqual(dispatch(args), 0)
        scene.assert_called_once_with(CHECK_MODULES_YAML)

    def test_dispatch_routes_check_modules_to_core(self) -> None:
        args = build_parser().parse_args(
            ["check", CHECK_MODULES_YAML]
        )
        core = mock.Mock()
        core.check_modules.return_value = RunOutcome(
            submitted=10, completed=10, failed=0, archived=10
        )
        with mock.patch("jarvishep2.core.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, EXIT_OK)
        core.load_task_yaml.assert_called_once_with(
            CHECK_MODULES_YAML,
            validate=True,
            strict=False,
            check_modules=True,
        )
        core.set_logging_options.assert_called_once_with(
            console_level="DEBUG", silence=False
        )
        core.check_modules.assert_called_once_with(timeout=None)
        core.run.assert_not_called()

    def test_dispatch_routes_resume_run_to_core(self) -> None:
        args = build_parser().parse_args(
            ["run", CHECK_MODULES_YAML, "--resume"]
        )
        core = mock.Mock()
        core.run.return_value = RunOutcome(
            submitted=10, completed=10, failed=0, archived=10
        )
        with mock.patch("jarvishep2.core.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, EXIT_OK)
        core.run.assert_called_once_with(resume=True)

    def test_dispatch_exit_codes(self) -> None:
        # (argv, what Core.run returns or raises, expected exit). None: Core
        # is not mocked; the command fails before or while loading the card.
        cases = (
            # Archived rows exist but every sample failed: must not exit 0.
            (
                ["run", CHECK_MODULES_YAML],
                RunOutcome(submitted=5, completed=0, failed=5, archived=5),
                EXIT_RUN_FAILED,
            ),
            (
                ["run", CHECK_MODULES_YAML],
                RunOutcome(submitted=5, completed=3, failed=2, archived=5),
                EXIT_RUN_FAILED,
            ),
            (["run", CHECK_MODULES_YAML], NotImplementedError("unsupported sampler"), EXIT_USAGE),
            # No subcommand and no path.
            ([], None, EXIT_USAGE),
            (["check", "/no/such/task.yaml"], None, EXIT_RUN_FAILED),
            (["--monitor", "--plot"], None, EXIT_USAGE),
        )
        for argv, result, expected in cases:
            with self.subTest(argv=argv, result=result):
                args = build_parser().parse_args(argv)
                if result is None:
                    self.assertEqual(dispatch(args), expected)
                    continue
                core = mock.Mock()
                if isinstance(result, Exception):
                    core.run.side_effect = result
                else:
                    core.run.return_value = result
                with mock.patch("jarvishep2.core.Jarvis2Core", return_value=core):
                    self.assertEqual(dispatch(args), expected)

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

    def test_operas_passthrough_forwards_the_complete_jopera_argv(self) -> None:
        from jarvishep2.client import main

        with mock.patch("jarvis_operas.cli.main", return_value=7) as operas_main:
            self.assertEqual(main(["operas", "list", "--namespace", "math"]), 7)
        operas_main.assert_called_once_with(["list", "--namespace", "math"])


class CliImportIsolationTests(unittest.TestCase):
    """D25.9: help/project chrome must not pull the runtime graph."""

    _BANNED = (
        "jarvishep2.worker",
        "jarvishep2.runtime.worker",
        "jarvishep2.core",
        "jarvishep2.runtime.core",
        "jarvishep2.factory",
        "jarvishep2.runtime.factory",
        "jarvishep2.redis_queue",
        "jarvishep2.queue.redis_queue",
    )

    def _assert_argv_stays_light(self, argv: list[str]) -> None:
        banned = ", ".join(repr(name) for name in self._BANNED)
        argv_literal = repr(list(argv))
        script = f"""
import sys
from jarvishep2.client import main
try:
    code = main({argv_literal})
except SystemExit as exc:
    code = 0 if exc.code is None else int(exc.code)
if code != 0:
    raise SystemExit("unexpected exit " + str(code))
loaded = [name for name in ({banned},) if name in sys.modules]
if loaded:
    raise SystemExit("unexpected modules: " + ",".join(loaded))
"""
        completed = subprocess.run(
            [sys.executable, "-c", script],
            cwd=os.path.dirname(TESTS_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr or completed.stdout)

    def test_help_does_not_import_the_runtime(self) -> None:
        for argv in (["project", "--help"], ["-h"]):
            with self.subTest(argv=argv):
                self._assert_argv_stays_light(argv)


if __name__ == "__main__":
    unittest.main()
