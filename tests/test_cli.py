#!/usr/bin/env python3
"""Jarvis2 CLI parse + dispatch routing tests."""

from __future__ import annotations

import os
import sys
import unittest
from unittest import mock

from jarvishep2.client import build_parser, dispatch, main


TESTS_ROOT = os.path.dirname(__file__)
CHECK_MODULES_YAML = os.path.join(TESTS_ROOT, "parity_project", "check_modules.yaml")


class CliParseTests(unittest.TestCase):
    def test_parse_run_with_resume(self) -> None:
        args = build_parser().parse_args([CHECK_MODULES_YAML, "--resume"])
        self.assertEqual(args.task_yaml, CHECK_MODULES_YAML)
        self.assertTrue(args.resume)
        self.assertFalse(args.check_modules)
        self.assertFalse(args.monitor)

    def test_parse_check_modules_flag(self) -> None:
        args = build_parser().parse_args([CHECK_MODULES_YAML, "--check-modules"])
        self.assertTrue(args.check_modules)

    def test_parse_monitor_with_redis_overrides(self) -> None:
        args = build_parser().parse_args(
            ["task.yaml", "--monitor", "--redis-host", "10.0.0.2", "--redis-port", "6380", "--redis-db", "3"]
        )
        self.assertTrue(args.monitor)
        self.assertEqual(args.redis_host, "10.0.0.2")
        self.assertEqual(args.redis_port, 6380)
        self.assertEqual(args.redis_db, 3)


class CliDispatchTests(unittest.TestCase):
    def test_dispatch_routes_check_modules_to_core(self) -> None:
        args = build_parser().parse_args([CHECK_MODULES_YAML, "--check-modules"])
        core = mock.Mock()
        core.check_modules.return_value = 10
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, 0)
        core.load_task_yaml.assert_called_once_with(CHECK_MODULES_YAML)
        core.check_modules.assert_called_once_with()
        core.run.assert_not_called()

    def test_dispatch_routes_resume_run_to_core(self) -> None:
        args = build_parser().parse_args([CHECK_MODULES_YAML, "--resume"])
        core = mock.Mock()
        core.run.return_value = 10
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, 0)
        core.run.assert_called_once_with(resume=True)

    def test_dispatch_missing_task_yaml_returns_usage_exit(self) -> None:
        args = build_parser().parse_args(["--check-modules"])
        code = dispatch(args)
        self.assertEqual(code, 2)

    def test_dispatch_missing_file_is_clear_error(self) -> None:
        args = build_parser().parse_args(["/no/such/task.yaml", "--check-modules"])
        code = dispatch(args)
        self.assertEqual(code, 1)

    def test_dispatch_not_implemented_returns_exit_2(self) -> None:
        args = build_parser().parse_args([CHECK_MODULES_YAML])
        core = mock.Mock()
        core.run.side_effect = NotImplementedError("unsupported sampler")
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, 2)

    def test_main_does_not_import_v1_client(self) -> None:
        self.assertNotIn("jarvishep.client", sys.modules)
        self.assertNotIn("jarvishep", {name.split(".")[0] for name in sys.modules if name.startswith("jarvishep")})


if __name__ == "__main__":
    unittest.main()