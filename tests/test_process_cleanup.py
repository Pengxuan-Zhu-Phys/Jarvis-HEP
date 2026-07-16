#!/usr/bin/env python3
"""Tests for leftover Jarvis process detection / cleanup CLI helpers."""

from __future__ import annotations

import signal
import unittest
from unittest import mock

from jarvishep2.process_cleanup import (
    JarvisProcess,
    _command_is_jarvis,
    format_process_table,
    list_jarvis_processes,
)


class CommandMatchTests(unittest.TestCase):
    def test_matches_known_titles(self) -> None:
        self.assertTrue(_command_is_jarvis("Jarvis2:EggBox_Bridson_V2"))
        self.assertTrue(_command_is_jarvis("Jarvis2-Worker-3:scan"))
        self.assertTrue(_command_is_jarvis("Jarvis2-Archiver:scan"))
        self.assertTrue(_command_is_jarvis("Jarvis2-FileOperation"))
        self.assertTrue(_command_is_jarvis("Jarvis-Redis:scan"))

    def test_rejects_unrelated(self) -> None:
        self.assertFalse(_command_is_jarvis("python3 -m pytest"))
        self.assertFalse(_command_is_jarvis("redis-server"))
        self.assertFalse(_command_is_jarvis(""))


class ListProcessesTests(unittest.TestCase):
    def test_parses_ps_output_and_skips_self(self) -> None:
        fake = (
            "  111 python3 -m pytest\n"
            "  222 Jarvis2:DemoScan\n"
            "  333 Jarvis2-Worker-0:DemoScan\n"
            "  444 Jarvis-Redis:DemoScan\n"
        )
        with mock.patch("jarvishep2.process_cleanup.os.getpid", return_value=999):
            with mock.patch("jarvishep2.process_cleanup.subprocess.run") as run:
                run.return_value = mock.Mock(returncode=0, stdout=fake, stderr="")
                procs = list_jarvis_processes()
        titles = [p.command for p in procs]
        self.assertEqual(
            titles,
            [
                "Jarvis2:DemoScan",
                "Jarvis2-Worker-0:DemoScan",
                "Jarvis-Redis:DemoScan",
            ],
        )
        self.assertEqual([p.pid for p in procs], [222, 333, 444])

    def test_format_empty(self) -> None:
        self.assertIn("No running Jarvis processes", format_process_table([]))

    def test_format_rows(self) -> None:
        text = format_process_table(
            [JarvisProcess(pid=42, command="Jarvis2:x")]
        )
        self.assertIn("42", text)
        self.assertIn("Jarvis2:x", text)


class KillOrderTests(unittest.TestCase):
    def test_kill_sends_term_then_kill_in_dependency_order(self) -> None:
        from jarvishep2.process_cleanup import kill_jarvis_processes

        procs = [
            JarvisProcess(pid=10, command="Jarvis-Redis:s"),
            JarvisProcess(pid=20, command="Jarvis2:s"),
            JarvisProcess(pid=30, command="Jarvis2-Worker-0:s"),
            JarvisProcess(pid=40, command="Jarvis2-Archiver:s"),
        ]
        calls: list[tuple[int, int]] = []

        def fake_kill(pid: int, sig: int) -> None:
            # 0 = existence probe
            if sig == 0:
                return
            calls.append((pid, int(sig)))

        with mock.patch("jarvishep2.process_cleanup.os.kill", side_effect=fake_kill):
            with mock.patch("jarvishep2.process_cleanup.time.sleep", return_value=None):
                result = kill_jarvis_processes(procs, sigterm_grace_sec=0.01, force=True)

        # TERM order: worker → archiver → control → redis
        term_pids = [pid for pid, sig in calls if sig == signal.SIGTERM]
        self.assertEqual(term_pids, [30, 40, 20, 10])
        self.assertEqual(result["signaled"], [30, 40, 20, 10])
        # All still "alive" under fake_kill → SIGKILL same set
        kill_pids = [pid for pid, sig in calls if sig == signal.SIGKILL]
        self.assertEqual(set(kill_pids), {30, 40, 20, 10})


class CliPsKillTests(unittest.TestCase):
    def test_ps_subcommand_list_only(self) -> None:
        from jarvishep2.client import main

        with mock.patch(
            "jarvishep2.process_cleanup.list_jarvis_processes",
            return_value=[],
        ):
            code = main(["ps"])
        self.assertEqual(code, 0)

    def test_kill_aborts_without_confirm(self) -> None:
        from jarvishep2.client import main

        fake = [JarvisProcess(pid=99, command="Jarvis2:x")]
        with mock.patch(
            "jarvishep2.process_cleanup.list_jarvis_processes",
            return_value=fake,
        ):
            with mock.patch(
                "jarvishep2.process_cleanup.confirm_kill",
                return_value=False,
            ):
                with mock.patch(
                    "jarvishep2.process_cleanup.kill_jarvis_processes"
                ) as kill_fn:
                    code = main(["kill"])
        self.assertEqual(code, 0)
        kill_fn.assert_not_called()

    def test_kill_yes_skips_confirm(self) -> None:
        from jarvishep2.client import main

        fake = [JarvisProcess(pid=99, command="Jarvis2:x")]
        with mock.patch(
            "jarvishep2.process_cleanup.list_jarvis_processes",
            side_effect=[fake, []],
        ):
            with mock.patch(
                "jarvishep2.process_cleanup.kill_jarvis_processes",
                return_value={
                    "signaled": [99],
                    "killed": [],
                    "missing": [],
                    "failed": [],
                },
            ) as kill_fn:
                code = main(["kill", "--yes"])
        self.assertEqual(code, 0)
        kill_fn.assert_called_once()

    def test_format_says_running_not_leftover(self) -> None:
        text = format_process_table([JarvisProcess(pid=1, command="Jarvis2:a")])
        self.assertIn("Running Jarvis processes", text)
        self.assertNotIn("leftover", text.lower())


if __name__ == "__main__":
    unittest.main()
