#!/usr/bin/env python3
"""Tests for leftover Jarvis process detection / cleanup CLI helpers."""

from __future__ import annotations

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
        self.assertIn("No leftover", format_process_table([]))

    def test_format_rows(self) -> None:
        text = format_process_table(
            [JarvisProcess(pid=42, command="Jarvis2:x")]
        )
        self.assertIn("42", text)
        self.assertIn("Jarvis2:x", text)


if __name__ == "__main__":
    unittest.main()
