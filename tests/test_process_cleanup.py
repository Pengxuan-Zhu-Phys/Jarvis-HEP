#!/usr/bin/env python3
"""Tests for leftover Jarvis process detection / cleanup CLI helpers."""

from __future__ import annotations

import io
import os
import signal
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from jarvishep2.process_cleanup import (
    JarvisProcess,
    _command_is_jarvis,
    _process_role_and_detail,
    format_scan_table,
    format_process_table,
    list_active_scans,
    list_jarvis_processes,
    list_running_scans,
    list_zombie_processes,
    resolve_scan_reference,
)


class _StickyRefTestCase(unittest.TestCase):
    """Isolate sticky R1/R2 registry so tests do not touch ~/.jarvis."""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self._registry = Path(self._tmp.name) / "scan_refs.json"
        self._env = mock.patch.dict(
            os.environ, {"JARVIS_SCAN_REF_REGISTRY": str(self._registry)}
        )
        self._env.start()

    def tearDown(self) -> None:
        self._env.stop()
        self._tmp.cleanup()


class CommandMatchTests(unittest.TestCase):
    def test_matches_known_titles(self) -> None:
        self.assertTrue(_command_is_jarvis("Jarvis:EggBox_Bridson_V2"))
        self.assertTrue(_command_is_jarvis("Jarvis-Worker-3:scan"))
        self.assertTrue(_command_is_jarvis("Jarvis-Archiver:scan"))
        self.assertTrue(_command_is_jarvis("Jarvis-FileOperation"))
        self.assertTrue(_command_is_jarvis("Jarvis-Redis:scan"))

    def test_rejects_retired_jarvis2_titles(self) -> None:
        """Old Jarvis2* process titles are no longer managed by ps/kill."""
        self.assertFalse(_command_is_jarvis("Jarvis2:EggBox_Bridson_V2"))
        self.assertFalse(_command_is_jarvis("Jarvis2-Worker-3:scan"))
        self.assertFalse(_command_is_jarvis("Jarvis2-Archiver:scan"))
        self.assertFalse(_command_is_jarvis("Jarvis2-FileOperation"))
        self.assertFalse(_command_is_jarvis("Jarvis2"))

    def test_rejects_unrelated(self) -> None:
        self.assertFalse(_command_is_jarvis("python3 -m pytest"))
        self.assertFalse(_command_is_jarvis("redis-server"))
        self.assertFalse(_command_is_jarvis(""))
        self.assertFalse(_command_is_jarvis("Jarvis-HEP@oldscan"))
        self.assertFalse(_command_is_jarvis("JarvisHelper"))

    def test_excludes_jarvis_lit_titles_and_paths(self) -> None:
        self.assertFalse(_command_is_jarvis("Jarvis-Lit:serve"))
        self.assertFalse(_command_is_jarvis("JarvisLit:serve"))
        self.assertFalse(
            _command_is_jarvis(
                "/Users/p.zhu/Jarvis-Workshop/Jarvis-Lit/.venv/bin/python "
                "/Users/p.zhu/Jarvis-Workshop/Jarvis-Lit/.venv/bin/jlit serve"
            )
        )
        self.assertFalse(
            _command_is_jarvis(
                "/Users/p.zhu/Library/Developer/Xcode/DerivedData/"
                "JarvisLit.app/Contents/MacOS/JarvisLit"
            )
        )


class ListProcessesTests(_StickyRefTestCase):
    def test_groups_running_processes_into_scan_references(self) -> None:
        scans = list_running_scans(
            [
                JarvisProcess(pid=20, command="Jarvis:Beta"),
                JarvisProcess(pid=10, command="Jarvis:Alpha"),
                JarvisProcess(pid=11, command="Jarvis-Worker-0:Alpha"),
                JarvisProcess(pid=12, command="Jarvis-Redis:Alpha"),
                JarvisProcess(pid=13, command="Jarvis-FileOperation:Alpha"),
            ]
        )
        # Sticky R1/R2 by first-seen control PID order (lower PID first on cold start).
        self.assertEqual(
            [(scan.reference, scan.name) for scan in scans],
            [("R1", "Alpha"), ("R2", "Beta")],
        )
        self.assertEqual([proc.pid for proc in scans[0].processes], [10, 11, 12, 13])
        self.assertIs(resolve_scan_reference("R1", scans), scans[0])
        self.assertIs(resolve_scan_reference("10", scans), scans[0])
        self.assertIs(resolve_scan_reference("Beta", scans), scans[1])
        self.assertIn("Jarvis kill R1", format_scan_table(scans))

    def test_sticky_ref_does_not_shift_when_another_scan_appears(self) -> None:
        """Regression: ephemeral dense R# renumbered and made kill R2 hit the wrong task."""
        beta_only = list_active_scans([JarvisProcess(pid=20, command="Jarvis:Beta")])
        self.assertEqual(beta_only[0].reference, "R1")
        both = list_active_scans(
            [
                JarvisProcess(pid=10, command="Jarvis:Alpha"),  # new; must not steal R1
                JarvisProcess(pid=20, command="Jarvis:Beta"),
            ]
        )
        by_name = {scan.name: scan.reference for scan in both}
        self.assertEqual(by_name["Beta"], "R1")
        self.assertEqual(by_name["Alpha"], "R2")
        self.assertEqual(resolve_scan_reference("R1", both).name, "Beta")
        self.assertEqual(resolve_scan_reference("20", both).name, "Beta")

    def test_file_operation_has_a_scan_role(self) -> None:
        role, detail = _process_role_and_detail(
            JarvisProcess(pid=13, command="Jarvis-FileOperation:Alpha"), "Alpha"
        )
        self.assertEqual(role, "FileOperation")
        self.assertIn("save / copy / delete", detail)

    def test_unscoped_components_use_zp_not_scan_references(self) -> None:
        processes = [
            JarvisProcess(pid=10, command="Jarvis"),
            JarvisProcess(pid=11, command="Jarvis-Worker-0"),
            JarvisProcess(pid=12, command="Jarvis-Archiver"),
            JarvisProcess(pid=13, command="Jarvis-FileOperation"),
            JarvisProcess(pid=14, command="Jarvis-Redis@6379/0"),
            JarvisProcess(pid=15, command="Jarvis-Worker-1:Alpha"),
            JarvisProcess(pid=16, command="JarvisHelper"),
        ]
        self.assertEqual(list_running_scans(processes)[0].name, "Alpha")
        self.assertEqual(list_active_scans(processes), [])
        zombies = list_zombie_processes(processes)
        self.assertEqual([proc.pid for proc in zombies], [10, 11, 12, 13, 14, 15])

    def test_orphan_with_stale_scan_suffix_is_zp_not_r_reference(self) -> None:
        processes = [JarvisProcess(pid=19264, command="Jarvis-Archiver:scan")]
        self.assertEqual(list_active_scans(processes), [])
        self.assertEqual(
            [proc.pid for proc in list_zombie_processes(processes)],
            [19264],
        )

    def test_stale_groups_do_not_consume_active_r_references(self) -> None:
        processes = [
            JarvisProcess(pid=10, command="Jarvis-Archiver:Alpha"),
            JarvisProcess(pid=20, command="Jarvis:Beta"),
            JarvisProcess(pid=21, command="Jarvis-Worker-0:Beta"),
        ]
        active = list_active_scans(processes)
        self.assertEqual([(scan.reference, scan.name) for scan in active], [("R1", "Beta")])
        self.assertEqual([proc.pid for proc in list_zombie_processes(processes)], [10])

    def test_zp_is_a_fixed_extra_row_in_process_table(self) -> None:
        zombies = [JarvisProcess(pid=13, command="Jarvis-FileOperation")]
        text = format_scan_table([], zombies)
        self.assertIn("ZP", text)
        self.assertIn("Zambie Process", text)
        self.assertIn("Jarvis kill ZP -y", text)

    def test_rejects_duplicate_live_scan_controller(self) -> None:
        from jarvishep2.process_cleanup import ensure_scan_name_available

        fake = [JarvisProcess(pid=10, command="Jarvis:Alpha")]
        with mock.patch("jarvishep2.process_cleanup.list_jarvis_processes", return_value=fake):
            with self.assertRaisesRegex(RuntimeError, "already running"):
                ensure_scan_name_available("Alpha")

    def test_requires_cleanup_when_only_orphan_worker_remains(self) -> None:
        from jarvishep2.process_cleanup import ensure_scan_name_available

        fake = [JarvisProcess(pid=10, command="Jarvis-Worker-0:Alpha")]
        with mock.patch("jarvishep2.process_cleanup.list_jarvis_processes", return_value=fake):
            with self.assertRaisesRegex(RuntimeError, "stale Jarvis runtime"):
                ensure_scan_name_available("Alpha")

    def test_resume_cleans_same_name_orphan_runtime(self) -> None:
        from jarvishep2.process_cleanup import ensure_scan_name_available

        fake = [JarvisProcess(pid=10, command="Jarvis-Worker-0:Alpha")]
        with (
            mock.patch(
                "jarvishep2.process_cleanup.list_jarvis_processes",
                side_effect=[fake, []],
            ),
            mock.patch("jarvishep2.process_cleanup.kill_jarvis_processes") as kill,
        ):
            kill.return_value = {
                "signaled": [10],
                "killed": [],
                "missing": [],
                "failed": [],
            }
            ensure_scan_name_available("Alpha", cleanup_stale=True)
        kill.assert_called_once()

    def test_parses_ps_output_and_skips_self(self) -> None:
        fake = (
            "  111 python3 -m pytest\n"
            "  222 Jarvis:DemoScan\n"
            "  333 Jarvis-Worker-0:DemoScan\n"
            "  444 Jarvis-Redis:DemoScan\n"
            "  555 /Users/p.zhu/Jarvis-Workshop/Jarvis-Lit/.venv/bin/jlit serve\n"
        )
        with mock.patch("jarvishep2.process_cleanup.os.getpid", return_value=999):
            with mock.patch("jarvishep2.process_cleanup.subprocess.run") as run:
                run.return_value = mock.Mock(returncode=0, stdout=fake, stderr="")
                procs = list_jarvis_processes()
        titles = [p.command for p in procs]
        self.assertEqual(
            titles,
            [
                "Jarvis:DemoScan",
                "Jarvis-Worker-0:DemoScan",
                "Jarvis-Redis:DemoScan",
            ],
        )
        self.assertEqual([p.pid for p in procs], [222, 333, 444])

    def test_format_empty(self) -> None:
        self.assertIn("No running Jarvis processes", format_process_table([]))

    def test_format_rows(self) -> None:
        text = format_process_table(
            [JarvisProcess(pid=42, command="Jarvis:x")]
        )
        self.assertIn("42", text)
        self.assertIn("Jarvis:x", text)


class KillOrderTests(unittest.TestCase):
    def test_kill_sends_term_then_kill_in_dependency_order(self) -> None:
        from jarvishep2.process_cleanup import kill_jarvis_processes

        procs = [
            JarvisProcess(pid=10, command="Jarvis-Redis:s"),
            JarvisProcess(pid=20, command="Jarvis:s"),
            JarvisProcess(pid=30, command="Jarvis-Worker-0:s"),
            JarvisProcess(pid=40, command="Jarvis-Archiver:s"),
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


class CliPsKillTests(_StickyRefTestCase):
    def test_ps_subcommand_lists_scan_choices(self) -> None:
        from jarvishep2.client import main

        fake = [JarvisProcess(pid=99, command="Jarvis:Demo")]
        with mock.patch(
            "jarvishep2.process_cleanup.list_jarvis_processes",
            return_value=fake,
        ):
            output = io.StringIO()
            with redirect_stdout(output):
                code = main(["ps"])
        self.assertEqual(code, 0)
        self.assertIn("R1", output.getvalue())
        self.assertIn("Demo", output.getvalue())

    def test_ps_reference_filters_to_one_scan(self) -> None:
        from jarvishep2.client import main

        fake = [
            JarvisProcess(pid=10, command="Jarvis:Alpha"),
            JarvisProcess(pid=20, command="Jarvis:Beta"),
        ]
        with mock.patch("jarvishep2.process_cleanup.list_jarvis_processes", return_value=fake):
            output = io.StringIO()
            with redirect_stdout(output):
                code = main(["ps", "R1"])
        self.assertEqual(code, 0)
        self.assertIn("Alpha", output.getvalue())
        self.assertNotIn("Beta", output.getvalue())

    def test_ps_lists_and_selects_zp(self) -> None:
        from jarvishep2.client import main

        fake = [
            JarvisProcess(pid=10, command="Jarvis:Alpha"),
            JarvisProcess(pid=20, command="Jarvis-FileOperation"),
        ]
        with mock.patch("jarvishep2.process_cleanup.list_jarvis_processes", return_value=fake):
            choices = io.StringIO()
            with redirect_stdout(choices):
                self.assertEqual(main(["ps"]), 0)
            selected = io.StringIO()
            with redirect_stdout(selected):
                self.assertEqual(main(["ps", "ZP"]), 0)
        self.assertIn("ZP", choices.getvalue())
        self.assertIn("20", selected.getvalue())
        self.assertIn("FileOperation", selected.getvalue())
        self.assertNotIn("Alpha", selected.getvalue())

    def test_kill_without_reference_only_lists_choices(self) -> None:
        from jarvishep2.client import main

        fake = [JarvisProcess(pid=99, command="Jarvis:Demo")]
        with mock.patch("jarvishep2.process_cleanup.list_jarvis_processes", return_value=fake):
            with mock.patch("jarvishep2.process_cleanup.kill_jarvis_processes") as kill_fn:
                self.assertEqual(main(["kill"]), 0)
        kill_fn.assert_not_called()

    def test_kill_aborts_without_confirm(self) -> None:
        from jarvishep2.client import main

        fake = [JarvisProcess(pid=99, command="Jarvis:x")]
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

        fake = [JarvisProcess(pid=99, command="Jarvis:x")]
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
                    code = main(["kill", "R1", "--yes"])
        self.assertEqual(code, 0)
        kill_fn.assert_called_once()

    def test_kill_zp_terminates_every_unscoped_component_only(self) -> None:
        from jarvishep2.client import main

        fake = [
            JarvisProcess(pid=10, command="Jarvis:Alpha"),
            JarvisProcess(pid=20, command="Jarvis-FileOperation"),
            JarvisProcess(pid=21, command="Jarvis-Worker-0"),
            JarvisProcess(pid=22, command="Jarvis-Redis"),
            JarvisProcess(pid=23, command="JarvisHelper"),
            JarvisProcess(pid=24, command="Jarvis-Archiver:stale"),
        ]
        with mock.patch(
            "jarvishep2.process_cleanup.list_jarvis_processes",
            side_effect=[fake, []],
        ):
            with mock.patch(
                "jarvishep2.process_cleanup.kill_jarvis_processes",
                return_value={
                    "signaled": [20, 21, 22, 24],
                    "killed": [],
                    "missing": [],
                    "failed": [],
                },
            ) as kill_fn:
                self.assertEqual(main(["kill", "ZP", "--yes"]), 0)
        targets = kill_fn.call_args.args[0]
        self.assertEqual([proc.pid for proc in targets], [20, 21, 22, 24])

    def test_empty_zp_is_a_successful_noop(self) -> None:
        from jarvishep2.client import main

        with mock.patch("jarvishep2.process_cleanup.list_jarvis_processes", return_value=[]):
            output = io.StringIO()
            with redirect_stdout(output):
                self.assertEqual(main(["kill", "zp", "-y"]), 0)
        self.assertIn("No ZP", output.getvalue())

    def test_format_says_running_not_leftover(self) -> None:
        text = format_process_table([JarvisProcess(pid=1, command="Jarvis:a")])
        self.assertIn("Running Jarvis processes", text)
        self.assertNotIn("leftover", text.lower())


if __name__ == "__main__":
    unittest.main()
