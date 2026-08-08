#!/usr/bin/env python3
"""FileOperationService: HEP-owned SAMPLE save/copy/delete process."""

from __future__ import annotations

import os
import select
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

from jarvishep2.file_operation_service import (
    FileOperationService,
    _file_operation_main,
    _signal_process_tree,
)
from jarvishep2.file_ops import apply_io_save_policy, save_io_copy
from jarvishep2.proc_title import file_operation_title


def _pid_exists(pid: int) -> bool:
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class SaveIoCopyTests(unittest.TestCase):
    def test_save_true_copies_into_sample_dir(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "runtime" / "input.json"
            src.parent.mkdir(parents=True)
            src.write_text('{"x": 1}\n', encoding="utf-8")
            sample = Path(tmp) / "SAMPLE" / "uuid-1"
            dest = save_io_copy(
                str(src),
                str(sample),
                module="EggBox",
                save=True,
            )
            self.assertIsNotNone(dest)
            assert dest is not None
            self.assertTrue(dest.endswith("input.json@EggBox"))
            self.assertTrue(os.path.isfile(dest))
            self.assertEqual(Path(dest).read_text(encoding="utf-8"), '{"x": 1}\n')

    def test_output_default_temp_when_save_false(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "runtime" / "output.json"
            src.parent.mkdir(parents=True)
            src.write_text('{"z": 2}\n', encoding="utf-8")
            sample = Path(tmp) / "SAMPLE" / "uuid-1"
            dest = apply_io_save_policy(
                source_path=str(src),
                sample_save_dir=str(sample),
                module="EggBox",
                spec={"save": False},
                direction="output",
            )
            self.assertIsNotNone(dest)
            assert dest is not None
            self.assertIn(os.path.join(".temp", "output.json@EggBox"), dest.replace("/", os.sep))
            self.assertTrue(os.path.isfile(dest))

    def test_input_no_copy_when_save_false(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "runtime" / "input.json"
            src.parent.mkdir(parents=True)
            src.write_text("{}", encoding="utf-8")
            sample = Path(tmp) / "SAMPLE" / "uuid-1"
            dest = apply_io_save_policy(
                source_path=str(src),
                sample_save_dir=str(sample),
                module="EggBox",
                spec={"save": False},
                direction="input",
            )
            self.assertIsNone(dest)


class FileOperationServiceTests(unittest.TestCase):
    def test_file_operation_title_includes_scan_name(self) -> None:
        self.assertEqual(file_operation_title(scan_name="Demo"), "Jarvis-FileOperation:Demo")

    def test_inline_save_and_delete(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "out.json"
            src.write_text('{"z": 1}\n', encoding="utf-8")
            sample = Path(tmp) / "SAMPLE"
            service = FileOperationService.start(mode="inline", delete_method="shutil")
            try:
                dest = service.save_io(
                    source_path=str(src),
                    sample_save_dir=str(sample),
                    module="Calc",
                    spec={"save": True, "name": "oup"},
                    direction="output",
                )
                self.assertIsNotNone(dest)
                assert dest is not None
                self.assertTrue(os.path.isfile(dest))
                service.delete([dest], missing_ok=False)
                self.assertFalse(os.path.lexists(dest))
            finally:
                service.shutdown()

    def test_process_mode_save_io(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "out.json"
            src.write_text('{"z": 3}\n', encoding="utf-8")
            sample = Path(tmp) / "SAMPLE"
            service = FileOperationService.start(mode="process", delete_method="shutil")
            try:
                if hasattr(os, "getpgid"):
                    assert service.pid is not None
                    deadline = time.monotonic() + 2.0
                    while (
                        time.monotonic() < deadline
                        and os.getpgid(service.pid) != service.pid
                    ):
                        time.sleep(0.02)
                    self.assertEqual(os.getpgid(service.pid), service.pid)
                dest = service.save_io(
                    source_path=str(src),
                    sample_save_dir=str(sample),
                    module="Calc",
                    spec={"save": True},
                    direction="output",
                )
                self.assertIsNotNone(dest)
                assert dest is not None
                self.assertTrue(os.path.isfile(dest))
                self.assertTrue(dest.endswith("out.json@Calc"))
            finally:
                service.shutdown()

    def test_orphan_process_exits_when_worker_parent_is_gone(self) -> None:
        request = mock.Mock()
        response = mock.Mock()
        with mock.patch("jarvishep2.file_operation_service.os.getppid", return_value=1):
            _file_operation_main(request, response, "shutil", "Demo", owner_pid=9876)
        request.get.assert_not_called()

    def test_client_fails_fast_when_file_operation_child_crashes(self) -> None:
        service = FileOperationService.start(
            mode="process",
            delete_method="shutil",
            scan_name="crash-test",
        )
        try:
            assert service._process is not None
            service._process.kill()
            service._process.join(timeout=2.0)
            started = time.monotonic()
            with self.assertRaisesRegex(RuntimeError, "exited before returning"):
                service.delete(["/path/that/will/not/be-used"])
            self.assertLess(time.monotonic() - started, 2.0)
        finally:
            service.shutdown(timeout=0.1)

    def test_forced_shutdown_signals_entire_file_operation_process_group(self) -> None:
        process = mock.Mock(pid=321)
        with (
            mock.patch("jarvishep2.file_operation_service.os.getpgid", return_value=321),
            mock.patch("jarvishep2.file_operation_service.os.killpg") as killpg,
            mock.patch("jarvishep2.file_operation_service.os.kill") as kill_one,
        ):
            _signal_process_tree(process, signal.SIGKILL)
        killpg.assert_called_once_with(321, signal.SIGKILL)
        kill_one.assert_not_called()

    @unittest.skipUnless(hasattr(os, "mkfifo") and hasattr(os, "kill"), "POSIX process test")
    def test_busy_file_operation_exits_when_worker_is_sigkilled(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fifo = os.path.join(tmp, "blocked-input")
            os.mkfifo(fifo)
            sample_dir = os.path.join(tmp, "SAMPLE")
            script = (
                "import sys; "
                "from jarvishep2.file_operation_service import FileOperationService; "
                "service=FileOperationService.start(mode='process', delete_method='shutil', "
                "scan_name='blocked-owner-test'); "
                "print(service.pid, flush=True); "
                "service.save_io(source_path=sys.argv[1], sample_save_dir=sys.argv[2], "
                "module='Calc', spec={'save': True}, direction='output')"
            )
            owner = subprocess.Popen(
                [sys.executable, "-c", script, fifo, sample_dir],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            file_operation_pid: int | None = None
            try:
                assert owner.stdout is not None
                readable, _, _ = select.select([owner.stdout], [], [], 10.0)
                self.assertTrue(readable, "owner did not publish FileOperation PID")
                file_operation_pid = int(owner.stdout.readline().strip())
                self.assertTrue(_pid_exists(file_operation_pid))
                os.kill(owner.pid, signal.SIGKILL)
                owner.wait(timeout=5.0)

                deadline = time.monotonic() + 8.0
                while time.monotonic() < deadline and _pid_exists(file_operation_pid):
                    time.sleep(0.1)
                self.assertFalse(
                    _pid_exists(file_operation_pid),
                    "busy FileOperation survived after its Worker owner was SIGKILLed",
                )
            finally:
                if owner.poll() is None:
                    owner.kill()
                    owner.wait(timeout=2.0)
                if file_operation_pid is not None and _pid_exists(file_operation_pid):
                    os.kill(file_operation_pid, signal.SIGKILL)
                if owner.stdout is not None:
                    owner.stdout.close()
                if owner.stderr is not None:
                    owner.stderr.close()


if __name__ == "__main__":
    unittest.main()
