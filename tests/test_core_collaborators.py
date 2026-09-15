#!/usr/bin/env python3
"""D25.3: Jarvis2Core is a façade over private runtime/scan/resume collaborators."""

from __future__ import annotations

import inspect
import json
import logging
import threading
import time
import unittest
from types import SimpleNamespace
from unittest import mock

from jarvishep2._resume_service import _ResumeService
from jarvishep2._runtime_supervisor import (
    REDIS_UNREACH_GRACE_SEC,
    _ARCHIVER_STALE_WARN_INTERVAL_SEC,
    _RuntimeSupervisor,
)
from jarvishep2._scan_driver import _ScanDriver
from jarvishep2.core import Jarvis2Core
from jarvishep2.redis_queue import make_fakeredis_queue


class CoreCollaboratorTests(unittest.TestCase):
    def test_core_holds_private_collaborators(self) -> None:
        core = Jarvis2Core()
        self.assertIsInstance(core._runtime, _RuntimeSupervisor)
        self.assertIsInstance(core._scan, _ScanDriver)
        self.assertIsInstance(core._resume, _ResumeService)
        self.assertIs(core._runtime._core, core)
        self.assertIs(core._scan._core, core)
        self.assertIs(core._resume._core, core)
        self.assertIsInstance(core._proc_board_lock, type(threading.Lock()))

    def test_new_without_init_still_lazily_builds_collaborators(self) -> None:
        core = Jarvis2Core.__new__(Jarvis2Core)
        self.assertIsInstance(core._get_scan(), _ScanDriver)
        self.assertIsInstance(core._get_runtime(), _RuntimeSupervisor)
        self.assertIsInstance(core._get_resume(), _ResumeService)

    def test_run_and_shutdown_signatures_unchanged(self) -> None:
        run_params = inspect.signature(Jarvis2Core.run).parameters
        self.assertEqual(
            list(run_params),
            [
                "self",
                "resume",
                "check_modules",
                "verify_golden",
                "write_run_summary",
                "check_timeout",
            ],
        )
        shutdown_params = inspect.signature(Jarvis2Core.shutdown).parameters
        self.assertEqual(list(shutdown_params), ["self", "wait", "write_run_summary"])

    def test_dead_archiver_process_is_restarted(self) -> None:
        class _DeadArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 4242

            def is_alive(self) -> bool:
                return False

            def join(self, timeout=None) -> None:
                return None

        core = Jarvis2Core()
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = logging.getLogger("test.archiver_restart")
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _DeadArchiver
        ):
            core.archiver = _DeadArchiver()
            core._runtime._ensure_archiver_alive()
        core.init_archiver.assert_called_once_with("/tmp/samples.hdf5")
        self.assertFalse(core._interrupt_requested)

    def test_archiver_restart_failure_requests_shutdown(self) -> None:
        class _DeadArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 7

            def is_alive(self) -> bool:
                return False

            def join(self, timeout=None) -> None:
                return None

        core = Jarvis2Core()
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = logging.getLogger("test.archiver_restart")
        core.init_archiver = mock.Mock(side_effect=RuntimeError("hdf5 locked"))
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _DeadArchiver
        ):
            core.archiver = _DeadArchiver()
            core._runtime._ensure_archiver_alive()
        self.assertTrue(core._interrupt_requested)

    def test_archiver_liveness_skips_live_process(self) -> None:
        class _LiveArchiver:
            db_path = "/tmp/samples.hdf5"

            def is_alive(self) -> bool:
                return True

        core = Jarvis2Core()
        core._shutdown_done = False
        core._interrupt_requested = False
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _LiveArchiver
        ):
            core.archiver = _LiveArchiver()
            core._runtime._ensure_archiver_alive()
        core.init_archiver.assert_not_called()

    def test_live_archiver_board_stays_non_empty(self) -> None:
        class _LiveArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 99

            def is_alive(self) -> bool:
                return True

        queue = make_fakeredis_queue()
        queue.publish_proc_board(
            "archiver",
            role="archiver",
            status="running",
            pid=99,
            records_written=3,
            ts=time.time(),
        )
        core = Jarvis2Core()
        core.redis = queue
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = logging.getLogger("test.archiver_board")
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _LiveArchiver
        ):
            core.archiver = _LiveArchiver()
            core._runtime._ensure_archiver_alive()
        board = queue.read_proc_board("archiver")
        self.assertTrue(board)
        self.assertEqual(int(board["pid"]), 99)
        self.assertEqual(board["status"], "running")
        core.init_archiver.assert_not_called()
        self.assertFalse(core._interrupt_requested)

    def test_stale_archiver_board_warns_but_does_not_kill(self) -> None:
        class _LiveArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 77

            def is_alive(self) -> bool:
                return True

        queue = make_fakeredis_queue()
        queue.publish_proc_board(
            "archiver",
            role="archiver",
            status="running",
            pid=77,
            ts=time.time() - 120.0,
        )
        core = Jarvis2Core()
        core.redis = queue
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = mock.Mock()
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _LiveArchiver
        ):
            core.archiver = _LiveArchiver()
            core._runtime._ensure_archiver_alive()
        core.init_archiver.assert_not_called()
        self.assertFalse(core._interrupt_requested)
        self.assertTrue(core.archiver.is_alive())
        core._logger.warning.assert_called()
        board = queue.read_proc_board("archiver")
        self.assertTrue(board)

    def test_missing_archiver_board_on_young_process_is_quiet(self) -> None:
        class _LiveArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 55

            def is_alive(self) -> bool:
                return True

        queue = make_fakeredis_queue()
        core = Jarvis2Core()
        core.redis = queue
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = mock.Mock()
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _LiveArchiver
        ):
            core.archiver = _LiveArchiver()
            core._runtime._ensure_archiver_alive()
        core.init_archiver.assert_not_called()
        self.assertFalse(core._interrupt_requested)
        core._logger.warning.assert_not_called()

    def test_missing_archiver_board_on_live_process_warns_without_kill(self) -> None:
        class _LiveArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 56

            def is_alive(self) -> bool:
                return True

        queue = make_fakeredis_queue()
        core = Jarvis2Core()
        core.redis = queue
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = mock.Mock()
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _LiveArchiver
        ):
            core.archiver = _LiveArchiver()
            limit = core._runtime._archiver_board_stale_limit_sec()
            core._archiver_observed = (56, time.monotonic() - limit - 1.0)
            core._runtime._ensure_archiver_alive()
        core.init_archiver.assert_not_called()
        self.assertFalse(core._interrupt_requested)
        self.assertTrue(core.archiver.is_alive())
        core._logger.warning.assert_called()

    def test_stale_archiver_warning_is_rate_limited(self) -> None:
        class _LiveArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 77

            def is_alive(self) -> bool:
                return True

        queue = make_fakeredis_queue()
        queue.publish_proc_board(
            "archiver",
            role="archiver",
            status="running",
            pid=77,
            ts=time.time() - 120.0,
        )
        core = Jarvis2Core()
        core.redis = queue
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = mock.Mock()
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _LiveArchiver
        ):
            core.archiver = _LiveArchiver()
            core._runtime._ensure_archiver_alive()
            core._runtime._ensure_archiver_alive()
            self.assertEqual(core._logger.warning.call_count, 1)
            core._archiver_stale_warned_at = (
                time.monotonic() - _ARCHIVER_STALE_WARN_INTERVAL_SEC
            )
            core._runtime._ensure_archiver_alive()
        self.assertEqual(core._logger.warning.call_count, 2)
        core.init_archiver.assert_not_called()

    def test_managed_redis_process_death_interrupts_without_empty_restart(self) -> None:
        queue = make_fakeredis_queue()
        queue.close = mock.Mock()
        process = mock.Mock()
        process.poll.return_value = 0
        process.pid = 4242
        managed = mock.Mock()
        managed.started_by_us = True
        managed.process = process
        managed.ensure = mock.Mock()
        managed.port = 6379
        managed.title = "Jarvis-Redis:scan"
        core = Jarvis2Core()
        core.redis = queue
        core._managed_redis = managed
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = mock.Mock()
        core.factory = mock.Mock()
        core._runtime._ensure_managed_redis_alive(now=0.0)
        self.assertTrue(core._interrupt_requested)
        managed.ensure.assert_not_called()
        queue.close.assert_called()
        core.factory.request_worker_shutdown.assert_called()
        core._logger.error.assert_called()
        logged = " ".join(str(call.args[0]) for call in core._logger.error.call_args_list)
        self.assertIn("managed redis-server died; refusing empty restart; use --resume", logged)

    def test_redis_ping_timeout_within_grace_warns_without_interrupt(self) -> None:
        queue = make_fakeredis_queue()
        queue.ping = mock.Mock(side_effect=TimeoutError("Timeout reading from socket"))
        queue.close = mock.Mock()
        process = mock.Mock()
        process.poll.return_value = None
        process.pid = 9
        managed = mock.Mock()
        managed.started_by_us = True
        managed.process = process
        managed.ensure = mock.Mock()
        managed.port = 6379
        managed.title = "Jarvis-Redis:scan"
        core = Jarvis2Core()
        core.redis = queue
        core._managed_redis = managed
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = mock.Mock()
        core._runtime._ensure_managed_redis_alive(now=0.0)
        core._runtime._ensure_managed_redis_alive(now=REDIS_UNREACH_GRACE_SEC - 0.1)
        self.assertFalse(core._interrupt_requested)
        queue.close.assert_not_called()
        managed.ensure.assert_not_called()
        self.assertTrue(core._logger.warning.called)
        board = queue.read_proc_board("core")
        self.assertNotEqual(board.get("scan_mode"), "stopping")

    def test_redis_ping_timeout_after_grace_interrupts_and_closes(self) -> None:
        queue = make_fakeredis_queue()
        queue.ping = mock.Mock(side_effect=TimeoutError("Timeout reading from socket"))
        queue.close = mock.Mock()
        process = mock.Mock()
        process.poll.return_value = None
        process.pid = 9
        managed = mock.Mock()
        managed.started_by_us = True
        managed.process = process
        managed.ensure = mock.Mock()
        managed.port = 6379
        managed.title = "Jarvis-Redis:scan"
        core = Jarvis2Core()
        core.redis = queue
        core._managed_redis = managed
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = mock.Mock()
        core._runtime._ensure_managed_redis_alive(now=0.0)
        self.assertFalse(core._interrupt_requested)
        core._runtime._ensure_managed_redis_alive(now=REDIS_UNREACH_GRACE_SEC)
        self.assertTrue(core._interrupt_requested)
        queue.close.assert_called()
        managed.ensure.assert_not_called()

    def test_lease_tick_does_not_overwrite_scan_mode_paused(self) -> None:
        queue = make_fakeredis_queue()
        core = Jarvis2Core()
        core.redis = queue
        core._shutdown_done = False
        core._interrupt_requested = False
        core._control_lock_owner = "owner-1"
        core._logger = logging.getLogger("test.lease_tick")
        core.sampler = SimpleNamespace(
            method="Random",
            _index=12,
            _maxp=100,
            _accepted_index=10,
            _seed=7,
        )
        queue.publish_proc_board("core", role="core", scan_mode="paused")
        with mock.patch.object(
            queue, "publish_proc_board", wraps=queue.publish_proc_board
        ) as published:
            core._runtime._publish_lease_proc_boards(redis_alive=1)
        for call in published.call_args_list:
            self.assertNotIn("scan_mode", call.kwargs)
        board = queue.read_proc_board("core")
        self.assertEqual(board.get("scan_mode"), "paused")
        self.assertEqual(int(board.get("redis_alive")), 1)
        sampler_status = json.loads(board["sampler_status"])
        self.assertEqual(sampler_status["method"], "Random")
        self.assertEqual(sampler_status["progress"]["current"], 12)

    def test_external_redis_board_has_empty_pid(self) -> None:
        queue = make_fakeredis_queue()
        managed = mock.Mock()
        managed.started_by_us = False
        managed.process = mock.Mock()
        managed.process.poll.return_value = 0
        managed.process.pid = 1
        managed.ensure = mock.Mock()
        managed.port = 6379
        managed.title = ""
        core = Jarvis2Core()
        core.redis = queue
        core._managed_redis = managed
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = logging.getLogger("test.external_redis")
        core._runtime._ensure_managed_redis_alive(now=0.0)
        self.assertFalse(core._interrupt_requested)
        managed.ensure.assert_not_called()
        board = queue.read_proc_board("redis")
        self.assertEqual(board.get("pid"), "")
        self.assertEqual(int(board.get("started_by_us")), 0)
        managed.process.poll.assert_not_called()

    def test_set_scan_mode_overlays_fuse_fields_without_rewriting_main(self) -> None:
        core = Jarvis2Core({"EnvReqs": {"V2": {"workers": 190}}})
        queue = make_fakeredis_queue()
        core.redis = queue
        core._runtime._publish_core_proc_board()
        before = queue.read_proc_board("core")
        self.assertEqual(int(before.get("workers_total") or 0), 190)
        self.assertNotIn("scan_mode", before)
        pid = before.get("pid")
        core._set_scan_mode(
            "paused",
            pause_reason="death_rate",
            workers_alive=150,
            workers_respawned=40,
            death_window_respawns=40,
            death_rate_1m=40,
            host="should-not-write",
            workers_total=1,
        )
        after = queue.read_proc_board("core")
        self.assertEqual(after["scan_mode"], "paused")
        self.assertEqual(after["pause_reason"], "death_rate")
        self.assertEqual(int(after.get("workers_alive") or 0), 150)
        self.assertEqual(int(after.get("workers_total") or 0), 190)
        self.assertEqual(after.get("pid"), pid)
        self.assertNotEqual(after.get("host"), "should-not-write")
        queue.publish_proc_board("core", redis_alive=1)
        still = queue.read_proc_board("core")
        self.assertEqual(still["scan_mode"], "paused")
        self.assertEqual(str(still.get("redis_alive")), "1")




if __name__ == "__main__":
    unittest.main()
