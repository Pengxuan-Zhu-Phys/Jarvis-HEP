"""Native spawn/Redis acceptance for the independent checkpoint component."""

import logging
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
from unittest import mock

import psutil
import pytest

from jarvishep2.process_cleanup import (
    JarvisProcess, _command_is_jarvis, _kill_order_key, _process_role_and_detail,
    list_jarvis_processes, list_running_scans, list_zombie_processes,
)
from jarvishep2.runtime.checkpoint_process import CheckpointProcess
from jarvishep2.sampling.sample_checkpoint import atomic_pickle, read_pickle


def until(predicate, *, timeout=10):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = predicate()
        if value:
            return value
        time.sleep(0.02)
    pytest.fail("checkpoint condition timed out")


class SerializationPid:
    def __getstate__(self):
        return {"serialized_pid": os.getpid()}


class ChildWriteFailure:
    def __init__(self):
        self.parent_pid = os.getpid()

    def __getstate__(self):
        if os.getpid() != self.parent_pid:
            raise OSError("injected child write failure")
        return vars(self)


class PausedChildWrite:
    def __init__(self, root):
        self.parent_pid = os.getpid()
        self.started = str(root / "started")
        self.release = str(root / "release")

    def __getstate__(self):
        if os.getpid() != self.parent_pid:
            Path(self.started).touch()
            deadline = time.monotonic() + 20
            while not Path(self.release).exists() and time.monotonic() < deadline:
                time.sleep(0.01)
        return vars(self)


class KilledChildWrite:
    def __init__(self):
        self.parent_pid = os.getpid()

    def __getstate__(self):
        if os.getpid() != self.parent_pid:
            os._exit(23)
        return vars(self)


def component(tmp_path, **kwargs):
    writer = CheckpointProcess(str(tmp_path / "runtime.pickle"),
                               logging.getLogger("test.checkpoint"), **kwargs)
    assert writer.start()
    return writer


def test_sidecar_bytes_land_beside_the_pickle_and_bad_names_do_not(tmp_path):
    from jarvishep2.runtime.checkpoint_process import _write_payload_sidecars

    payload = {
        "sampler_state": {"native_sampler_blob": b"engine-bytes", "next_sample_index": 3},
        "sidecars": [
            {"name": "nested_engine.pkl", "field": "sampler_state.native_sampler_blob"},
            {"name": "../escape.pkl", "field": "sampler_state.native_sampler_blob"},
            {"name": "empty.pkl", "field": "missing"},
        ],
    }
    _write_payload_sidecars(str(tmp_path / "runtime.pickle"), payload, logging.getLogger("test.checkpoint"))
    assert (tmp_path / "nested_engine.pkl").read_bytes() == b"engine-bytes"
    assert not (tmp_path.parent / "escape.pkl").exists()
    assert not (tmp_path / "empty.pkl").exists()

    writer = component(tmp_path, interval=30)
    try:
        assert writer.submit(payload)
        assert writer.wait()
        landed = read_pickle(writer.path)
        assert landed["sampler_state"]["native_sampler_blob"] == b"engine-bytes"
        assert (tmp_path / "nested_engine.pkl").read_bytes() == b"engine-bytes"
    finally:
        assert writer.finish()


def test_spawn_child_owns_serialization_and_periodic_requests(tmp_path):
    writer = component(tmp_path, interval=0.05, scan_name="acceptance")
    try:
        assert writer.pid != os.getpid()
        assert writer._process._start_method == "spawn"
        assert writer._process.daemon is False
        until(writer.requested)
        assert not Path(writer.path).exists()
        assert writer.submit({"probe": SerializationPid(), "sampler_state": {"next_sample_index": 7}})
        assert writer.wait()
        payload = read_pickle(writer.path)
        assert payload["probe"].serialized_pid == writer.pid
        assert writer._written.value == 1
        until(writer.requested)
    finally:
        assert writer.finish()
    assert not writer.is_alive()
    assert writer.finish()  # Explicit shutdown is idempotent.


def test_busy_process_never_accumulates_or_blocks_new_snapshots(tmp_path):
    writer = component(tmp_path)
    try:
        assert writer.submit({"probe": PausedChildWrite(tmp_path)})
        until(lambda: (tmp_path / "started").exists())
        start = time.monotonic()
        for i in range(100):
            assert not writer.submit({"index": i})
        assert time.monotonic() - start < 1
        assert not writer.available()
        (tmp_path / "release").touch()
        assert writer.wait()
        assert writer._written.value == 1
    finally:
        (tmp_path / "release").touch()
        assert writer.finish()


def test_write_failure_preserves_previous_file_and_process_stays_available(tmp_path):
    path = tmp_path / "runtime.pickle"
    atomic_pickle(str(path), {"previous": 1})
    original = path.read_bytes()
    on_saved = mock.Mock()
    writer = component(tmp_path, on_saved=on_saved)
    barrier = {"generation": 2, "index": 8, "next_sample_index": 7,
               "session": {"retained": 0}}
    try:
        assert writer.submit({"probe": ChildWriteFailure(), "sampler_state": barrier})
        assert writer.wait()
        assert writer._failures.value == 1
        assert path.read_bytes() == original
        assert writer.last_saved_index is None
        on_saved.assert_not_called()
        assert writer.available()
        assert writer.submit({"sampler_state": barrier})
        assert writer.wait()
        assert read_pickle(str(path)) == {"sampler_state": barrier}
        assert writer.last_saved_index == 7
        on_saved.assert_called_once_with(str(path), "generation barrier", progress="generation=2 proposals=8")
    finally:
        assert writer.finish()


def test_dead_child_does_not_stall_submit_or_shutdown(tmp_path):
    writer = component(tmp_path)
    try:
        assert writer.submit({"probe": KilledChildWrite()})
        until(lambda: not writer.is_alive())
        start = time.monotonic()
        assert not writer.available()
        assert not writer.submit({"next": 2})
        assert not writer.wait(timeout=0.1)
        assert writer.finish(timeout=0.1)
        assert time.monotonic() - start < 2
    finally:
        writer.finish(timeout=0.1)


def test_stalled_io_is_killed_with_bounded_shutdown(tmp_path):
    writer = component(tmp_path)
    try:
        assert writer.submit({"probe": PausedChildWrite(tmp_path)})
        until(lambda: (tmp_path / "started").exists())
        start = time.monotonic()
        assert writer.finish(timeout=0.05)
        assert time.monotonic() - start < 4
        assert not writer.is_alive()
        assert not Path(writer.path).exists()
    finally:
        writer.finish(timeout=0.1)


def test_start_failure_does_not_raise_or_leave_a_child(tmp_path):
    writer = CheckpointProcess(str(tmp_path / "runtime.pickle"), mock.Mock())
    with mock.patch.object(writer._process, "start", side_effect=OSError("spawn failed")):
        assert writer.start() is False
    assert not writer.available()
    assert writer.finish()


@pytest.fixture
def native_redis(tmp_path):
    from jarvishep2.redis_queue import RedisQueue

    executable = shutil.which("redis-server")
    if executable is None:
        pytest.skip("redis-server is required for native telemetry acceptance")
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    server = subprocess.Popen([
        executable, "--bind", "127.0.0.1", "--port", str(port),
        "--save", "", "--appendonly", "no", "--dir", str(tmp_path),
        "--logfile", str(tmp_path / "redis.log"),
    ], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    queue = RedisQueue({"host": "127.0.0.1", "port": port, "db": 0,
                        "socket_connect_timeout": 0.2})
    try:
        def ready():
            try:
                queue.connect()
                return True
            except Exception:
                queue.close()
                return False
        until(ready)
        yield queue
    finally:
        queue.close()
        server.terminate()
        try:
            server.wait(timeout=5)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait(timeout=5)


def test_native_redis_board_and_process_discovery(tmp_path, native_redis, monkeypatch):
    from jarvishep2.dashboard import MonitorView

    monkeypatch.setenv("JARVIS_SCAN_REF_REGISTRY", str(tmp_path / "scan-refs.json"))
    owner = "checkpoint-acceptance-owner"
    assert native_redis.claim_control_lock(owner)
    writer = component(tmp_path, interval=0.05, scan_name="checkpoint-acceptance",
                       redis_config=native_redis.connection_config(), control_owner=owner,
                       logs_dir=str(tmp_path / "logs"))
    try:
        board = until(lambda: native_redis.read_proc_board("checkpoint"))
        assert int(board["pid"]) == writer.pid
        assert board["role"] == "checkpoint"
        assert board["owner"] == owner
        assert native_redis.r.ttl("hep:proc:checkpoint") > 0
        until(writer.requested)
        assert writer.submit({"sampler_state": {"next_sample_index": 13}})
        assert writer.wait()
        until(lambda: int(native_redis.read_proc_board("checkpoint").get("snapshots_written", 0)) == 1)
        board = native_redis.read_proc_board("checkpoint")
        assert int(board["last_index"]) == 13
        assert float(board["last_saved"]) > 0
        assert not board["last_error"]
        assert not native_redis.r.exists("hep:task_queue", "hep:archive_queue")
        view = MonitorView.from_snapshot({"proc_checkpoint": board})
        assert view.proc_checkpoint == board
        processes = until(lambda: [p for p in list_jarvis_processes() if p.pid == writer.pid])
        assert processes[0].command.startswith("Jarvis-Checkpoint:checkpoint-acceptance")
        assert (tmp_path / "logs/checkpoint.log").is_file()
        assert "Checkpoint saved to" in (tmp_path / "logs/checkpoint.log").read_text()
        assert list_running_scans(processes)[0].name == "checkpoint-acceptance"
    finally:
        assert writer.finish()
    assert native_redis.read_proc_board("checkpoint")["status"] == "stopped"


def test_stolen_lease_stops_only_the_checkpoint_component(tmp_path, native_redis):
    assert native_redis.claim_control_lock("first-owner")
    writer = component(tmp_path, redis_config=native_redis.connection_config(), control_owner="first-owner")
    try:
        until(lambda: native_redis.read_proc_board("checkpoint"))
        native_redis.r.set("hep:control:lock", "new-owner")
        until(lambda: not writer.is_alive())
        assert not writer.available()
        assert native_redis.get_control_lock_owner() == "new-owner"
    finally:
        assert writer.finish()


def test_ps_kill_recognizes_scoped_and_orphan_checkpoint_components(tmp_path, monkeypatch):
    monkeypatch.setenv("JARVIS_SCAN_REF_REGISTRY", str(tmp_path / "scan-refs.json"))
    process = JarvisProcess(pid=99, command="Jarvis-Checkpoint:scan")
    assert _command_is_jarvis(process.command)
    assert list_running_scans([process])[0].name == "scan"
    assert list_zombie_processes([process]) == [process]
    orphan = JarvisProcess(pid=100, command="Jarvis-Checkpoint")
    assert list_zombie_processes([orphan]) == [orphan]
    assert _process_role_and_detail(process, "scan")[0] == "Checkpoint"
    assert _kill_order_key(process) < _kill_order_key(JarvisProcess(pid=1, command="Jarvis:scan"))
    assert _kill_order_key(process) < _kill_order_key(JarvisProcess(pid=2, command="Jarvis-Redis:scan"))


def test_sampler_continues_after_checkpoint_process_crashes(tmp_path):
    from test_sample_checkpoint import config, make_sampler, ReverseQueue

    sampler = make_sampler(config(tmp_path, "Random"))
    sampler.configure_checkpoint(checkpoint_file=str(tmp_path / "runtime.pickle"))
    writer = sampler._pickle_writer
    try:
        assert writer.wait()
        assert writer.submit({"probe": KilledChildWrite()})
        until(lambda: not writer.is_alive())
        sampler.set_redis(ReverseQueue(20, sampler))
        assert sampler.run_distributed(timeout=2) == 20
        assert sampler._completed_count + sampler._failed_count == 20
        assert not sampler._unfinished
    finally:
        sampler.shutdown_checkpointing()


def test_control_owns_and_stops_the_checkpoint_component(tmp_path):
    from jarvishep2.runtime.core import Jarvis2Core
    from test_sample_checkpoint import config, make_sampler

    core = Jarvis2Core(config(tmp_path, "Random"))
    core.sampler = make_sampler(core.config)
    try:
        core.start_runtime_checkpoint()
        writer = core.checkpoint_process
        assert writer is core.sampler._pickle_writer
        assert writer.is_alive()
        assert core.sampler._owns_checkpoint is False
        core.start_runtime_checkpoint()
        assert core.checkpoint_process is writer
    finally:
        core.shutdown()
    assert core.checkpoint_process is None
    assert not writer.is_alive()


def test_control_crash_does_not_leave_a_live_checkpoint_process(tmp_path):
    script = tmp_path / "crash_control.py"
    script.write_text('''
import logging, os, sys, time
from pathlib import Path
from jarvishep2.runtime.checkpoint_process import CheckpointProcess
if __name__ == "__main__":
    writer = CheckpointProcess(sys.argv[1], logging.getLogger("test"))
    writer.start()
    writer.request_snapshot()
    time.sleep(0.3)
    Path(sys.argv[2]).write_text(str(writer.pid))
    os._exit(0)
''')
    result = subprocess.run([sys.executable, str(script), str(tmp_path / "runtime.pickle"),
                             str(tmp_path / "child.pid")], timeout=10, capture_output=True)
    assert result.returncode == 0, result.stderr.decode()
    pid = int((tmp_path / "child.pid").read_text())
    def stopped():
        try:
            return psutil.Process(pid).status() == psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:
            return True
    try:
        until(stopped)
    finally:
        if not stopped():
            psutil.Process(pid).kill()
