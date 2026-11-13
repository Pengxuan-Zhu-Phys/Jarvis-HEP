"""Archive admission invariants, including concurrent claims on real Redis."""
from __future__ import annotations

import shutil
import subprocess
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from jarvishep2.queue._redis_dispatch import DISPATCH_WAIT, DISPATCH_WAKE
from jarvishep2.queue.redis_queue import RedisQueue
from jarvishep2.runtime.archive_backpressure import ArchiveWatermarks
from jarvishep2.runtime.dispatch_wait import DispatchWait
from jarvishep2.runtime_config import get_archive_backpressure_config


def test_archive_admission_policy(subtests):
    for overrides, enabled in [
        ({}, True),
        ({"store_samples": False}, False),
        ({"archiver": {"pack_buckets": False}}, False),
        ({"sample_directory": {"enabled": False}}, False),
        ({"sample_directory": {"pack": False}}, False),
    ]:
        with subtests.test(overrides=overrides):
            v2 = {"factory": {"archive_backpressure": {"enabled": True}}, **overrides}
            card = {"EnvReqs": {"V2": v2}}
            assert get_archive_backpressure_config(card)["enabled"] is enabled
            card["_check_modules_sample_layout"] = True
            assert not get_archive_backpressure_config(card)["enabled"]

    policy = ArchiveWatermarks(10, 3)
    stats = {"pack_owner": "a", "pack_progress": "1", "ready_buckets": "8", "packer_slots": "2"}
    assert not policy.observe(0, stats)
    assert policy.observe(.5, stats)
    stats.update(ready_buckets="2")
    assert policy.observe(1, stats)
    assert policy.observe(1.5, stats)
    assert not policy.observe(2, stats)
    stats.update(ready_buckets="8")
    assert not policy.observe(2.5, stats)
    assert not policy.observe(4, stats)
    assert policy.observe(4.5, stats)
    assert not policy.observe(60.5, stats)
    assert policy.reason == "packer_stalled"
    stats.update(pack_owner="b", pack_progress="2", ready_buckets="0")
    assert not policy.observe(61, stats)

    policy = ArchiveWatermarks(10, 3)
    stats = {
        "pack_owner": "rate", "pack_progress": "1", "ready_buckets": "6",
        "outstanding_points": "1800", "completed_buckets": "0", "completed_points": "0",
    }
    assert not policy.observe(0, stats)
    stats.update(completed_buckets="4", completed_points="120", pack_progress="2")
    assert not policy.observe(1, stats)
    assert policy.observe(1.5, stats)
    assert policy.metric == "seconds"
    stats.update(outstanding_points="120")
    assert policy.observe(2, stats)
    assert policy.observe(2.5, stats)
    assert not policy.observe(3, stats)


@pytest.fixture
def admission_redis():
    import redis

    exe = shutil.which("redis-server")
    if exe is None:
        pytest.skip("requires real Redis")
    # macOS Unix sockets have a short path limit.
    with tempfile.TemporaryDirectory(prefix="jhep-admit-", dir="/tmp") as directory:
        socket = str(Path(directory) / "r.sock")
        server = subprocess.Popen(
            [exe, "--port", "0", "--unixsocket", socket, "--save", "", "--appendonly", "no"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        try:
            client = redis.Redis(unix_socket_path=socket, decode_responses=True)
            deadline = time.monotonic() + 5
            while True:
                try:
                    client.ping()
                    break
                except redis.ConnectionError:
                    if time.monotonic() >= deadline:
                        raise
                    time.sleep(.02)
            yield socket
        finally:
            server.terminate()
            server.communicate(timeout=5)


@pytest.mark.slow
@pytest.mark.parametrize("codec", ["json", "msgpack"])
def test_atomic_admission_ownership_wake_and_per_task_wait(admission_redis, codec):
    import redis

    client = redis.Redis(unix_socket_path=admission_redis, decode_responses=codec == "json")
    queue = RedisQueue({"codec": codec}, client=client)
    queue.init_archive_admission("run-a", 4)
    queue.push_task({"uuid": "running", "u_coords": [.1]})
    assert queue.pull_task_to_inflight("0")["uuid"] == "running"
    assert queue.set_archive_admission("close", lease_sec=2)
    queue.push_many_tasks([
        {"uuid": f"queued-{i}", "u_coords": [.2], "context": b"opaque" if codec == "msgpack" else ""}
        for i in range(8)
    ])
    queued_wait = DispatchWait(queue, .10, ["queued-0"])
    mixed_wait = DispatchWait(queue, .10, ["queued-0", "running"])
    with ThreadPoolExecutor(max_workers=4) as executor:
        assert list(executor.map(queue.pull_task_to_inflight, map(str, range(1, 5)))) == [None] * 4
    assert queue.pull_task_to_inflight("0")["uuid"] == "running"
    assert queue.fetch_sample_stats()["running"] == 1
    assert client.llen(DISPATCH_WAKE) == 0
    time.sleep(.15)
    assert queued_wait.remaining(["queued-0"]) > .05
    assert mixed_wait.remaining(["queued-0", "running"]) < 0

    # Closing again renews the same pause. Opening wakes all available Workers.
    assert queue.set_archive_admission("close", lease_sec=2)
    assert queue.set_archive_admission("open")
    assert client.llen(DISPATCH_WAKE) == 4
    ready = threading.Barrier(4)
    def claim(worker):
        ready.wait(timeout=5)
        return queue.pull_task_to_inflight(str(worker))
    with ThreadPoolExecutor(max_workers=4) as executor:
        claimed = list(executor.map(claim, range(1, 5)))
    assert {task["uuid"] for task in claimed} == {f"queued-{i}" for i in range(4)}

    # A retry retains the first attempt's queue credit and gains credit only
    # after its own requeue. A stopped owner cannot be reopened.
    reclaimed = queue.reclaim_inflight_task("0")
    assert len(reclaimed) == 1
    queue.set_archive_admission("close", lease_sec=.12)
    queue.lpush_tasks(reclaimed)
    assert queue.pull_task_to_inflight("5") is None
    time.sleep(.16)
    assert queue.pull_task_to_inflight("5")["uuid"] == "running"
    assert .07 < queue.dispatch_wait_credits(["running"])["running"] < .14
    queue.publish_feedback({"uuid": "running", "logL": 0.})
    assert queue.ack_inflight_task("5", "running")
    assert not queue.ack_inflight_task("5", "running")
    assert queue.pull_feedback()["uuid"] == "running"
    assert not client.hexists(DISPATCH_WAIT, "running")
    queue.set_archive_admission("stop")
    assert not queue.set_archive_admission("open")
    assert queue.pull_task_to_inflight("6") is None
    queue.clear_archive_admission()
    assert not client.exists(DISPATCH_WAIT)
    queue.reconcile_resume_ephemeral(completed=0)
    queue.init_archive_admission("run-b", 4)
    old = RedisQueue({"codec": codec}, client=client)
    old._dispatch_epoch = "run-a"
    old._dispatch_has_paused = True
    assert old.dispatch_wait_credits(["running"]) == {}
    assert not old.set_archive_admission("close")
    assert queue.pull_task_to_inflight("6") is None
    woke = threading.Event()
    def wait_for_new_work():
        queue._wait_for_admission()
        woke.set()
    waiter = threading.Thread(target=wait_for_new_work)
    waiter.start()
    queue.push_task({"uuid": "new", "u_coords": [.3]})
    assert woke.wait(.5)
    waiter.join()
    assert queue.pull_task_to_inflight("6")["uuid"] == "new"


@pytest.mark.slow
def test_bucket_pressure_terminal_and_packer_replacement(admission_redis, tmp_path):
    import redis

    client = redis.Redis(unix_socket_path=admission_redis, decode_responses=True)
    queue = RedisQueue(client=client)
    queue.init_sample_buckets(sample_root=str(tmp_path), limit=1, observe=True)
    queue._pack_transition("start", "old", arg=2)
    for bucket in range(1, 4):
        assert queue.allocate_sample_bucket()["bucket_id"] == bucket
        queue.finish_sample_bucket(bucket)
        queue.note_sample_archived(bucket)
    stats = queue.sample_pack_pressure()
    assert stats["ready_buckets"] == "2"
    assert stats["outstanding_points"] == "2"
    first = queue.pull_ready_bucket(owner="old")
    assert first["bucket_id"] == 1
    assert queue.sample_pack_pressure()["active_buckets"] == "1"
    queue._pack_transition("finish", "old", "1", 1)
    queue._pack_transition("finish", "old", "1", 1)
    assert queue.sample_pack_pressure()["completed_buckets"] == "1"
    assert queue.pull_ready_bucket(owner="old")["bucket_id"] == 2
    queue._pack_transition("start", "new", arg=2)
    assert not queue._pack_transition("finish", "old", "2", 1)
    stats = queue.sample_pack_pressure()
    assert (stats["active_buckets"], stats["outstanding_points"], stats["failed_buckets"]) == ("0", "0", "1")
    assert queue.pull_ready_bucket(owner="new") is None


@pytest.mark.slow
@pytest.mark.skipif(not shutil.which("redis-server"), reason="requires real Redis")
@pytest.mark.parametrize("method", [
    "Grid", "Random", "Bridson", "CSV", "AdaptiveBridson", "Dynesty", "MultiNest",
    "ToyMCMC", "MCMC", "AMMCMC", "DRAM", "DEMCMC", "EnsembleMCMC", "PTMCMC", "PTEnsemble",
    "RecoverIndexed",
])
def test_spawn_sampler_survives_closed_admission(tmp_path, method):
    """Exercise each sampler's actual submission/feedback path with a closed gate."""
    import json
    import os
    import sys

    import yaml
    from test_resume_protocol import _card

    finite = {"Grid", "Random", "Bridson", "CSV"}
    source, card = _card(tmp_path, method if method in finite else "Random")
    card["Sampling"]["Method"] = method
    plugin = None
    if method in {"Dynesty", "MultiNest"}:
        options = {"maxiter": 12, "print_progress": False}
        if method == "Dynesty":
            options["maxbatch"] = 1
        card["Sampling"]["Bounds"] = {"nlive": 8, "seed": 11, "run_nested": options}
    elif method == "AdaptiveBridson":
        card["Sampling"]["Bounds"] = {
            "seed": 11, "target_expression": "z", "target_value": 10.0,
            "initial_radius": .12, "max_generations": 2, "max_points": 32,
            "max_new_per_generation": 16,
        }
    elif method == "RecoverIndexed":
        from test_plugin_resume_runtime import SOURCE

        plugin = tmp_path / "RecoverIndexed.py"
        plugin.write_text(SOURCE[:SOURCE.index("def checkpoint")] + SOURCE[SOURCE.index("def draw_block"):])
    elif method not in finite:
        card["Sampling"]["Bounds"] = {
            "seed": 11, "num_chains": 4, "num_iters": 3, "proposal_scale": .1,
        }
        if method in {"PTMCMC", "PTEnsemble"}:
            card["Sampling"]["Bounds"]["temperature_ladder"] = [1.0, 2.0, 4.0, 8.0]
    card["EnvReqs"]["V2"].update(
        factory={
            "watchdog": {"enabled": method == "Random", "poll_interval_sec": .1},
            "archive_backpressure": {"enabled": True},
        },
        archiver={"mode": "thread" if method == "CSV" else "process"},
    )
    source.write_text(yaml.safe_dump(card))
    driver = """
import json, os, signal, sys, threading, time
from pathlib import Path
from jarvishep2.queue.redis_queue import RedisQueue
from jarvishep2.runtime.core import Jarvis2Core

enqueue = RedisQueue._enqueue_with_admission
evidence = []
timer = None
def gated(self, direction, tasks, encoded):
    global timer
    if timer is not None:
        return enqueue(self, direction, tasks, encoded)
    replace_worker = core.config["Sampling"]["Method"] == "Random"
    assert self.set_archive_admission("close", lease_sec=15 if replace_worker else 3)
    submitted = enqueue(self, direction, tasks, encoded)
    uuid = str(tasks[0]["uuid"])
    if replace_worker:
        os.kill(core.factory.workers[0].pid, signal.SIGKILL)
    def release():
        if replace_worker:
            deadline = time.monotonic() + 8
            while not core.factory._respawn_count and time.monotonic() < deadline:
                time.sleep(.02)
        evidence.append((self.dispatch_wait_credits([uuid])[uuid],
                         self.fetch_sample_stats()["running"], core.factory._respawn_count))
        assert self.set_archive_admission("open")
    timer = threading.Timer(.35, release)
    timer.start()
    return submitted

RedisQueue._enqueue_with_admission = gated
if len(sys.argv) == 3:
    from jarvis_hep_sampler import SamplerRegistry
    SamplerRegistry().load(Path(sys.argv[2]))
core = Jarvis2Core()
core.load_task_yaml(sys.argv[1])
core.config.update(skip_draw_flowchart=True, skip_emit_plot_scenes=True, skip_render_plots=True)
outcome = core.run(write_run_summary=False)
assert timer is not None
timer.join(2)
Path("outcome.json").write_text(json.dumps([outcome.to_dict(), evidence]))
if len(sys.argv) == 3:
    assert not list(Path("checkpoints").rglob("runtime.pickle"))
raise SystemExit(outcome.exit_code)
"""
    args = [str(source)]
    if plugin is not None:
        args.append(str(plugin))
    result = subprocess.run(
        [sys.executable, "-c", driver, *args], cwd=tmp_path, text=True,
        capture_output=True, timeout=90,
        env=dict(os.environ, JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "refs.json"),
                 JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"),
                 JHS_TEST_RESTORE=str(tmp_path / "prepare.txt"),
                 JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curve-cache")),
    )
    assert result.returncode == 0, (result.stdout + result.stderr)[-8000:]
    outcome, evidence = json.loads((tmp_path / "outcome.json").read_text())
    assert evidence and evidence[0][0] >= .25 and evidence[0][1] == 0
    if method == "Random":
        assert evidence[0][2] == 1
    assert outcome["archived"] == outcome["submitted"] > 0


@pytest.mark.slow
@pytest.mark.skipif(not shutil.which("redis-server"), reason="requires real Redis")
def test_interrupt_during_gate_pause_resumes_queued_checkpoint(tmp_path):
    import os
    import sys

    from jarvishep2.database import SimpleHDF5Writer
    from jarvishep2.sampling.sample_checkpoint import read_pickle
    from test_resume_protocol import _card

    source, card = _card(tmp_path, "Random", count=16)
    env = dict(os.environ, JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "refs.json"))
    driver = """
import sys
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.sampling.completion_sampler import CompletionSampler

core = Jarvis2Core()
core.load_task_yaml(sys.argv[1])
core.config.update(skip_draw_flowchart=True, skip_emit_plot_scenes=True, skip_render_plots=True)
resume = sys.argv[2] == 'resume'
if not resume:
    original = CompletionSampler.flush_batch
    def interrupt(self, batch):
        if not batch:
            return original(self, batch)
        assert self.redis.set_archive_admission('close', lease_sec=30)
        original(self, batch)
        assert self.redis.fetch_sample_stats()['running'] == 0
        assert self.save_checkpoint(force=True)
        assert self._pickle_writer.wait()
        raise KeyboardInterrupt
    CompletionSampler.flush_batch = interrupt
outcome = core.run(resume=resume, write_run_summary=False)
raise SystemExit(outcome.exit_code)
"""
    name = card["Scan"]["name"]
    for phase, code in (("interrupt", 130), ("resume", 0)):
        result = subprocess.run(
            [sys.executable, "-c", driver, str(source), phase], cwd=tmp_path,
            env=env, text=True, capture_output=True, timeout=45,
        )
        assert result.returncode == code, (result.stdout + result.stderr)[-8000:]
        if phase == "interrupt":
            checkpoint = tmp_path / f"checkpoints/{name}/Random/runtime.pickle"
            state = read_pickle(str(checkpoint))["sampler_state"]
            assert len(state["unfinished"]) == 1
            queued = next(iter(state["unfinished"].values()))["uuid"]
    rows = SimpleHDF5Writer(str(tmp_path / f"outputs/{name}/DATABASE/samples.hdf5")).read_records()
    assert len(rows) == 16
    assert sum(row["uuid"] == queued for row in rows) == 1
