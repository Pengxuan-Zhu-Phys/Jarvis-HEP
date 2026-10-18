"""Host and Redis failures must not invalidate the other monitor source."""

from unittest.mock import Mock, patch

import pytest

from jarvishep2.monitor.collector import Collector
from jarvishep2.monitor.live_overview import LiveOverviewProjector
from jarvishep2.monitor.scans import ScanChoice
from jarvishep2.monitor.status_pages import HostPane
from jarvishep2.redis_queue import SAMPLE_STATS, make_fakeredis_queue


@pytest.fixture
def sources():
    queue = make_fakeredis_queue()
    queue.r.hset(SAMPLE_STATS, mapping={"completed": 9})
    host = Mock(return_value={
        "available": True, "cpu_percent": None,
        "memory_used": 1024**3, "memory_total": 8 * 1024**3,
        "processes": [],
    })
    return queue, host, Collector(queue, host_snapshotter=host)


def projector():
    return LiveOverviewProjector(ScanChoice(
        reference="R1", name="scan", control_pid=1,
        process_count=1, pids=(1,), simulated=False,
    ))


def host_text(frame):
    pane = HostPane()
    left, right = Mock(), Mock()
    with patch.object(pane, "query_one", side_effect=[left, right]):
        pane.set_frame(frame)
    return left.update.call_args.args[0]


@pytest.mark.parametrize("warm", [False, True])
def test_host_failures_keep_advancing_redis_then_recover(sources, warm, caplog):
    queue, host, collector = sources
    if warm:
        collector.tick("overview")
    host.side_effect = PermissionError("host counters denied")
    view = projector()
    for completed in (10, 11):
        queue.r.hset(SAMPLE_STATS, "completed", completed)
        frame = collector.tick("overview")
        assert frame.sample_stats["completed"] == completed
        assert not frame.stale
        assert frame.error is None
        assert frame.host_stale
        assert frame.host_error == "host counters denied"
        assert frame.host == {
            "available": False, "reason": "host counters denied", "processes": [],
        }
        assert collector.last_frame == frame
        projected = view.project(frame)
        assert projected.done == completed
        assert not projected.resources_available
        assert projected.mode != "stale"
        assert "host counters denied" in host_text(frame)
    assert caplog.text.count("Monitor Host collection failed") == 1
    host.side_effect = None
    recovered = collector.tick("overview")
    assert not recovered.host_stale
    assert recovered.host_error is None
    assert projector().project(recovered).resources_available
    host.side_effect = PermissionError("host counters denied")
    collector.tick("overview")
    assert caplog.text.count("Monitor Host collection failed") == 2


@pytest.mark.parametrize("warm", [False, True])
def test_redis_failure_keeps_refreshing_host_then_recovers(sources, warm, caplog):
    queue, host, collector = sources
    good = collector.tick("overview") if warm else collector.last_frame
    with patch("jarvishep2.monitor.collector.attach_reader", side_effect=ConnectionError("Redis offline")):
        for memory in (2, 3):
            host.return_value["memory_used"] = memory * 1024**3
            frame = collector.tick("overview")
            assert frame.stale
            assert frame.error == "Redis offline"
            assert frame.sample_stats == good.sample_stats
            assert collector.last_frame == good
            assert frame.host["memory_used"] == memory * 1024**3
            assert not frame.host_stale
            assert frame.host_error is None
            assert "[STALE]" not in host_text(frame)
            assert projector().project(frame).resources_available
    assert caplog.text.count("Monitor Redis collection failed") == 1
    queue.r.hset(SAMPLE_STATS, "completed", 12)
    frame = collector.tick("overview")
    assert not frame.stale
    assert frame.error is None
    assert frame.sample_stats["completed"] == 12


@pytest.mark.parametrize("first_recovery", ["host", "redis"])
def test_both_sources_fail_and_recover_independently(sources, first_recovery):
    queue, host, collector = sources
    collector.tick("overview")
    host.side_effect = OSError("host unavailable")
    with patch("jarvishep2.monitor.collector.attach_reader") as reader:
        reader.side_effect = ConnectionError("Redis offline")
        failed = collector.tick("overview")
        assert failed.stale and failed.host_stale
        assert failed.sample_stats["completed"] == 9
        assert failed.error == "Redis offline"
        assert failed.host_error == "host unavailable"
        if first_recovery == "host":
            host.side_effect = None
            frame = collector.tick("overview")
            assert frame.stale and not frame.host_stale
    if first_recovery == "redis":
        queue.r.hset(SAMPLE_STATS, "completed", 14)
        frame = collector.tick("overview")
        assert not frame.stale and frame.host_stale
        assert frame.sample_stats["completed"] == 14
    host.side_effect = None
    recovered = collector.tick("overview")
    assert not recovered.stale and not recovered.host_stale


def test_failed_want_clear_does_not_block_host_page(sources):
    queue, host, collector = sources
    collector.tick("samples")
    with patch.object(queue, "clear_monitor_want", side_effect=ConnectionError("want clear failed")):
        frame = collector.tick("host")
    host.assert_called_once()
    assert frame.stale
    assert frame.error == "want clear failed"
    assert frame.host["available"]
    assert not frame.host_stale
    assert "[STALE]" not in host_text(frame)


def test_unavailable_host_result_and_partial_metrics_are_not_zero(sources):
    queue, host, collector = sources
    host.return_value = {"available": False, "reason": "psutil is not installed"}
    frame = collector.tick("overview")
    assert not frame.stale
    assert frame.host_stale
    assert frame.host_error == "psutil is not installed"
    assert frame.sample_stats["completed"] == 9
    host.return_value = {"available": True, "memory_total": 8 * 1024**3}
    frame = collector.tick("overview")
    assert not frame.host_stale
    assert not projector().project(frame).resources_available
    text = host_text(frame)
    assert "CPU  —" in text
    assert "MEM  ≈—G" in text
    assert "SWAP  —/—G" in text
    assert "0.0" not in text


@pytest.mark.parametrize("page", ["workers", "factory", "sampler", "samples", "calculators"])
def test_other_pages_do_not_request_host_metrics(sources, page):
    _, host, collector = sources
    host.side_effect = AssertionError("unexpected host read")
    frame = collector.tick(page)
    host.assert_not_called()
    assert not frame.stale
    assert not frame.host_stale
    assert frame.host == {}
