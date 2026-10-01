"""Monitor resource figures: this scan's process tree, read from the OS by PID."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
import unittest
from dataclasses import replace

from jarvishep2.monitor.collector import Collector, _CpuMeter
from jarvishep2.monitor.live_overview import LiveOverviewProjector
from jarvishep2.monitor.collector import CollectorFrame
from jarvishep2.monitor.overview import render_resource_subtitle
from jarvishep2.monitor.scans import ScanChoice
from jarvishep2.monitor.simu import SimuEngine
from jarvishep2.redis_queue import make_fakeredis_queue

# A "worker" that keeps starting short CPU-bound children, like fast
# calculators that begin and end between two monitor refreshes.
_TREE = r"""
import subprocess, sys, time
if len(sys.argv) > 1:
    end = time.time() + 6
    while time.time() < end:
        subprocess.run([sys.executable, "-c",
                        "import time\nt=time.time()\nwhile time.time()-t<0.2: pass"])
else:
    subprocess.Popen([sys.executable, __file__, "worker"]).wait()
"""


class CpuMeterTests(unittest.TestCase):
    def test_first_reading_is_unknown_then_cores_per_second(self) -> None:
        meter = _CpuMeter()
        self.assertIsNone(meter.rate("scan", 10.0, now=100.0))
        self.assertAlmostEqual(meter.rate("scan", 13.0, now=102.0), 1.5)

    def test_unmeasurable_readings_repeat_the_last_value(self) -> None:
        meter = _CpuMeter()
        meter.rate("scan", 10.0, now=100.0)
        self.assertAlmostEqual(meter.rate("scan", 12.0, now=101.0), 2.0)
        # Refresh that came too soon: keep the previous baseline and value.
        self.assertAlmostEqual(meter.rate("scan", 12.01, now=101.01), 2.0)
        # A process left the group, so the total went down.
        self.assertAlmostEqual(meter.rate("scan", 5.0, now=102.0), 2.0)
        # Measuring resumes from the new baseline.
        self.assertAlmostEqual(meter.rate("scan", 6.0, now=103.0), 1.0)


class ScanTreeCpuTests(unittest.TestCase):
    def test_short_lived_children_are_counted(self) -> None:
        with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as handle:
            handle.write(_TREE)
        self.addCleanup(os.unlink, handle.name)
        root = subprocess.Popen([sys.executable, handle.name])
        try:
            time.sleep(0.8)
            collector = Collector(
                make_fakeredis_queue(),
                process_inventory=[{"role": "core", "pid": root.pid}],
            )
            first = collector._collect_host_snapshot()
            self.assertIsNone(first["cpu_cores_used"])
            self.assertIsNone(first["cpu_percent"])
            readings = []
            for _ in range(3):
                time.sleep(0.6)
                readings.append(collector._collect_host_snapshot())
        finally:
            import psutil

            for child in psutil.Process(root.pid).children(recursive=True):
                child.kill()
            root.kill()
            root.wait()
        cores = [reading["cpu_cores_used"] for reading in readings]
        # Each child lives 0.2 s, so sampling only live processes reads ~0.
        # Cumulative CPU time including reaped children sees about one core.
        self.assertGreater(min(cores), 0.3, cores)
        self.assertLess(max(cores), 1.6, cores)
        last = readings[-1]
        self.assertAlmostEqual(
            last["cpu_percent"],
            100.0 * last["cpu_cores_used"] / last["cpu_cores_available"],
        )
        self.assertEqual(last["processes"][0]["role"], "core")
        self.assertGreater(last["memory_used"], 0)


class ResourceDisplayTests(unittest.TestCase):
    def _project(self, host: dict) -> object:
        choice = ScanChoice(
            reference="R1",
            name="scan",
            control_pid=1,
            process_count=1,
            pids=(1,),
            simulated=False,
        )
        return LiveOverviewProjector(choice).project(CollectorFrame(page="overview", host=host))

    def test_first_refresh_shows_memory_and_unknown_cpu(self) -> None:
        frame = self._project(
            {
                "available": True,
                "cpu_percent": None,
                "cpu_cores_used": None,
                "cpu_cores_available": 8,
                "memory_used": 2 * 1024**3,
                "memory_total": 16 * 1024**3,
                "processes": [],
            }
        )
        self.assertTrue(frame.resources_available)
        self.assertIsNone(frame.cpu)
        self.assertEqual(render_resource_subtitle(frame), "— cores · fds —")

    def test_cores_and_fd_closest_to_its_limit(self) -> None:
        frame = self._project(
            {
                "available": True,
                "cpu_percent": 42.5,
                "cpu_cores_used": 3.4,
                "cpu_cores_available": 8,
                "memory_used": 2 * 1024**3,
                "memory_total": 16 * 1024**3,
                "processes": [
                    {"fds": 500, "fd_limit": 4096},
                    {"fds": 900, "fd_limit": 1024},
                    {"fds": 10, "fd_limit": None},
                ],
            }
        )
        self.assertEqual(frame.cpu, 42.5)
        self.assertEqual(
            render_resource_subtitle(frame), "3.4 / 8 cores · fds 900/1024"
        )

    def test_simulated_frame_carries_core_counts(self) -> None:
        frame = SimuEngine().tick((10, 10))
        self.assertEqual(frame.cpu_cores_available, 64)
        self.assertIn("/ 64 cores", render_resource_subtitle(replace(frame)))


if __name__ == "__main__":
    unittest.main()
