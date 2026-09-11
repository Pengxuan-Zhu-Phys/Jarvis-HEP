"""In-process fake scan used by the hidden splash ``S`` binding."""

from __future__ import annotations

import math
import time
from dataclasses import dataclass

from jarvishep2.monitor.scans import ScanChoice, simulated_choice

_SPARK_SAMPLES = (
    3, 4, 5, 6, 8, 10, 9, 11, 12, 14, 13, 12, 15, 16, 14, 13,
    15, 17, 18, 16, 15, 14, 16, 18, 19, 17, 16, 18, 20, 19, 18, 17,
)
_SPARK_QUEUE = (
    8, 9, 12, 14, 18, 22, 20, 16, 14, 11, 9, 10, 13, 17, 24, 22,
    19, 15, 12, 10, 8, 9, 12, 16, 20, 24, 21, 18, 14, 12, 10, 11,
)
_SPARK_CPU = (
    22, 28, 35, 40, 44, 38, 33, 41, 48, 52, 47, 43, 39, 36, 42, 49,
    55, 51, 46, 41, 37, 34, 39, 45, 50, 47, 41, 38, 36, 40, 44, 41,
)


@dataclass(frozen=True)
class OverviewFrame:
    scan: str
    mode: str
    elapsed: str
    ref: str
    redis: str
    hz: str
    target: int
    done: int
    running: int
    failed: int
    rate: str
    eta: str
    avg: str
    task_q: int
    archive_q: int
    feedback_q: int
    workers_alive: int
    workers_total: int
    stale: int
    busy: int
    cpu: float
    mem_g: float
    mem_total_g: float
    method: str
    run_id: str
    host: str
    started: str
    calculators: tuple[tuple[str, int, int], ...]
    spark_samples: tuple[int, ...]
    spark_queue: tuple[int, ...]
    spark_cpu: tuple[int, ...]
    fds: int
    fds_limit: int


class SimuEngine:
    """Mutating fake scan. ``tick()`` yields a new snapshot each interval."""

    def __init__(self) -> None:
        self._t0 = time.monotonic()
        self.done = 1842
        self.failed = 7
        self.running = 12
        self.target = 3000
        self.task_q = 24
        self.archive_q = 6
        self.workers_alive = 47
        self.workers_total = 48
        self.stale = 1
        self.busy = 12
        self.cpu = 41.0
        self.mem_g = 18.4
        self.fds = 1240
        self.fds_limit = 10240
        self._phase = 0
        self._samples = list(_SPARK_SAMPLES)
        self._queue = list(_SPARK_QUEUE)
        self._cpu = list(_SPARK_CPU)
        self._calc_busy = [4, 0, 2, 1, 0, 2, 0, 1]
        self._calc_total = [16, 8, 12, 12, 16, 10, 4, 7]
        self._calc_names = (
            "SoftSUSY",
            "micrOMEGAs",
            "HiggsBounds",
            "HiggsSignals",
            "SModelS",
            "HiggsTools",
            "FeynHiggs",
            "SPheno",
        )

    def tick(self) -> OverviewFrame:
        self._phase += 1
        elapsed_s = max(0, int(time.monotonic() - self._t0))
        wave = math.sin(self._phase / 4.0)
        if self.done + self.running < self.target and self._phase % 2 == 0:
            self.done += 1
            if self._phase % 17 == 0:
                self.failed += 1
                self.done = max(0, self.done - 1)
        self.running = 8 + (self._phase % 7)
        self.busy = self.running
        self.task_q = max(0, int(18 + 8 * wave))
        self.archive_q = max(0, int(4 + 4 * math.sin(self._phase / 5.0)))
        self.cpu = max(8.0, min(92.0, 41.0 + 18.0 * math.sin(self._phase / 6.0)))
        self.mem_g = max(12.0, min(40.0, 18.4 + 3.0 * math.sin(self._phase / 9.0)))
        self.fds = int(max(400, min(8000, 1240 + 280 * math.sin(self._phase / 7.0))))
        self.stale = 1 if self._phase % 11 else 0
        self.workers_alive = self.workers_total - self.stale
        for index, total in enumerate(self._calc_total):
            wander = int((self._phase + index * 3) % (total // 2 + 1))
            self._calc_busy[index] = max(0, min(total, wander))
        rate = max(4.0, 12.4 + 2.4 * math.sin(self._phase / 5.0))
        remain = max(0, self.target - self.done - self.running)
        eta_s = int(remain / max(rate / 60.0, 0.05))
        self._samples.append(int(rate))
        self._queue.append(self.task_q)
        self._cpu.append(int(self.cpu))
        self._samples = self._samples[-64:]
        self._queue = self._queue[-64:]
        self._cpu = self._cpu[-64:]
        hours, rem = divmod(elapsed_s, 3600)
        minutes, seconds = divmod(rem, 60)
        eta_h, eta_r = divmod(eta_s, 3600)
        eta_m, eta_sec = divmod(eta_r, 60)
        calcs = tuple(
            (name, busy, total)
            for name, busy, total in zip(
                self._calc_names, self._calc_busy, self._calc_total
            )
        )
        return OverviewFrame(
            scan="simu-iDM_Vector_V1",
            mode="running",
            elapsed=f"{hours:02d}:{minutes:02d}:{seconds:02d}",
            ref="SIM",
            redis="127.0.0.1:6379",
            hz="2.0 Hz",
            target=self.target,
            done=self.done,
            running=self.running,
            failed=self.failed,
            rate=f"{rate:0.1f} /min",
            eta=f"{eta_h}:{eta_m:02d}:{eta_sec:02d}",
            avg="3.82 s",
            task_q=self.task_q,
            archive_q=self.archive_q,
            feedback_q=0,
            workers_alive=self.workers_alive,
            workers_total=self.workers_total,
            stale=self.stale,
            busy=self.busy,
            cpu=self.cpu,
            mem_g=self.mem_g,
            mem_total_g=64.0,
            method="AdaptiveBridson",
            run_id="run-20260910-120104",
            host="nersc-login-02",
            started="2026-09-10 12:01:04 UTC",
            calculators=calcs,
            spark_samples=tuple(self._samples),
            spark_queue=tuple(self._queue),
            spark_cpu=tuple(self._cpu),
            fds=self.fds,
            fds_limit=self.fds_limit,
        )


def simu_choice() -> ScanChoice:
    return simulated_choice()


def simu_overview() -> OverviewFrame:
    return SimuEngine().tick()


__all__ = [
    "OverviewFrame",
    "SimuEngine",
    "simu_choice",
    "simu_overview",
]
