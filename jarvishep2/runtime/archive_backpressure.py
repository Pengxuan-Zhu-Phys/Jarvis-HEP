"""Factory policy for throttling new claims while SAMPLE packing falls behind."""

from __future__ import annotations

import threading
import time
from collections import deque
from dataclasses import dataclass, field
from typing import Mapping


@dataclass
class ArchiveWatermarks:
    high_sec: float
    low_sec: float
    paused: bool = False
    metric: str = "buckets"
    estimate: float | None = None
    reason: str = "below_high_water"
    _owner: str = ""
    _progress: str = ""
    _last_progress_at: float = 0.0
    _last_open_at: float = float("-inf")
    _observations: int = 0
    _rate_window: deque[tuple[float, int, int]] = field(default_factory=deque)

    def _open(self, now: float, reason: str) -> bool:
        if self.paused:
            self._last_open_at = now
        self._observations = 0
        self.paused = False
        self.reason = reason
        return False

    def observe(self, now: float, stats: Mapping[str, str]) -> bool:
        owner = stats.get("pack_owner", "")
        ready = int(stats.get("ready_buckets", 0))
        active = int(stats.get("active_buckets", 0))
        buckets = ready + active
        self.estimate = None
        if owner != self._owner:
            self._owner = owner
            self._progress = ""
            self._rate_window.clear()
            self._open(now, "packer_restarted")
        if not owner or not buckets:
            self._rate_window.clear()
            self._progress = ""
            return self._open(now, "empty" if owner else "packer_unavailable")

        progress = stats.get("pack_progress", "")
        if progress != self._progress:
            self._last_progress_at = now
            self._progress = progress
        if not progress or now - self._last_progress_at >= 60.0:
            self._rate_window.clear()
            return self._open(now, "packer_stalled")

        points = int(stats.get("outstanding_points", 0))
        completed = int(stats.get("completed_points", 0))
        finished = int(stats.get("completed_buckets", 0))
        while self._rate_window and now - self._rate_window[0][0] > 15.0:
            self._rate_window.popleft()
        if self._rate_window:
            old_time, old_points, old_buckets = self._rate_window[0]
            if finished - old_buckets >= 4 and completed > old_points:
                self.estimate = points * (now - old_time) / (completed - old_points)
        self._rate_window.append((now, completed, finished))

        slots = int(stats.get("packer_slots", 1))
        metric = "seconds" if self.estimate is not None else "buckets"
        if self.paused and self.metric == "buckets":
            metric = "buckets"
        if self.metric != metric:
            self.metric = metric
            self._observations = 0
        if metric == "seconds":
            high = self.estimate >= self.high_sec
            low = self.estimate <= self.low_sec
        else:
            high = buckets >= max(8, 4 * slots)
            low = buckets <= max(2, slots)

        if self.paused:
            self._observations = self._observations + 1 if low else 0
            if self._observations >= 3:
                return self._open(now, "low_water")
        elif now - self._last_open_at >= 2.0:
            self._observations = self._observations + 1 if high else 0
            if self._observations >= 2:
                self.paused = True
                self.reason = "high_water"
                self._observations = 0
        return self.paused


class ArchiveAdmissionLoop:
    """Owns the renewable performance gate for a single Factory run."""

    def __init__(self, redis, logger, *, high_water_sec: float, low_water_sec: float):
        self.redis = redis
        self.logger = logger
        self.policy = ArchiveWatermarks(high_water_sec, low_water_sec)
        self._stop = threading.Event()
        self._thread = threading.Thread(
            target=self._run, name="Jarvis-ArchiveAdmission", daemon=True
        )

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> None:
        if self._stop.is_set():
            return
        self._stop.set()
        # Redis serializes stop against a concurrent lease renewal. Once
        # stopped, neither a renewal nor an automatic open can take effect.
        try:
            self.redis.set_archive_admission("stop")
        except Exception as exc:
            self.logger.warning("could not stop archive admission gate -> %s", exc)

    def join(self) -> None:
        if self._thread.ident is not None:
            self._thread.join(timeout=2.0)

    def _run(self) -> None:
        paused = False
        started = 0.0
        error_logged = False
        while not self._stop.is_set():
            now = time.monotonic()
            try:
                stats = self.redis.sample_pack_pressure()
                desired = self.policy.observe(now, stats)
                if desired or paused:
                    if not self.redis.set_archive_admission("close" if desired else "open"):
                        return
                if desired != paused:
                    if desired:
                        started = now
                    self.logger.info(
                        "archive admission %s: reason=%s metric=%s ready=%s active=%s "
                        "failed=%s estimated_drain_sec=%s revision=%d pause_sec=%.3f",
                        "paused" if desired else "open", self.policy.reason, self.policy.metric,
                        stats.get("ready_buckets", "0"), stats.get("active_buckets", "0"),
                        stats.get("failed_buckets", "0"),
                        round(self.policy.estimate, 2) if self.policy.estimate is not None else "unknown",
                        self.redis._dispatch_revision,
                        0.0 if desired else now - started,
                    )
                paused = desired
                self.redis.publish_proc_board(
                    "core",
                    archive_admission="paused" if paused else "open",
                    archive_admission_reason=self.policy.reason,
                    archive_admission_revision=self.redis._dispatch_revision,
                    archive_ready_buckets=stats.get("ready_buckets", "0"),
                    archive_active_buckets=stats.get("active_buckets", "0"),
                    archive_failed_buckets=stats.get("failed_buckets", "0"),
                    archive_estimated_drain_sec=(
                        self.policy.estimate if self.policy.estimate is not None else "unknown"
                    ),
                )
                error_logged = False
            except Exception as exc:
                # No valid observation means no renewal. Redis expires the
                # performance pause independently of this thread.
                self.policy._open(now, "observation_unavailable")
                if not error_logged:
                    self.logger.warning("archive admission observation unavailable -> %s", exc)
                    error_logged = True
            self._stop.wait(0.5)
