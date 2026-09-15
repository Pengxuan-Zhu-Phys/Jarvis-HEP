"""In-process fake scan used by the hidden splash ``S`` binding."""

from __future__ import annotations

import math
import time
from dataclasses import dataclass

from jarvishep2.monitor.calculators import CalculatorPoolData, CalculatorsBlockData
from jarvishep2.monitor.metrics import format_sample_rate
from jarvishep2.monitor.queues import (
    QueueBlockData,
    QueueDirection,
    QueueDirectionWindow,
)
from jarvishep2.monitor.samples import (
    BUILTIN_SAMPLERS,
    SamplerDisplay,
    SamplerProgress,
)
from jarvishep2.monitor.workers import WorkerBlockData
from jarvishep2.monitor.scans import ScanChoice, simulated_choice

SIMU_HZ = 2.0
SIMU_INTERVAL_SEC = 1.0 / SIMU_HZ
QUEUE_TASK_PERIOD_STEPS = 16.0
QUEUE_ARCHIVE_PERIOD_STEPS = 20.0
_SHARDED_FEEDBACK_SAMPLERS = frozenset({"MCMC", "ToyMCMC", "AMMCMC", "DRAM"})
_SHARED_FEEDBACK_SAMPLERS = frozenset(
    {"EnsembleMCMC", "DEMCMC", "PTMCMC", "PTEnsemble"}
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
    queue_block: QueueBlockData
    worker_block: WorkerBlockData
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
    calculator_block: CalculatorsBlockData
    spark_samples: tuple[float | None, ...]
    spark_queue: tuple[float | None, ...]
    fds: int
    fds_limit: int
    sampler_generation: int
    sampler_max_generations: int
    sampler_radius: float
    sampler_core: int
    sampler_open: int
    sampler_state: str
    sampler_saved: int
    sampler: SamplerDisplay
    resources_available: bool = True
    health_items: tuple[tuple[str, str, str, str], ...] = ()


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
        self.feedback_q = 0
        self.workers_alive = 48
        self.workers_total = 48
        self.stale = 0
        self.busy = 12
        self.cpu = 41.0
        self.mem_g = 18.4
        self.fds = 1240
        self.fds_limit = 10240
        self.sampler_generation = 8
        self.sampler_max_generations = 25
        self.sampler_radius = 0.025
        self.sampler_core = 320
        self.sampler_open = 4
        self.sampler_state = "partial"
        self.sampler_saved = 1818
        self._sampler_index = BUILTIN_SAMPLERS.index("AdaptiveBridson")
        self._phase = 0
        # Histories are newest-first.  Once the pane reports its actual column
        # width, the arrays are exactly that wide and start with empty slots.
        self._history_widths = (0, 0)
        self._samples: list[float | None] = []
        self._queue: list[float | None] = []
        self._queue_directions = QueueDirectionWindow()
        self._feedback_shards: dict[str, int] = {}
        self._last_rate_per_sec: float | None = None
        self._last_eta_s: int | None = None
        self._last_elapsed_s: int | None = None
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
        self._calc_active_packs = {
            name: {
                1 + (index * 3 + offset * 5) % total
                for offset in range(self._calc_busy[index])
            }
            for index, (name, total) in enumerate(zip(self._calc_names, self._calc_total))
        }
        self._update_feedback_depths()

    def tick(
        self, history_width: int | tuple[int, int] | None = None
    ) -> OverviewFrame:
        if history_width is not None:
            if isinstance(history_width, tuple):
                self.set_history_widths(*history_width)
            else:
                self.set_history_width(history_width)
        elif self._history_widths == (0, 0):
            # Direct callers without a TUI width still get one real column;
            # the workspace always supplies the measured display width.
            self.set_history_width(1)
        self._phase += 1
        elapsed_s = max(0, int(time.monotonic() - self._t0))
        # One history column represents one display interval: 1 / HZ seconds.
        # The values below are interval means, not instantaneous samples at the
        # right edge of the column.
        task_wave = self._interval_sine_mean(QUEUE_TASK_PERIOD_STEPS)
        archive_wave = self._interval_sine_mean(QUEUE_ARCHIVE_PERIOD_STEPS)
        if self.done + self.running < self.target and self._phase % 2 == 0:
            self.done += 1
            if self._phase % 17 == 0:
                self.failed += 1
                self.done = max(0, self.done - 1)
        self.running = 8 + (self._phase % 7)
        self.busy = self.running
        self.task_q = max(0, int(round(18 + 8 * task_wave)))
        self.archive_q = max(0, int(round(4 + 4 * archive_wave)))
        self._update_feedback_depths()
        self.cpu = max(8.0, min(92.0, 41.0 + 18.0 * self._interval_sine_mean(6.0)))
        self.mem_g = max(12.0, min(40.0, 18.4 + 3.0 * self._interval_sine_mean(9.0)))
        self.fds = int(max(400, min(8000, 1240 + 280 * math.sin(self._phase / 7.0))))
        if self._phase % 24 == 0 and self.sampler_generation < self.sampler_max_generations:
            self.sampler_generation += 1
            self.sampler_radius = max(0.002, self.sampler_radius * 0.82)
        self.sampler_core = 320 + int(round(8 * self._interval_sine_mean(7.0)))
        self.sampler_open = max(0, 4 + int(round(2 * self._interval_sine_mean(4.0))))
        self.sampler_saved = max(self.sampler_saved, max(0, self.done - 24))
        # The normal simulation is healthy; expose a short stale-worker warning
        # periodically so HEALTH also demonstrates its degraded state.
        self.stale = 1 if self._phase % 11 == 0 else 0
        self.workers_alive = self.workers_total - self.stale
        for index, total in enumerate(self._calc_total):
            wander = int((self._phase + index * 3) % (total // 2 + 1))
            self._calc_busy[index] = max(0, min(total, wander))
        self._sync_calculator_activity()
        rate_per_min = max(4.0, 12.4 + 2.4 * self._interval_sine_mean(5.0))
        rate_per_sec = rate_per_min / 60.0
        remain = max(0, self.target - self.done - self.running)
        eta_s = int(remain / max(rate_per_sec, 0.05))
        self._push_history(self._samples, rate_per_sec, self._history_widths[0])
        self._push_history(
            self._queue,
            18 + 8 * task_wave,
            self._history_widths[1],
        )
        self._last_rate_per_sec = rate_per_sec
        self._last_eta_s = eta_s
        self._last_elapsed_s = elapsed_s
        return self._build_frame(
            rate_per_sec=rate_per_sec,
            eta_s=eta_s,
            elapsed_s=elapsed_s,
            queue_observed_at=self._phase * SIMU_INTERVAL_SEC,
        )

    def empty_frame(self) -> OverviewFrame:
        """Return the pre-tick frame: real counters, but no spark columns."""
        return self._build_frame()

    def set_history_width(self, width: int) -> None:
        """Make each history exactly as wide as one rendered SPARKS column."""
        width = max(1, int(width))
        self.set_history_widths(width, width)

    def set_history_widths(self, samples_width: int, queue_width: int) -> None:
        """Set the independent display widths for samples and queue histories."""
        widths = (max(1, int(samples_width)), max(1, int(queue_width)))
        if widths == self._history_widths:
            return
        self._history_widths = widths
        self._samples = self._resize_history(self._samples, widths[0])
        self._queue = self._resize_history(self._queue, widths[1])

    def current_frame(self) -> OverviewFrame:
        """Return the current counters after a width-only history resize."""
        return self._build_frame(
            rate_per_sec=self._last_rate_per_sec,
            eta_s=self._last_eta_s,
            elapsed_s=self._last_elapsed_s,
        )

    @property
    def sampler_method(self) -> str:
        return BUILTIN_SAMPLERS[self._sampler_index]

    def set_sampler_method(self, method: str) -> str:
        """Select one built-in sampler fixture without resetting the scan."""
        self._sampler_index = BUILTIN_SAMPLERS.index(str(method))
        self._queue_directions.reset("feedback")
        self._update_feedback_depths()
        return self.sampler_method

    def cycle_sampler(self, delta: int = 1) -> str:
        """Cycle sampler-specific SAMPLES canvases in the hidden SIMU."""
        self._sampler_index = (
            self._sampler_index + int(delta)
        ) % len(BUILTIN_SAMPLERS)
        self._queue_directions.reset("feedback")
        self._update_feedback_depths()
        return self.sampler_method

    def _interval_sine_mean(self, period_steps: float) -> float:
        """Average a synthetic signal over the current 1 / SIMU_HZ interval."""
        previous = (self._phase - 1) / period_steps
        current = self._phase / period_steps
        # Integral of sin(t / period_steps) over one simulation step.  One
        # step is exactly SIMU_INTERVAL_SEC, so the result is a true bucket
        # mean rather than a right-edge instantaneous sample.
        return period_steps * (math.cos(previous) - math.cos(current))

    def _push_history(
        self, history: list[float | None], value: float, width: int
    ) -> None:
        """Insert a new right-edge value and discard the oldest ring entry."""
        history.insert(0, float(value))
        del history[width:]

    def _feedback_mode(self) -> str:
        if self.sampler_method in _SHARDED_FEEDBACK_SAMPLERS:
            return "sharded"
        if self.sampler_method in _SHARED_FEEDBACK_SAMPLERS:
            return "shared"
        return "inactive"

    def _feedback_chain_count(self) -> int:
        return 4

    def _update_feedback_depths(self) -> None:
        """Produce bounded, method-accurate feedback fixtures for QUEUES."""
        mode = self._feedback_mode()
        if mode == "inactive":
            self.feedback_q = 0
            self._feedback_shards = {}
            return
        if mode == "shared":
            self.feedback_q = max(
                0, int(round(2.0 + 1.5 * math.sin(self._phase / 8.0)))
            )
            self._feedback_shards = {}
            return
        self._feedback_shards = {
            str(index): max(
                0,
                int(round(1.0 + 1.2 * math.sin(self._phase / 6.0 + index))),
            )
            for index in range(self._feedback_chain_count())
        }
        self.feedback_q = sum(self._feedback_shards.values())

    def _queue_block(self, observed_at: float | None) -> QueueBlockData:
        """Build the same bounded hand-off snapshot required by live Monitor."""
        mode = self._feedback_mode()
        if observed_at is None:
            task_direction = QueueDirection()
            archive_direction = QueueDirection()
            feedback_direction = QueueDirection()
        else:
            task_direction = self._queue_directions.observe(
                "task", self.task_q, observed_at=observed_at
            )
            archive_direction = self._queue_directions.observe(
                "archive", self.archive_q, observed_at=observed_at
            )
            feedback_direction = self._queue_directions.observe(
                "feedback",
                None if mode == "inactive" else self.feedback_q,
                observed_at=observed_at,
            )
        return QueueBlockData(
            task_depth=self.task_q,
            archive_depth=self.archive_q,
            feedback_depth=None if mode == "inactive" else self.feedback_q,
            feedback_mode=mode,
            feedback_shards=dict(self._feedback_shards),
            task_direction=task_direction,
            archive_direction=archive_direction,
            feedback_direction=feedback_direction,
        )

    def _sync_calculator_activity(self) -> None:
        """Apply simulated acquires/releases while retaining PackID positions."""
        for index, (name, slots, target) in enumerate(
            zip(self._calc_names, self._calc_total, self._calc_busy)
        ):
            active = self._calc_active_packs.setdefault(name, set())
            while len(active) > target:
                active.remove(max(active))
            candidate = 1 + (self._phase * 3 + index * 5) % slots
            while len(active) < target:
                if candidate not in active:
                    active.add(candidate)
                candidate = 1 + candidate % slots

    def _calculator_block(self) -> CalculatorsBlockData:
        return CalculatorsBlockData(
            pools=tuple(
                CalculatorPoolData(
                    name=name,
                    slots=slots,
                    busy_packs=tuple(sorted(self._calc_active_packs.get(name, set()))),
                )
                for name, slots in zip(self._calc_names, self._calc_total)
            )
        )

    def _worker_block(self) -> WorkerBlockData:
        """Represent one FileOperator service per process-mode simulated Worker."""
        missing = 1 if self._phase and self._phase % 17 == 0 else 0
        alive = max(0, self.workers_alive - missing)
        return WorkerBlockData(
            expected=self.workers_total,
            present=self.workers_alive,
            busy=self.busy,
            idle=max(0, self.workers_alive - self.busy),
            assigned=self.busy,
            file_ops_expected=self.workers_alive,
            file_ops_alive=alive,
            file_ops_missing=missing,
            heartbeat_recent=max(0, self.workers_alive - self.stale),
            heartbeat_stale=self.stale,
        )

    @staticmethod
    def _resize_history(
        history: list[float | None], width: int
    ) -> list[float | None]:
        """Keep newest values and add empty slots on the old/left side."""
        return [*history[:width], *([None] * max(0, width - len(history)))]

    def _build_frame(
        self,
        *,
        rate_per_sec: float | None = None,
        eta_s: int | None = None,
        elapsed_s: int | None = None,
        queue_observed_at: float | None = None,
    ) -> OverviewFrame:
        elapsed_s = (
            max(0, int(time.monotonic() - self._t0))
            if elapsed_s is None
            else elapsed_s
        )
        hours, rem = divmod(elapsed_s, 3600)
        minutes, seconds = divmod(rem, 60)
        if eta_s is None:
            eta = "—"
        else:
            eta_h, eta_r = divmod(eta_s, 3600)
            eta_m, eta_sec = divmod(eta_r, 60)
            eta = f"{eta_h}:{eta_m:02d}:{eta_sec:02d}"
        calcs = tuple(
            (name, busy, total)
            for name, busy, total in zip(
                self._calc_names, self._calc_busy, self._calc_total
            )
        )
        sampler = self._build_sampler_display(eta)
        return OverviewFrame(
            scan="simu-iDM_Vector_V1",
            mode="running",
            elapsed=f"{hours:02d}:{minutes:02d}:{seconds:02d}",
            ref="simu-iDM_Vector_V1",
            redis="127.0.0.1:6379",
            hz=f"{SIMU_HZ:0.1f} Hz",
            target=self.target,
            done=self.done,
            running=self.running,
            failed=self.failed,
            rate=format_sample_rate(rate_per_sec),
            eta=eta,
            avg="—" if rate_per_sec is None else "3.82 s",
            task_q=self.task_q,
            archive_q=self.archive_q,
            feedback_q=self.feedback_q,
            queue_block=self._queue_block(queue_observed_at),
            worker_block=self._worker_block(),
            workers_alive=self.workers_alive,
            workers_total=self.workers_total,
            stale=self.stale,
            busy=self.busy,
            cpu=self.cpu,
            mem_g=self.mem_g,
            mem_total_g=64.0,
            method=sampler.method,
            run_id="run-20260910-120104",
            host="nersc-login-02",
            started="2026-09-10 12:01:04 UTC",
            calculators=calcs,
            calculator_block=self._calculator_block(),
            spark_samples=tuple(self._samples),
            spark_queue=tuple(self._queue),
            fds=self.fds,
            fds_limit=self.fds_limit,
            sampler_generation=self.sampler_generation,
            sampler_max_generations=self.sampler_max_generations,
            sampler_radius=self.sampler_radius,
            sampler_core=self.sampler_core,
            sampler_open=self.sampler_open,
            sampler_state=self.sampler_state,
            sampler_saved=self.sampler_saved,
            sampler=sampler,
        )

    def _build_sampler_display(self, eta: str) -> SamplerDisplay:
        """Return representative, method-specific telemetry for visual QA."""
        method = self.sampler_method
        saved = self.sampler_saved
        finite = {
            "Random": ("finite", 1830, 3000, {"accepted": 1540, "seed": 7}),
            "Grid": (
                "finite",
                1840,
                3000,
                {"shape": "10×10×30", "dimensions": 3},
            ),
            "CSV": (
                "finite",
                1922,
                3000,
                {"accepted": 1840, "basename": "points.csv"},
            ),
            "Bridson": (
                "finite",
                1840,
                3000,
                {"radius": 0.10, "max_attempt": 30},
            ),
        }
        if method in finite:
            kind, current, target, metrics = finite[method]
            return SamplerDisplay(
                method=method,
                family="simple",
                state="running",
                progress=SamplerProgress(
                    kind=kind,
                    current=current,
                    target=target,
                    eta=eta,
                ),
                metrics=metrics,
                saved=saved,
            )
        if method == "AdaptiveBridson":
            return SamplerDisplay(
                method=method,
                family="adaptive",
                state=self.sampler_state,
                progress=SamplerProgress(
                    kind="generation",
                    current=self.sampler_generation,
                    target=self.sampler_max_generations,
                    eta=eta,
                ),
                metrics={
                    "radius": self.sampler_radius,
                    "core": self.sampler_core,
                    "open": self.sampler_open,
                },
                saved=saved,
            )
        if method in {"Dynesty", "MultiNest"}:
            return SamplerDisplay(
                method=method,
                family="nested",
                state="running",
                progress=SamplerProgress(kind="evidence"),
                metrics={
                    "nlive": 500,
                    "niter": 1842,
                    "ncall": 12034,
                    "efficiency": 0.153 if method == "Dynesty" else 0.187,
                    "dlogz_target": 0.14 if method == "Dynesty" else 0.42,
                },
                saved=saved,
            )

        role_count = 4
        kind = "iteration"
        metrics: dict[str, object] = {
            "chains": role_count,
            "accept_rate": 0.31,
        }
        if method == "MCMC":
            metrics["rhat"] = 1.02
        elif method == "ToyMCMC":
            metrics.update(accept_rate=0.34, proposal_scale=0.20)
        elif method == "AMMCMC":
            metrics.update(accept_rate=0.29, adapt_state="active")
        elif method == "DRAM":
            metrics.update(accept_rate=0.44, retry_rate=0.18)
        elif method == "EnsembleMCMC":
            kind = "generation"
            metrics = {"walkers": 16, "accept_rate": 0.42, "stretch_a": 2.0}
        elif method == "DEMCMC":
            kind = "generation"
            metrics = {"walkers": 16, "accept_rate": 0.38, "de_gamma": None}
        elif method == "PTMCMC":
            metrics = {
                "chains": 8,
                "accept_rate": 0.31,
                "swap_accepts": 18,
                "swap_attempts": 24,
            }
        elif method == "PTEnsemble":
            kind = "generation"
            metrics = {
                "chains": 8,
                "accept_rate": 0.36,
                "swap_accepts": 20,
                "swap_attempts": 28,
            }
        return SamplerDisplay(
            method=method,
            family="MCMC",
            state="running",
            progress=SamplerProgress(
                kind=kind,
                current=184,
                target=500,
                eta=eta,
            ),
            metrics=metrics,
            saved=saved,
        )


def simu_choice() -> ScanChoice:
    return simulated_choice()


def simu_overview() -> OverviewFrame:
    return SimuEngine().empty_frame()


__all__ = [
    "OverviewFrame",
    "SimuEngine",
    "simu_choice",
    "simu_overview",
]
