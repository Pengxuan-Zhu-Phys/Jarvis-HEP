"""Finite samplers recover from terminal feedback, independently of archiving."""

from __future__ import annotations

import os
import time
from abc import ABC, abstractmethod
from collections.abc import Callable, Mapping
from copy import deepcopy
from typing import Any
from uuid import uuid4

import numpy as np

from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_checkpoint_config, get_runtime_block
from jarvishep2.sample import Sample
from jarvishep2.sampling.sampler import SamplingVirtual
from jarvishep2.sampling.sample_checkpoint import (
    FORMAT, VERSION, atomic_pickle, config_identity, read_pickle,
)
from jarvishep2.runtime.checkpoint_process import CheckpointProcess
from jarvishep2.sampling.stateless_batch import deterministic_sampler_uuid
from jarvishep2.sampling.variables import load_variables


class CompletionSampler(SamplingVirtual, ABC):
    """One owner thread generates points, consumes feedback and takes snapshots.

    Every accepted proposal is registered before it can leave the sampler.
    [0, next_sample_index) minus unfinished is exactly the completed set.
    Only unfinished entries consume memory; no lifetime UUID collection.
    """

    uses_sample_checkpoint = True
    uuid_prefix = "sample"

    def __init__(self) -> None:
        super().__init__()
        self.vars = []
        self._mapper_pipeline = None
        self._selectionexp = None
        self._seed = 0
        self._index = 0
        self._accepted_index = 0
        self._batch_size = 16
        self._max_inflight = 16
        self._unfinished: dict[int, dict[str, Any]] = {}
        self._uuid_to_unfinished: dict[str, int] = {}
        self._completed_count = 0
        self._failed_count = 0
        self._resume_pending: list[int] = []
        self._interrupt_check: Callable[[], bool] | None = None
        self._control_next_check = 0.0
        self._feedback_timeout = 3600.0
        self._last_completion_at = 0.0
        self._pickle_writer: CheckpointProcess | None = None
        self._owns_checkpoint = False
        self._config_id = ""
        self._dataset_id: str | None = None
        self._logger = get_jarvis_logger("sampler")

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = self.config.get("Sampling") or {}
        bounds = sampling.get("Bounds") or {}
        runtime = get_runtime_block(self.config)
        csv_source = str(getattr(self, "method", "")).lower() == "csv"
        self.vars = [] if csv_source else load_variables(self.config)
        if not csv_source:
            from jarvishep2.mapper import MapperPipeline

            self._mapper_pipeline = MapperPipeline.from_config(self.config)
        self._selectionexp = sampling.get("selection")
        self._seed = int(bounds.get("seed", 0) or 0)
        workers = max(1, int(runtime.get("workers", 1) or 1))
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        self._max_inflight = max(1, int(
            runtime.get("max_inflight", bounds.get("max_worker"))
            or max(self._batch_size, workers * 4)
        ))
        self._config_id = config_identity(self.config)

    def set_interrupt_check(self, check: Callable[[], bool]) -> None:
        self._interrupt_check = check

    def check_runtime_control(self) -> None:
        if getattr(self, "_interrupt_requested", False) or (
            self._interrupt_check and self._interrupt_check()
        ):
            raise KeyboardInterrupt
        now = time.monotonic()
        if self.redis is None or now < self._control_next_check:
            return
        self._control_next_check = now + 0.1
        board = self.redis.read_proc_board("core") or {}
        if isinstance(board, Mapping) and board.get("scan_mode") in {"paused", "stopping", "draining"}:
            raise RuntimeError(f"sampling aborted: scan_mode={board['scan_mode']}")

    def _uuid_for_accepted_index(self, index: int) -> str:
        return deterministic_sampler_uuid(prefix=self.uuid_prefix, seed=self._seed, sample_index=index)

    def _remember(self, sample: Sample, *, point_index: int | None = None) -> Sample:
        index = sample.sample_index
        if index is None or index != self._accepted_index - 1:
            raise ValueError("proposal must have the next contiguous sample index")
        if sample.uuid in self._uuid_to_unfinished:
            raise ValueError(f"duplicate unfinished sample UUID: {sample.uuid}")
        entry: dict[str, Any] = {"uuid": sample.uuid}
        if point_index is None:
            entry.update(u_coords=sample.u_coords.copy(), opera_params=dict(sample.opera_params))
        else:
            entry["point_index"] = point_index
        self._unfinished[index] = entry
        self._uuid_to_unfinished[sample.uuid] = index
        return sample

    def _rebuild_unfinished(self, index: int) -> Sample:
        entry = self._unfinished[index]
        if "point_index" in entry:
            coords = self._point_coords(entry["point_index"])
        else:
            coords = entry["u_coords"]
        sample = self._build_sample(coords)
        sample.uuid = entry["uuid"]
        sample.sample_index = index
        sample.opera_params = dict(entry.get("opera_params") or {})
        return sample

    def accept_completion(self, record: Mapping[str, Any]) -> bool:
        """Only a matched terminal feedback removes an unfinished point."""
        uuid = str(record.get("uuid") or "")
        index = self._uuid_to_unfinished.get(uuid)
        if index is None:
            return False  # Duplicate or unrelated feedback cannot advance progress.
        if record.get("sample_index") != index:
            raise ValueError("completion feedback sample index differs from its UUID")
        status = record.get("terminal_status")
        if status not in {"Completed", "Failed"}:
            raise ValueError("completion feedback requires a terminal status")
        if status == "Failed":
            self._failed_count += 1
        else:
            self._completed_count += 1
        del self._unfinished[index]
        del self._uuid_to_unfinished[uuid]
        return True

    def _drain_completions(self) -> int:
        # One bounded Redis round trip per tick, even for very cheap samples.
        return sum(self.accept_completion(record) for record in
                   self.redis.poll_feedback_batch(self._max_inflight))

    def _sampling_tick(self) -> None:
        self.check_runtime_control()
        if self.redis is not None:
            if self._drain_completions():
                self._last_completion_at = time.monotonic()
            elif self._unfinished and time.monotonic() - self._last_completion_at > self._feedback_timeout:
                raise TimeoutError("sampler timed out waiting for terminal feedback")
        self.save_checkpoint()

    def flush_batch(self, batch: list[Sample]) -> int:
        if not batch:
            return 0
        # propose_next already registered every entry before dispatch.
        if len(batch) == 1:
            self._submit(batch[0])
        else:
            self._submit_group(batch)
        return len(batch)

    def run_distributed(self, *, timeout: float = 3600.0) -> int:
        if self.redis is None:
            raise RuntimeError("sampler requires a Redis queue")
        self._feedback_timeout = max(1.0, float(timeout))
        self._last_completion_at = time.monotonic()
        pushed = 0
        # Resume can reduce worker/window counts. Requeue in bounded batches
        # without generating any new work until the saved window is drained.
        pending = iter(self._resume_pending)
        replay = next(pending, None)
        active: set[int] = set()
        while replay is not None or active:
            self._sampling_tick()
            active.intersection_update(self._unfinished)
            batch = []
            while replay is not None and len(active) < self._max_inflight and len(batch) < self._batch_size:
                batch.append(self._rebuild_unfinished(replay))
                active.add(replay)
                replay = next(pending, None)
            pushed += self.flush_batch(batch)
            if not batch and active:
                time.sleep(0.01)
        self._resume_pending = []

        exhausted = False
        while not exhausted or self._unfinished:
            self._sampling_tick()
            room = self._max_inflight - len(self._unfinished)
            batch = []
            while not exhausted and len(batch) < min(self._batch_size, room):
                sample = self.propose_next()
                if sample is None:
                    exhausted = True
                    break
                batch.append(sample)
            pushed += self.flush_batch(batch)
            if not batch and self._unfinished:
                time.sleep(0.01)
        self.save_checkpoint(force=True)
        return pushed

    def checkpoint_state(self) -> dict[str, Any]:
        """A detached O(unfinished) snapshot, called only between transitions."""
        return deepcopy({
            "next_sample_index": self._accepted_index,
            "source_state": self._source_state(),
            "dataset_id": self._dataset_id,
            "unfinished": self._unfinished,
            "completed_count": self._completed_count,
            "failed_count": self._failed_count,
        })

    def restore_checkpoint(self, state: Mapping[str, Any]) -> None:
        next_index = state["next_sample_index"]
        completed = state["completed_count"]
        failed = state["failed_count"]
        unfinished = state["unfinished"]
        if any(type(value) is not int or value < 0 for value in (next_index, completed, failed)):
            raise ValueError("invalid sampler completion counters")
        if not isinstance(unfinished, dict) or completed + failed + len(unfinished) != next_index:
            raise ValueError("sampler completion counters do not partition the generated points")
        uuids = {}
        for index, entry in unfinished.items():
            if type(index) is not int or not 0 <= index < next_index:
                raise ValueError("unfinished index is outside the generated range")
            if not isinstance(entry, dict):
                raise ValueError("unfinished sample must be a mapping")
            uuid = entry.get("uuid")
            if not isinstance(uuid, str) or not uuid or uuid in uuids:
                raise ValueError("invalid/duplicate unfinished UUID")
            uuids[uuid] = index
        if state["dataset_id"] != self._dataset_id:
            raise ValueError("all.pickle and runtime.pickle refer to different point sets")
        self._restore_source_state(state["source_state"])
        if next_index > self._source_position():
            raise ValueError("sample count is ahead of the source cursor")
        for entry in unfinished.values():
            if isinstance(self, PointSetSampler):
                point_index = entry.get("point_index")
                if type(point_index) is not int or not 0 <= point_index < self._index:
                    raise ValueError("unfinished point is outside the generated point set")
                if entry["uuid"] != self._point_uuids[point_index]:
                    raise ValueError("unfinished UUID differs from all.pickle")
            elif (not isinstance(entry.get("u_coords"), np.ndarray)
                    or entry["u_coords"].shape != (len(self.vars),)
                    or not isinstance(entry.get("opera_params"), dict)):
                raise ValueError("invalid unfinished sample data")
        self._accepted_index = next_index
        self._completed_count = completed
        self._failed_count = failed
        self._unfinished = deepcopy(unfinished)
        self._uuid_to_unfinished = uuids
        self._resume_pending = sorted(unfinished)

    def configure_checkpoint(self, *, checkpoint_file: str, writer: CheckpointProcess | None = None) -> None:
        policy = get_checkpoint_config(self.config)
        if not policy["enabled"] or self.config.get("_check_modules_sample_layout"):
            return
        if not checkpoint_file:
            raise ValueError("sampler checkpoint path is required")
        self._prepare_checkpoint_data(checkpoint_file)
        self._owns_checkpoint = writer is None
        if writer is None:
            try:
                writer = CheckpointProcess(checkpoint_file, self._logger,
                                           interval=policy["heartbeat"],
                                           scan_name=str(self.config.get("scan_name") or ""))
            except Exception as exc:
                self._logger.warning("Checkpoint process unavailable -> %s", exc)
                return
            writer.start()
        self._pickle_writer = writer
        self.save_checkpoint(force=True)

    def save_checkpoint(self, *, force: bool = False) -> bool:
        writer = self._pickle_writer
        if writer is None or (not force and not writer.requested()):
            return False
        if force:
            writer.wait()
        if not writer.available():
            return False
        state = self.checkpoint_state()
        return writer.submit({
            "format": FORMAT,
            "version": VERSION,
            "method": self.method.lower(),
            "config_id": self._config_id,
            "integrity": {"unit_execution": (self.config.get("_unit_execution") or {}).get("execution_identity")},
            "sampler_state": state,
        })

    def shutdown_checkpointing(self) -> None:
        if self._pickle_writer is not None:
            if self._owns_checkpoint:
                self._pickle_writer.finish()
            self._pickle_writer = None

    def _prepare_checkpoint_data(self, checkpoint_file: str) -> None:
        pass

    @abstractmethod
    def propose_next(self) -> Sample | None:
        """Generate and register one accepted point before returning it."""

    @abstractmethod
    def _source_state(self) -> dict[str, Any]:
        """Read only the generator cursor and any necessary RNG state."""

    @abstractmethod
    def _restore_source_state(self, state: Mapping[str, Any]) -> None:
        """Restore the source after separately loading its fixed data."""

    def _source_position(self) -> int:
        return self._index

    def log_settings_rows(self) -> list[tuple[str, Any]]:
        return [("selection", self._selectionexp)] if self._selectionexp else []

    def log_summary_rows(self) -> list[tuple[str, Any]]:
        return [("points generated", self._accepted_index),
                ("completed", self._completed_count), ("failed", self._failed_count),
                ("unfinished", len(self._unfinished))]


class PointSetSampler(CompletionSampler):
    """Immutable Grid/Bridson coordinates and UUIDs live in all.pickle once."""

    def __init__(self) -> None:
        super().__init__()
        self._point_uuids: tuple[str, ...] | None = None

    def _point_uuid(self, index: int) -> str:
        if self._point_uuids is not None:
            return self._point_uuids[index]
        return self._uuid_for_accepted_index(index)

    def _prepare_checkpoint_data(self, checkpoint_file: str) -> None:
        if self._dataset_id is not None:
            return  # Already loaded from all.pickle during resume.
        self._ensure_grid()
        self._dataset_id = uuid4().hex
        self._point_uuids = tuple(self._point_uuid(i) for i in range(len(self._P)))
        atomic_pickle(os.path.join(os.path.dirname(checkpoint_file), "all.pickle"), {
            "format": FORMAT, "version": VERSION, "config_id": self._config_id,
            "dataset_id": self._dataset_id, "coordinates": self._P,
            "uuids": self._point_uuids,
        })

    def prepare_resume_points(self, checkpoint_file: str, state: Mapping[str, Any]) -> None:
        payload = read_pickle(os.path.join(os.path.dirname(checkpoint_file), "all.pickle"))
        if (payload.get("format") != FORMAT or payload.get("version") != VERSION
                or payload.get("config_id") != self._config_id
                or not payload.get("dataset_id") or payload["dataset_id"] != state.get("dataset_id")):
            raise ValueError("all.pickle does not belong to this runtime/configuration")
        points, uuids = payload["coordinates"], payload["uuids"]
        if (not isinstance(points, np.ndarray) or points.ndim != 2
                or points.shape[1] != len(self.vars) or len(uuids) != len(points)):
            raise ValueError("all.pickle has invalid coordinates/UUIDs")
        self._P, self._point_uuids = points, tuple(uuids)
        self._dataset_id = payload["dataset_id"]
        self.info["NSamples"] = len(points)

    def _source_state(self) -> dict[str, Any]:
        return {"index": self._index}

    def _restore_source_state(self, state: Mapping[str, Any]) -> None:
        index = state["index"]
        if type(index) is not int or not 0 <= index <= len(self._P):
            raise ValueError("point-set cursor is outside all.pickle")
        self._index = index

    def _point_coords(self, index: int) -> np.ndarray:
        return np.asarray(self._P[index], dtype=np.float64)
