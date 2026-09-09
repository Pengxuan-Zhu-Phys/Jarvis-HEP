#!/usr/bin/env python3
"""FeedbackSampler base unit tests (D13.1)."""

from __future__ import annotations

import threading
import time
import unittest
from collections.abc import Mapping, Sequence
from typing import Any
from unittest import mock

import numpy as np

from jarvishep2.Sampling.adaptive_bridson import AdaptiveBridsonSampler
from jarvishep2.Sampling.feedback_sampler import FeedbackSampler
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.sample import Sample


class _ToyFeedbackSampler(FeedbackSampler):
    """Minimal FeedbackSampler: two generations of fixed u-points."""

    method = "ToyFeedback"

    def __init__(self) -> None:
        super().__init__()
        self._init_seed_sequence(42)
        self._batch_size = 4
        self._plan: list[list[str]] = [["a", "b"], ["c"]]
        self.proposals_seen: list[list[str]] = []
        self.absorbed: list[list[str]] = []
        self._step = 0

    def propose_generation(self) -> Sequence[Sample] | None:
        if self._step >= len(self._plan):
            return None
        uuids = list(self._plan[self._step])
        self.proposals_seen.append(uuids)
        samples = []
        for uid in uuids:
            sample = self._build_sample(np.array([0.1, 0.2], dtype=np.float64))
            sample.uuid = uid
            samples.append(sample)
        self._step += 1
        return samples

    def absorb_generation(self, results: Sequence[Mapping[str, Any]]) -> None:
        self.absorbed.append([str(r.get("uuid", "")) for r in results])


class FeedbackSamplerUnitTests(unittest.TestCase):
    def test_als_is_feedback_sampler(self) -> None:
        self.assertTrue(issubclass(AdaptiveBridsonSampler, FeedbackSampler))
        sampler = AdaptiveBridsonSampler()
        self.assertIsInstance(sampler, FeedbackSampler)
        self.assertTrue(sampler.at_safe_barrier())

    def test_generation_rng_is_deterministic(self) -> None:
        a = _ToyFeedbackSampler()
        b = _ToyFeedbackSampler()
        va = a._generation_rng(3).random(5)
        vb = b._generation_rng(3).random(5)
        np.testing.assert_array_equal(va, vb)
        # Different generations diverge
        v0 = a._generation_rng(0).random(3)
        v1 = a._generation_rng(1).random(3)
        self.assertFalse(np.array_equal(v0, v1))

    def test_generic_loop_propose_publish_absorb(self) -> None:
        queue = make_fakeredis_queue()
        sampler = _ToyFeedbackSampler()
        sampler.set_config({"Runtime": {"mode": "redis", "workers": 1}})
        sampler.set_redis(queue)

        # Worker mock: drain tasks and push feedback in a side thread.
        stop = threading.Event()

        def worker() -> None:
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is None:
                    continue
                queue.publish_feedback(
                    {
                        "uuid": task["uuid"],
                        "status": "Completed",
                        "observables": {"y": 1.0},
                    }
                )
                queue.add_archived_uuids("scan", [str(task["uuid"])])

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        try:
            n = sampler.run_adaptive(timeout=30.0)
        finally:
            stop.set()
            thread.join(timeout=2.0)

        self.assertEqual(n, 3)
        self.assertEqual(sampler.proposals_seen, [["a", "b"], ["c"]])
        # Absorb sees both generations (order within a gen may vary)
        self.assertEqual(len(sampler.absorbed), 2)
        self.assertEqual(sorted(sampler.absorbed[0]), ["a", "b"])
        self.assertEqual(sampler.absorbed[1], ["c"])
        self.assertTrue(sampler.at_safe_barrier())
        self.assertEqual(sampler._generation, 2)

    def test_wait_timeout_raises(self) -> None:
        queue = make_fakeredis_queue()
        sampler = _ToyFeedbackSampler()
        sampler.set_config({"Runtime": {"mode": "redis"}})
        sampler.set_redis(queue)
        sampler._register_pending("stuck-uuid")
        with self.assertRaises(TimeoutError):
            sampler.wait_for_generation(timeout=1.0)

    def test_unmatched_feedback_is_logged(self) -> None:
        queue = make_fakeredis_queue()
        sampler = _ToyFeedbackSampler()
        sampler.set_config({"Runtime": {"mode": "redis"}})
        sampler.set_redis(queue)
        sampler._register_pending("expected")

        class CapturingLogger:
            def __init__(self) -> None:
                self.warnings: list[str] = []

            def warning(self, message: str, *args: Any) -> None:
                self.warnings.append(message % args if args else message)

        logger = CapturingLogger()
        sampler._logger = logger
        queue.publish_feedback({"uuid": "stray", "logL": -1.0})
        queue.publish_feedback({"uuid": "expected", "logL": -2.0})
        sampler.wait_for_generation(timeout=2.0)
        self.assertTrue(any("unmatched hep:feedback" in item for item in logger.warnings))

    def test_feedback_export_import_roundtrip(self) -> None:
        sampler = _ToyFeedbackSampler()
        sampler._generation = 4
        sampler._pending_uuids = {"u1", "u2"}
        sampler._submitted_uuids = ["u0", "u1", "u2"]
        sampler.set_persisted_uuids({"u0"})
        sampler._on_failure = "halt"
        blob = sampler._feedback_export_state()
        restored = _ToyFeedbackSampler()
        restored._feedback_import_state(blob)
        self.assertEqual(restored._generation, 4)
        self.assertEqual(restored._pending_uuids, {"u1", "u2"})
        self.assertEqual(restored._on_failure, "halt")
        self.assertEqual(restored._seed, 42)

    def test_checkpoint_barrier_uses_scan_block_name_without_loader_alias(self) -> None:
        queue = make_fakeredis_queue()
        sampler = _ToyFeedbackSampler()
        sampler.set_config(
            {"Scan": {"name": "raw-api-scan"}, "Runtime": {"mode": "redis"}}
        )
        sampler.set_redis(queue)
        sampler._submitted_uuids = ["u0"]
        queue.add_archived_uuids("raw-api-scan", ["u0"])
        sampler._save_checkpoint_callback = lambda **_kwargs: True

        self.assertTrue(sampler.checkpoint_at_barrier())
        self.assertEqual(sampler._persisted_uuids, {"u0"})

    def test_never_published_scan_mode_does_not_abort_generation(self) -> None:
        queue = make_fakeredis_queue()
        sampler = _ToyFeedbackSampler()
        sampler.set_config({"Runtime": {"mode": "redis"}})
        sampler.set_redis(queue)
        sampler._register_pending("expected")
        clock = {"t": 0.0}

        def fake_mono() -> float:
            return clock["t"]

        def fake_pull(timeout=1, **_kwargs):
            clock["t"] += float(timeout)
            if clock["t"] >= 90.0:
                return {"uuid": "expected", "logL": -1.0}
            return None

        with mock.patch(
            "jarvishep2.sampling.feedback_sampler.time.monotonic", fake_mono
        ):
            sampler.redis.pull_feedback = fake_pull  # type: ignore[method-assign]
            results = sampler.wait_for_generation(timeout=120.0)
        self.assertEqual(results[0]["uuid"], "expected")

    def test_stale_sec_120_unknown_fail_grace_is_240_not_60(self) -> None:
        queue = make_fakeredis_queue()
        sampler = _ToyFeedbackSampler()
        sampler.set_config(
            {
                "EnvReqs": {
                    "V2": {"factory": {"watchdog": {"stale_sec": 120}}}
                }
            }
        )
        sampler.set_redis(queue)
        sampler._register_pending("stuck")
        queue.publish_proc_board("core", scan_mode="running")
        clock = {"t": 0.0}
        seen = {"n": 0}

        def fake_mono() -> float:
            return clock["t"]

        def fake_pull(timeout=1, **_kwargs):
            clock["t"] += float(timeout)
            return None

        real_read = queue.read_proc_board

        def fake_read(role, **kwargs):
            seen["n"] += 1
            if seen["n"] == 1:
                return real_read(role, **kwargs)
            return {}

        with mock.patch(
            "jarvishep2.sampling.feedback_sampler.time.monotonic", fake_mono
        ):
            sampler.redis.pull_feedback = fake_pull  # type: ignore[method-assign]
            sampler.redis.read_proc_board = fake_read  # type: ignore[method-assign]
            with self.assertRaises(RuntimeError) as ctx:
                sampler.wait_for_generation(timeout=1000.0)
        self.assertIn("unknown", str(ctx.exception).lower())
        self.assertGreaterEqual(clock["t"], 240.0)
        self.assertLess(clock["t"], 300.0)

    def test_paused_scan_mode_ends_wait_for_generation_within_one_second(self) -> None:
        queue = make_fakeredis_queue()
        sampler = _ToyFeedbackSampler()
        sampler.set_config({"Runtime": {"mode": "redis"}})
        sampler.set_redis(queue)
        sampler._register_pending("stuck")
        queue.publish_proc_board(
            "core", scan_mode="paused", pause_reason="death_rate"
        )
        started = time.monotonic()
        with self.assertRaises(RuntimeError) as ctx:
            sampler.wait_for_generation(timeout=30.0)
        elapsed = time.monotonic() - started
        self.assertLess(elapsed, 1.5)
        self.assertIn("paused", str(ctx.exception).lower())
        self.assertIn("death_rate", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
