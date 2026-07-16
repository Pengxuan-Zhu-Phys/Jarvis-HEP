#!/usr/bin/env python3
"""D13.5 Dynesty + RedisEvaluationPool + vendored UUID channel tests."""

from __future__ import annotations

import threading
import unittest
from typing import Any

import numpy as np

from jarvishep2.Sampling.Source.Dynesty.py.dynesty.jarvis_uuid import (
    looks_uuid_augmented,
    split_vid_batch,
)
from jarvishep2.Sampling.dynesty_sampler import DynestySampler, _jarvis_prior_transform
from jarvishep2.Sampling.redis_evaluation_pool import RedisEvaluationPool
from jarvishep2.distributor import Distributor, STATELESS_METHODS
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.sample import Sample


class JarvisUuidHelpersTests(unittest.TestCase):
    def test_prior_transform_appends_uuid(self) -> None:
        u = np.array([0.1, 0.2, 0.3])
        out = _jarvis_prior_transform(u)
        self.assertTrue(looks_uuid_augmented(out))
        self.assertEqual(len(out), 4)
        np.testing.assert_allclose(np.asarray(out[:-1], dtype=float), u)

    def test_split_vid_batch(self) -> None:
        rows = [_jarvis_prior_transform(np.array([0.5, 0.5])) for _ in range(3)]
        v, uid, full = split_vid_batch(rows)
        self.assertEqual(v.shape, (3, 2))
        self.assertIsNotNone(uid)
        self.assertEqual(len(uid), 3)
        self.assertEqual(len(full), 3)


class RedisPoolUnitTests(unittest.TestCase):
    def test_prior_transform_is_local(self) -> None:
        queue = make_fakeredis_queue()
        calls = {"n": 0}

        def build_sample(payload, uuid):
            calls["n"] += 1
            return Sample(uuid=uuid, u_coords=np.asarray(payload, dtype=float))

        pool = RedisEvaluationPool(
            queue, build_sample=build_sample, batch_size=4, seed=1
        )
        items = [np.array([0.1, 0.2]), np.array([0.3, 0.4])]
        out = pool.map(_jarvis_prior_transform, items)
        self.assertEqual(len(out), 2)
        self.assertTrue(looks_uuid_augmented(out[0]))
        self.assertEqual(calls["n"], 0)  # no Sample built for prior

    def test_logl_batch_via_redis(self) -> None:
        queue = make_fakeredis_queue()

        def build_sample(payload, uuid):
            s = Sample(uuid=uuid, u_coords=np.asarray(payload, dtype=float))
            return s

        pool = RedisEvaluationPool(
            queue, build_sample=build_sample, batch_size=4, seed=7, timeout=10.0
        )
        stop = threading.Event()

        def worker() -> None:
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is None:
                    continue
                u = np.asarray(task.get("u_coords") or [], dtype=float)
                logl = -float(np.sum((u - 0.5) ** 2))
                queue.publish_feedback(
                    {
                        "uuid": task["uuid"],
                        "status": "Completed",
                        "observables": {"LogL": logl},
                    }
                )

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        try:
            # Fake loglike name so pool treats as remote
            def loglikelihood(x):
                raise AssertionError("should not call local loglike")

            loglikelihood.__name__ = "loglikelihood"
            items = [
                _jarvis_prior_transform(np.array([0.5, 0.5])),
                _jarvis_prior_transform(np.array([0.0, 0.0])),
            ]
            out = pool.map(loglikelihood, items)
            self.assertEqual(len(out), 2)
            self.assertAlmostEqual(float(out[0].val), 0.0, places=6)
            self.assertLess(float(out[1].val), 0.0)
        finally:
            stop.set()
            thread.join(timeout=2.0)


class DistributorDynestyTests(unittest.TestCase):
    def test_registered(self) -> None:
        self.assertNotIn("Dynesty", STATELESS_METHODS)
        s = Distributor.set_method("Dynesty")
        self.assertIsInstance(s, DynestySampler)


class JarvisLoggerBridgeTests(unittest.TestCase):
    def test_progress_and_warnings_use_injected_logger(self) -> None:
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import jarvis_logging as jl
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty.utils import print_fn

        class CapturingLogger:
            def __init__(self) -> None:
                self.infos: list[str] = []
                self.warnings: list[str] = []
                self.extra = {"jarvis_module": "Dynesty.Inner"}

            def bind(self, **kwargs):
                return self

            def info(self, msg, *a, **k):
                self.infos.append(str(msg))

            def warning(self, msg, *a, **k):
                self.warnings.append(str(msg))

        capt = CapturingLogger()
        jl.set_dynesty_logger(capt)
        try:
            # Minimal results-like object for print_fn
            results = type(
                "R",
                (),
                {
                    "loglstar": -1.0,
                    "logz": -10.0,
                    "logzvar": 0.01,
                    "delta_logz": 0.5,
                    "bounditer": 0,
                    "nc": 3,
                    "eff": 12.0,
                },
            )()
            print_fn(results, niter=10, ncall=50, dlogz=0.01, logger=capt)
            self.assertTrue(capt.infos or capt.warnings)
            text = "\n".join(capt.infos + capt.warnings)
            self.assertIn("Progress", text)
            self.assertIn("logz", text.lower())

            jl.emit_warning("test plateau warning", UserWarning)
            self.assertTrue(any("plateau" in w or "test" in w for w in capt.warnings))
        finally:
            jl.set_dynesty_logger(None)


class DynestyNestedSmokeTests(unittest.TestCase):
    def test_static_nested_with_redis_pool(self) -> None:
        """Small NestedSampler run (nlive small) through RedisEvaluationPool."""
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import NestedSampler

        queue = make_fakeredis_queue()
        stop = threading.Event()

        def worker() -> None:
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is None:
                    continue
                u = np.asarray(task.get("u_coords") or [], dtype=float)
                # Peak at centre of unit cube
                logl = -0.5 * float(np.sum((u - 0.5) ** 2) / 0.05**2)
                queue.publish_feedback(
                    {
                        "uuid": task["uuid"],
                        "status": "Completed",
                        "observables": {"LogL": logl},
                    }
                )

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()

        def build_sample(payload, uuid):
            return Sample(uuid=uuid, u_coords=np.asarray(payload, dtype=float))

        pool = RedisEvaluationPool(
            queue, build_sample=build_sample, batch_size=8, seed=3, timeout=60.0
        )

        def loglike(x):
            # Should not be used when pool handles logl
            arr = np.asarray(x, dtype=object).reshape(-1)
            u = np.asarray(arr[:-1], dtype=float) if looks_uuid_augmented(arr) else np.asarray(arr, dtype=float)
            return -0.5 * float(np.sum((u - 0.5) ** 2) / 0.05**2)

        try:
            sampler = NestedSampler(
                loglikelihood=loglike,
                prior_transform=_jarvis_prior_transform,
                ndim=2,
                nlive=20,
                pool=pool,
                rstate=np.random.default_rng(0),
                queue_size=4,
            )
            sampler.run_nested(dlogz=1.0, maxiter=40, print_progress=False)
            res = sampler.results
            self.assertIn("logz", res.keys())
            self.assertGreater(len(res["logz"]), 0)
            # UUID channel should produce samples_uid when dead points saved
            if "samples_uid" in res.keys():
                uids = list(res["samples_uid"])
                self.assertTrue(any(str(u) for u in uids))
        finally:
            stop.set()
            thread.join(timeout=2.0)


if __name__ == "__main__":
    unittest.main()
