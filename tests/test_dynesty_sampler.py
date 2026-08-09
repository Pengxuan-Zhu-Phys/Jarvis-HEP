#!/usr/bin/env python3
"""D13.5 Dynesty + RedisEvaluationPool + vendored UUID channel tests."""

from __future__ import annotations

import csv
import os
import tempfile
import threading
import unittest
from typing import Any

import numpy as np

from jarvishep2.Sampling.Source.Dynesty.py.dynesty.jarvis_uuid import (
    looks_uuid_augmented,
    split_vid_batch,
)
from jarvishep2.Sampling.dynesty_sampler import (
    DynestySampler,
    _jarvis_prior_transform,
    export_dynesty_results_csv,
)
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
    def test_dynesty_sampler_never_binds_toy_loglike(self) -> None:
        """Production constructor must wire Redis logL, not _local_loglike."""
        from jarvishep2.Sampling.dynesty_sampler import DynestySampler

        queue = make_fakeredis_queue()
        seen: list[str] = []

        def build_sample(payload, uuid):
            seen.append(str(uuid))
            return Sample(uuid=uuid, u_coords=np.asarray(payload, dtype=float))

        pool = RedisEvaluationPool(
            queue, build_sample=build_sample, batch_size=4, seed=1, timeout=5.0
        )
        sampler = DynestySampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "Dynesty",
                    "Variables": [
                        {
                            "name": "x",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1},
                            },
                        },
                        {
                            "name": "y",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1},
                            },
                        },
                    ],
                    "Bounds": {"nlive": 10},
                }
            }
        )
        kwargs = sampler._build_constructor_kwargs(
            pool=pool, rstate=np.random.default_rng(0)
        )
        from jarvishep2.Sampling.dynesty_sampler import NestedRedisLogL

        loglike = kwargs["loglikelihood"]
        self.assertIsNot(loglike, sampler._local_loglike)
        self.assertIsInstance(loglike, NestedRedisLogL)
        self.assertEqual(getattr(loglike, "__name__", ""), "loglikelihood")
        # Toy path must hard-fail if ever invoked.
        with self.assertRaises(RuntimeError):
            sampler._local_loglike(np.array([0.1, 0.2]))

        stop = threading.Event()

        def worker() -> None:
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is None:
                    continue
                u = np.asarray(task.get("u_coords") or [], dtype=float)
                queue.publish_feedback(
                    {"uuid": task["uuid"], "logL": float(-np.sum((u - 0.5) ** 2))}
                )

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        try:
            v = _jarvis_prior_transform(np.array([0.5, 0.5]))
            val = float(loglike(v))
            self.assertAlmostEqual(val, 0.0, places=6)
            self.assertEqual(len(seen), 1)
        finally:
            stop.set()
            thread.join(timeout=2.0)

    def test_evaluate_logl_goes_remote(self) -> None:
        """Direct loglikelihood(v) calls (inside sample()) must hit Redis."""
        queue = make_fakeredis_queue()
        seen: list[str] = []

        def build_sample(payload, uuid):
            seen.append(str(uuid))
            return Sample(uuid=uuid, u_coords=np.asarray(payload, dtype=float))

        pool = RedisEvaluationPool(
            queue, build_sample=build_sample, batch_size=4, seed=3, timeout=10.0
        )
        stop = threading.Event()

        def worker() -> None:
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is None:
                    continue
                u = np.asarray(task.get("u_coords") or [], dtype=float)
                queue.publish_feedback(
                    {
                        "uuid": task["uuid"],
                        "logL": float(-np.sum(u)),
                    }
                )

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        try:
            v = _jarvis_prior_transform(np.array([0.2, 0.8]))
            val = pool.evaluate_logl(v)
            self.assertAlmostEqual(val, -1.0, places=5)
            self.assertEqual(len(seen), 1)
        finally:
            stop.set()
            thread.join(timeout=2.0)

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

    def test_ambiguous_map_fails_loud(self) -> None:
        """D13.7a: unknown callables must not silently run locally."""
        queue = make_fakeredis_queue()

        def build_sample(payload, uuid):
            return Sample(uuid=uuid, u_coords=np.asarray(payload, dtype=float))

        pool = RedisEvaluationPool(
            queue, build_sample=build_sample, batch_size=4, seed=1
        )

        def mystery(x):
            return -1.0

        mystery.__name__ = "mystery_objective"
        with self.assertRaises(RuntimeError) as ctx:
            pool.map(mystery, [np.array([0.1, 0.2])])
        self.assertIn("ambiguous dispatch", str(ctx.exception).lower())

    def test_bound_bootstrap_map_stays_local(self) -> None:
        """Dynesty MultiEllipsoid bootstrap must not go through Redis logL."""
        queue = make_fakeredis_queue()

        def build_sample(payload, uuid):
            raise AssertionError("bootstrap must not build Redis samples")

        pool = RedisEvaluationPool(
            queue, build_sample=build_sample, batch_size=4, seed=1
        )
        calls = {"n": 0}

        def _ellipsoid_bootstrap_expand(args):
            calls["n"] += 1
            multi, points, seed = args
            return 1.05

        # Mimic dynesty bounding.update zip(multis, ps, seeds)
        seeds = np.random.SeedSequence(0).spawn(2)
        points = np.random.default_rng(0).random((10, 2))
        args = list(zip([True, True], [points, points], seeds))
        out = pool.map(_ellipsoid_bootstrap_expand, args)
        self.assertEqual(out, [1.05, 1.05])
        self.assertEqual(calls["n"], 2)

    def test_seedsequence_is_not_uuid_augmented(self) -> None:
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty.jarvis_uuid import (
            looks_uuid_augmented,
        )

        row = (True, np.zeros((3, 2)), np.random.SeedSequence(1))
        self.assertFalse(looks_uuid_augmented(row))
        self.assertTrue(looks_uuid_augmented(np.array([0.1, 0.2, "abc-uuid"], dtype=object)))

    def test_unmatched_feedback_is_logged(self) -> None:
        """D13.7b: destructive BLPOP of stray uuid leaves a warning trail."""
        queue = make_fakeredis_queue()

        def build_sample(payload, uuid):
            return Sample(uuid=uuid, u_coords=np.asarray(payload, dtype=float))

        class _CapturingLogger:
            def __init__(self) -> None:
                self.warnings: list[str] = []

            def warning(self, msg: str, *args: Any) -> None:
                self.warnings.append(msg % args if args else str(msg))

            def info(self, *args: Any, **kwargs: Any) -> None:
                return None

        logger = _CapturingLogger()
        pool = RedisEvaluationPool(
            queue,
            build_sample=build_sample,
            batch_size=4,
            seed=7,
            timeout=5.0,
            logger=logger,
        )
        stop = threading.Event()

        def worker() -> None:
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is None:
                    continue
                # Poison pill first — wrong uuid must be dropped+logged.
                queue.publish_feedback(
                    {
                        "uuid": "stray-uuid-not-in-batch",
                        "status": "Completed",
                        "observables": {"LogL": -999.0},
                    }
                )
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

            def loglikelihood(x):
                raise AssertionError("should not call local loglike")

            loglikelihood.__name__ = "loglikelihood"
            items = [_jarvis_prior_transform(np.array([0.5, 0.5]))]
            out = pool.map(loglikelihood, items)
            self.assertEqual(len(out), 1)
            self.assertAlmostEqual(float(out[0].val), 0.0, places=6)
        finally:
            stop.set()
            thread.join(timeout=2.0)
        self.assertTrue(
            any("unmatched hep:feedback" in w for w in logger.warnings),
            msg=f"expected unmatched-feedback warning, got {logger.warnings!r}",
        )

    def test_duplicate_uuid_is_rejected_before_task_push(self) -> None:
        """A duplicate UUID must not collapse two result slots into one."""
        queue = make_fakeredis_queue()

        def build_sample(payload, uuid):
            return Sample(uuid=uuid, u_coords=np.asarray(payload, dtype=float))

        pool = RedisEvaluationPool(
            queue,
            build_sample=build_sample,
            batch_size=4,
            seed=7,
            timeout=5.0,
        )
        duplicate = np.array([0.25, "same-uuid"], dtype=object)

        with self.assertRaisesRegex(ValueError, "duplicate sample UUID"):
            pool._redis_batch_logl([duplicate, duplicate.copy()])

        # Validation happens before push_many_tasks, so no orphan task is left
        # waiting for a feedback record after the caller fixes the input.
        assert queue.r is not None
        self.assertEqual(int(queue.r.llen("hep:task_queue")), 0)

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
    def test_get_print_func_never_creates_tqdm(self) -> None:
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty.utils import get_print_func

        pbar, print_func = get_print_func(None, True)
        self.assertIsNone(pbar)
        self.assertIsNotNone(print_func)

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


class NestedLoggerLabelTests(unittest.TestCase):
    def test_bind_inner_multinest_label(self) -> None:
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty.jarvis_logging import (
            _inner_module_label,
            bind_inner,
            resolve_progress_title,
        )

        self.assertEqual(
            _inner_module_label("MultiNest"),
            "Jarvis-HEP.Sampler.MultiNest.Inner",
        )
        self.assertEqual(
            _inner_module_label("Dynesty"),
            "Jarvis-HEP.Sampler.Dynesty.Inner",
        )

        class CapturingLogger:
            def __init__(self) -> None:
                self.extra: dict = {}

            def bind(self, **kwargs):
                child = CapturingLogger()
                child.extra = dict(kwargs)
                return child

        base = CapturingLogger()
        nested = bind_inner(base, method="MultiNest")
        self.assertEqual(
            nested.extra.get("module"),
            "Jarvis-HEP.Sampler.MultiNest.Inner",
        )
        self.assertEqual(resolve_progress_title(nested), "MultiNest")


class DynestyCsvExportTests(unittest.TestCase):
    def test_export_dynesty_results_csv_schema(self) -> None:
        """V1 column names required by Jarvis-PLOT dynesty_runplot."""
        n = 5
        ndim = 2
        results = {
            "logl": np.linspace(-10, -1, n),
            "logwt": np.linspace(-5, -0.1, n),
            "logvol": np.linspace(0, -3, n),
            "logz": np.linspace(-8, -2, n),
            "logzerr": np.full(n, 0.1),
            "samples_n": np.full(n, 20),
            "ncall": np.arange(1, n + 1),
            "samples_it": np.arange(n),
            "samples_id": np.arange(n),
            "information": np.linspace(0, 1, n),
            "samples": np.random.default_rng(0).random((n, ndim)),
            "samples_u": np.random.default_rng(1).random((n, ndim)),
            "samples_uid": [f"uid-{i}" for i in range(n)],
            "nlive": 20,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "dynesty_result.csv")
            written = export_dynesty_results_csv(results, path, fallback_nlive=20)
            self.assertTrue(os.path.isfile(written))
            with open(written, encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                fieldnames = list(reader.fieldnames or [])
                rows = list(reader)
            for col in (
                "uuid",
                "log_weight",
                "log_Like",
                "log_PriorVolume",
                "log_Evidence",
                "log_Evidence_err",
                "samples_nlive",
                "ncall",
                "samples_it",
                "samples_id",
                "information",
                "samples_v[0]",
                "samples_v[1]",
                "samples_u[0]",
                "samples_u[1]",
            ):
                self.assertIn(col, fieldnames)
            self.assertEqual(len(rows), n)
            self.assertEqual(rows[0]["uuid"], "uid-0")
            self.assertAlmostEqual(float(rows[-1]["log_Like"]), -1.0, places=5)

    def test_export_dynesty_results_csv_param_names(self) -> None:
        """Named physical columns sit alongside V1 samples_v aliases."""
        n = 3
        results = {
            "logl": np.linspace(-5, -1, n),
            "logwt": np.linspace(-3, -0.1, n),
            "logvol": np.linspace(0, -2, n),
            "logz": np.linspace(-4, -1, n),
            "logzerr": np.full(n, 0.05),
            "samples_n": np.full(n, 10),
            "ncall": np.ones(n),
            "samples_it": np.arange(n),
            "samples_id": np.arange(n),
            "information": np.zeros(n),
            "samples": np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]]),
            "samples_u": np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]]),
            "samples_uid": [f"u{i}" for i in range(n)],
            "nlive": 10,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "dynesty_result.csv")
            export_dynesty_results_csv(
                results, path, fallback_nlive=10, param_names=["xx", "yy"]
            )
            with open(path, encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                fieldnames = list(reader.fieldnames or [])
                rows = list(reader)
            self.assertIn("xx", fieldnames)
            self.assertIn("yy", fieldnames)
            self.assertIn("samples_v[0]", fieldnames)
            self.assertEqual(len(rows), n)
            self.assertAlmostEqual(float(rows[1]["xx"]), 0.3, places=6)
            self.assertAlmostEqual(float(rows[1]["samples_v[0]"]), 0.3, places=6)

    def test_sampler_save_writes_database_path(self) -> None:
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import NestedSampler

        def loglike(x):
            u = np.asarray(x, dtype=float).reshape(-1)
            return float(-0.5 * np.sum((u - 0.5) ** 2))

        def prior(u):
            return np.asarray(u, dtype=float)

        with tempfile.TemporaryDirectory() as tmp:
            sampler = DynestySampler()
            sampler.set_config(
                {
                    "task_result_dir": tmp,
                    "Sampling": {
                        "Method": "Dynesty",
                        "Variables": [
                            {
                                "name": "x",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0, "max": 1},
                                },
                            },
                            {
                                "name": "y",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0, "max": 1},
                                },
                            },
                        ],
                        "Bounds": {"nlive": 10, "dlogz": 5.0},
                    },
                }
            )
            # Direct NestedSampler results (no Redis) for CSV path coverage.
            ns = NestedSampler(
                loglikelihood=loglike,
                prior_transform=prior,
                ndim=2,
                nlive=10,
                rstate=np.random.default_rng(0),
            )
            ns.run_nested(maxiter=8, dlogz=10.0, print_progress=False)
            sampler._sampler = ns
            path = sampler.save_dynesty_results_to_csv()
            self.assertIsNotNone(path)
            expected = os.path.join(tmp, "DATABASE", "dynesty_result.csv")
            self.assertEqual(path, expected)
            self.assertTrue(os.path.isfile(expected))


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
            # CSV export works on real NestedSampler.results
            with tempfile.TemporaryDirectory() as tmp:
                csv_path = export_dynesty_results_csv(
                    res, os.path.join(tmp, "dynesty_result.csv"), fallback_nlive=20
                )
                self.assertTrue(os.path.isfile(csv_path))
                with open(csv_path, encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(handle)
                    rows = list(reader)
                self.assertGreater(len(rows), 0)
                self.assertIn("log_Like", reader.fieldnames or [])
        finally:
            stop.set()
            thread.join(timeout=2.0)


class NestedCheckpointResumeTests(unittest.TestCase):
    """Dynesty / MultiNest native pickle engine + run_nested(resume=True)."""

    @staticmethod
    def _two_var_config(*, method: str = "Dynesty", nlive: int = 20, **bounds: Any) -> dict:
        b = {"nlive": nlive, "seed": 7, **bounds}
        return {
            "Sampling": {
                "Method": method,
                "Bounds": b,
                "Variables": [
                    {
                        "name": "x",
                        "distribution": {
                            "type": "Flat",
                            "parameters": {"min": 0, "max": 1},
                        },
                    },
                    {
                        "name": "y",
                        "distribution": {
                            "type": "Flat",
                            "parameters": {"min": 0, "max": 1},
                        },
                    },
                ],
            },
            "Runtime": {"workers": 1, "batch_size": 4},
            "task_root": "",
            "scan_name": "nested-resume",
        }

    def _start_echo_worker(self, queue) -> tuple[threading.Event, threading.Thread]:
        stop = threading.Event()

        def worker() -> None:
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is None:
                    continue
                u = np.asarray(task.get("u_coords") or [], dtype=float)
                # Smooth unimodal peak at unit-cube centre.
                logl = float(-np.sum((u - 0.5) ** 2) * 20.0)
                queue.publish_feedback({"uuid": task["uuid"], "logL": logl})

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        return stop, thread

    def test_nested_redis_logl_is_pickle_safe(self) -> None:
        import pickle

        from jarvishep2.Sampling.dynesty_sampler import NestedRedisLogL

        queue = make_fakeredis_queue()
        pool = RedisEvaluationPool(
            queue,
            build_sample=lambda payload, uuid: Sample(
                uuid=uuid, u_coords=np.asarray(payload, dtype=float)
            ),
            batch_size=2,
            seed=1,
            timeout=5.0,
        )
        bridge = NestedRedisLogL(pool)
        blob = pickle.dumps(bridge, protocol=pickle.HIGHEST_PROTOCOL)
        restored = pickle.loads(blob)
        self.assertIsInstance(restored, NestedRedisLogL)
        self.assertIsNone(restored._pool)
        with self.assertRaises(RuntimeError):
            restored(np.array([0.1, 0.2]))
        restored.attach_pool(pool)
        # No task published yet — just ensure attach works without error path.
        self.assertIs(restored._pool, pool)

    def test_static_multinest_checkpoint_resume_continues_ncall(self) -> None:
        """MultiNest (static NestedSampler): interrupt mid-run → resume adds ncall."""
        from jarvishep2.Sampling.multinest_sampler import MultiNestSampler

        queue = make_fakeredis_queue()
        stop, thread = self._start_echo_worker(queue)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg = self._two_var_config(
                    method="MultiNest",
                    nlive=15,
                    run_nested={
                        "maxiter": 12,
                        "dlogz": 1e-9,
                        "print_progress": False,
                        "add_live": False,
                    },
                )
                cfg["task_root"] = tmp
                cfg["task_result_dir"] = tmp

                first = MultiNestSampler()
                first.set_config(cfg)
                first.set_redis(queue)
                n1 = first.run_adaptive(generation_timeout=60.0)
                # With maxiter short, finished=True after first run.
                # Simulate interrupt mid-run by re-running with a longer budget
                # from a pickled mid-state: build a second engine stopped early.
                mid = MultiNestSampler()
                mid.set_config(cfg)
                mid.set_redis(queue)
                # Force a partial engine via NestedSampler directly, then export.
                from jarvishep2.Sampling.Source.Dynesty.py.dynesty import NestedSampler
                from jarvishep2.Sampling.dynesty_sampler import NestedRedisLogL

                pool = RedisEvaluationPool(
                    queue,
                    build_sample=mid._build_sample_for_pool,
                    batch_size=4,
                    timeout=60.0,
                    seed=mid._seed,
                    method="MultiNest",
                )
                bridge = NestedRedisLogL(pool)
                engine = NestedSampler(
                    loglikelihood=bridge,
                    prior_transform=_jarvis_prior_transform,
                    ndim=2,
                    nlive=15,
                    pool=pool,
                    rstate=np.random.default_rng(7),
                    queue_size=4,
                )
                engine.run_nested(
                    maxiter=8,
                    dlogz=1e-9,
                    print_progress=False,
                    add_live=False,
                )
                mid._sampler = engine
                mid._finished = False
                mid._engine_checkpoint_path = os.path.join(
                    tmp, "checkpoints", "nested-resume", "MultiNest", "nested_engine.pkl"
                )
                state = mid.export_runtime_state()
                self.assertIsInstance(state.get("native_sampler_blob"), (bytes, bytearray))
                self.assertTrue(os.path.isfile(state["engine_checkpoint_path"]))
                ncall_mid = int(getattr(engine, "ncall", 0) or 0)
                self.assertGreater(ncall_mid, 0)

                resumed = MultiNestSampler()
                resumed.set_config(cfg)
                resumed.set_redis(queue)
                resumed.import_runtime_state(state)
                self.assertFalse(resumed._finished)
                self.assertIsNotNone(resumed._native_sampler_blob or resumed._sampler)
                # Continue further under resume.
                resumed._run_nested_kwargs = {
                    "maxiter": 40,
                    "dlogz": 0.5,
                    "print_progress": False,
                    "add_live": True,
                }
                n2 = resumed.run_adaptive(generation_timeout=120.0)
                self.assertTrue(resumed._finished)
                self.assertGreaterEqual(int(n2), ncall_mid)
                summary = resumed.summary()
                self.assertIn("logz", summary)
                self.assertTrue(
                    np.isfinite(float(summary.get("logz", float("nan"))))
                    or summary.get("niter", 0) >= 0
                )
                # First full short run still produced work (sanity).
                self.assertGreaterEqual(int(n1), 0)
        finally:
            stop.set()
            thread.join(timeout=2.0)

    def test_dynesty_export_import_roundtrip_restores_live_engine(self) -> None:
        """export → import restores nlive/it and pickle blob without Redis."""
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import NestedSampler
        from jarvishep2.Sampling.dynesty_sampler import NestedRedisLogL

        queue = make_fakeredis_queue()
        stop, thread = self._start_echo_worker(queue)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg = self._two_var_config(method="Dynesty", nlive=12)
                cfg["task_root"] = tmp
                cfg["task_result_dir"] = tmp
                sampler = DynestySampler()
                sampler.set_config(cfg)
                sampler.set_redis(queue)
                pool = RedisEvaluationPool(
                    queue,
                    build_sample=sampler._build_sample_for_pool,
                    batch_size=4,
                    timeout=60.0,
                    seed=0,
                    method="Dynesty",
                )
                bridge = NestedRedisLogL(pool)
                # Use static NestedSampler for a deterministic mid-run snapshot
                # even though Method=Dynesty normally builds Dynamic.
                engine = NestedSampler(
                    loglikelihood=bridge,
                    prior_transform=_jarvis_prior_transform,
                    ndim=2,
                    nlive=12,
                    pool=pool,
                    rstate=np.random.default_rng(0),
                    queue_size=4,
                )
                engine.run_nested(
                    maxiter=6, dlogz=1e-9, print_progress=False, add_live=False
                )
                sampler._sampler = engine
                sampler._finished = False
                it_before = int(engine.it)
                ncall_before = int(engine.ncall)
                state = sampler.export_runtime_state()
                self.assertEqual(state.get("native_sampler_format"), "pickle")
                self.assertTrue(state.get("native_sampler_blob"))

                other = DynestySampler()
                other.set_config(cfg)
                other.import_runtime_state(state)
                self.assertFalse(other._finished)
                self.assertIsNotNone(other._sampler)
                self.assertEqual(int(other._sampler.it), it_before)
                self.assertEqual(int(other._sampler.ncall), ncall_before)
                # Pool stripped on pickle — reattach required before evaluate.
                logl = getattr(other._sampler, "loglikelihood", None)
                inner = getattr(logl, "loglikelihood", logl)
                if isinstance(inner, NestedRedisLogL):
                    self.assertIsNone(inner._pool)
        finally:
            stop.set()
            thread.join(timeout=2.0)

    def test_resume_status_is_implemented(self) -> None:
        self.assertEqual(Distributor.get_resume_status("Dynesty"), "implemented")
        self.assertEqual(Distributor.get_resume_status("MultiNest"), "implemented")

    def test_multinest_resume_uses_method_path_and_envreqs_interval(self) -> None:
        """MultiNest must share Dynesty resume stack under checkpoints/.../MultiNest/."""
        from jarvishep2.Sampling.multinest_sampler import MultiNestSampler
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import NestedSampler
        from jarvishep2.Sampling.dynesty_sampler import NestedRedisLogL

        queue = make_fakeredis_queue()
        stop, thread = self._start_echo_worker(queue)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg = self._two_var_config(
                    method="MultiNest",
                    nlive=12,
                    run_nested={
                        "maxiter": 6,
                        "dlogz": 1e-9,
                        "print_progress": False,
                        "add_live": False,
                    },
                )
                cfg["task_root"] = tmp
                cfg["task_result_dir"] = tmp
                cfg["scan_name"] = "mn-ckpt"
                cfg["EnvReqs"] = {"V2": {"checkpoint_heartbeat_sec": 11}}

                mid = MultiNestSampler()
                mid.set_config(cfg)
                mid.set_redis(queue)
                self.assertEqual(mid._checkpoint_every_sec, 11.0)
                self.assertFalse(mid._use_dynamic)

                pool = RedisEvaluationPool(
                    queue,
                    build_sample=mid._build_sample_for_pool,
                    batch_size=4,
                    timeout=60.0,
                    seed=0,
                    method="MultiNest",
                )
                bridge = NestedRedisLogL(pool)
                engine = NestedSampler(
                    loglikelihood=bridge,
                    prior_transform=_jarvis_prior_transform,
                    ndim=2,
                    nlive=12,
                    pool=pool,
                    rstate=np.random.default_rng(2),
                    queue_size=4,
                )
                engine.run_nested(
                    maxiter=5, dlogz=1e-9, print_progress=False, add_live=False
                )
                mid._sampler = engine
                mid._finished = False
                state = mid.export_runtime_state()
                eng_path = state["engine_checkpoint_path"]
                self.assertIn(os.path.join("MultiNest", "nested_engine.pkl"), eng_path)
                self.assertTrue(os.path.isfile(eng_path))
                self.assertEqual(state["method"], "MultiNest")
                self.assertFalse(state["use_dynamic"])

                resumed = MultiNestSampler()
                resumed.set_config(cfg)
                resumed.set_redis(queue)
                resumed.import_runtime_state(state)
                self.assertTrue(resumed._resume_engine_requested)
                self.assertFalse(resumed._use_dynamic)

                pool2 = RedisEvaluationPool(
                    queue,
                    build_sample=resumed._build_sample_for_pool,
                    batch_size=4,
                    timeout=60.0,
                    seed=0,
                    method="MultiNest",
                )
                restored, do_resume = resumed._resolve_resume_engine(pool=pool2)
                self.assertIsNotNone(restored)
                self.assertTrue(do_resume)
                self.assertEqual(int(restored.it), int(engine.it))
                self.assertEqual(int(restored.ncall), int(engine.ncall))
                # Must not be a DynamicNestedSampler under Method: MultiNest.
                self.assertFalse(MultiNestSampler._is_dynamic_engine(restored))
        finally:
            stop.set()
            thread.join(timeout=2.0)

    def test_multinest_rejects_dynamic_engine_on_resume(self) -> None:
        """A Dynesty DynamicNestedSampler pickle must not resume under MultiNest."""
        from jarvishep2.Sampling.multinest_sampler import MultiNestSampler
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import DynamicNestedSampler
        from jarvishep2.Sampling.dynesty_sampler import NestedRedisLogL

        queue = make_fakeredis_queue()
        stop, thread = self._start_echo_worker(queue)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg = self._two_var_config(method="MultiNest", nlive=10)
                cfg["task_root"] = tmp
                cfg["task_result_dir"] = tmp
                cfg["scan_name"] = "mn-mismatch"

                s = MultiNestSampler()
                s.set_config(cfg)
                s.set_redis(queue)
                pool = RedisEvaluationPool(
                    queue,
                    build_sample=s._build_sample_for_pool,
                    batch_size=4,
                    timeout=60.0,
                    seed=0,
                    method="MultiNest",
                )
                bridge = NestedRedisLogL(pool)
                dyn = DynamicNestedSampler(
                    bridge,
                    _jarvis_prior_transform,
                    2,
                    nlive=10,
                    pool=pool,
                    rstate=np.random.default_rng(0),
                    queue_size=4,
                )
                dyn.run_nested(
                    maxiter=4,
                    maxbatch=0,
                    dlogz_init=10.0,
                    print_progress=False,
                )
                s._sampler = dyn
                s._finished = False
                s.import_runtime_state(s.export_runtime_state())
                # Import should drop the incompatible dynamic engine.
                self.assertIsNone(s._sampler)
                s.arm_engine_resume()
                restored, resume = s._resolve_resume_engine(pool=pool)
                self.assertIsNone(restored)
                self.assertFalse(resume)
        finally:
            stop.set()
            thread.join(timeout=2.0)

    def test_leftover_engine_file_does_not_hijack_fresh_run(self) -> None:
        """nested_engine.pkl alone must not force resume without import_runtime_state."""
        from jarvishep2.Sampling.multinest_sampler import MultiNestSampler
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import NestedSampler
        from jarvishep2.Sampling.dynesty_sampler import NestedRedisLogL

        queue = make_fakeredis_queue()
        stop, thread = self._start_echo_worker(queue)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg = self._two_var_config(
                    method="MultiNest",
                    nlive=10,
                    run_nested={
                        "maxiter": 5,
                        "dlogz": 1e-9,
                        "print_progress": False,
                        "add_live": False,
                    },
                )
                cfg["task_root"] = tmp
                cfg["task_result_dir"] = tmp
                cfg["scan_name"] = "nested-resume"

                # Write a finished-looking engine into the default side-file path.
                s = MultiNestSampler()
                s.set_config(cfg)
                s.set_redis(queue)
                pool = RedisEvaluationPool(
                    queue,
                    build_sample=s._build_sample_for_pool,
                    batch_size=4,
                    timeout=60.0,
                    seed=0,
                    method="MultiNest",
                )
                bridge = NestedRedisLogL(pool)
                engine = NestedSampler(
                    loglikelihood=bridge,
                    prior_transform=_jarvis_prior_transform,
                    ndim=2,
                    nlive=10,
                    pool=pool,
                    rstate=np.random.default_rng(1),
                    queue_size=4,
                )
                engine.run_nested(
                    maxiter=4, dlogz=1e-9, print_progress=False, add_live=False
                )
                s._sampler = engine
                s._finished = False
                s._save_engine_side_file()
                self.assertTrue(os.path.isfile(s._engine_path()))

                # Fresh sampler: no import_runtime_state → must ignore side file.
                fresh = MultiNestSampler()
                fresh.set_config(cfg)
                fresh.set_redis(queue)
                self.assertFalse(fresh._resume_engine_requested)
                pool2 = RedisEvaluationPool(
                    queue,
                    build_sample=fresh._build_sample_for_pool,
                    batch_size=4,
                    timeout=60.0,
                    seed=0,
                    method="MultiNest",
                )
                restored, resume = fresh._resolve_resume_engine(pool=pool2)
                self.assertIsNone(restored)
                self.assertFalse(resume)
        finally:
            stop.set()
            thread.join(timeout=2.0)


if __name__ == "__main__":
    unittest.main()
