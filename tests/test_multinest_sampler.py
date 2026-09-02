#!/usr/bin/env python3
"""MultiNest (static NestedSampler) + plot CSV path tests."""

from __future__ import annotations

import csv
import os
import tempfile
import threading
import time
import unittest
from typing import Any

import numpy as np

from jarvishep2.Sampling.dynesty_sampler import (
    NestedRedisLogL,
    _jarvis_prior_transform,
    export_dynesty_results_csv,
)
from jarvishep2.Sampling.multinest_sampler import MultiNestSampler, create_multinest
from jarvishep2.Sampling.redis_evaluation_pool import RedisEvaluationPool
from jarvishep2.distributor import STATELESS_METHODS, Distributor
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.sample import Sample


class MultiNestRegistrationTests(unittest.TestCase):
    def test_registered_as_feedback(self) -> None:
        self.assertNotIn("MultiNest", STATELESS_METHODS)
        sampler = Distributor.set_method("MultiNest")
        self.assertIsInstance(sampler, MultiNestSampler)
        self.assertEqual(sampler.method, "MultiNest")

    def test_create_factory(self) -> None:
        self.assertIsInstance(create_multinest(), MultiNestSampler)


class MultiNestConfigTests(unittest.TestCase):
    def test_always_static(self) -> None:
        sampler = MultiNestSampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "MultiNest",
                    "Variables": [
                        {
                            "name": "x",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1},
                            },
                        },
                    ],
                    "Bounds": {
                        "nlive": 20,
                        "dlogz": 1.0,
                        "seed": 3,
                    },
                },
            }
        )
        self.assertFalse(sampler._use_dynamic)
        self.assertEqual(sampler._nlive, 20)
        self.assertEqual(sampler._dim, 1)


class MultiNestCsvExportTests(unittest.TestCase):
    def test_save_writes_multinest_database_path(self) -> None:
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import NestedSampler

        def loglike(x):
            u = np.asarray(x, dtype=float).reshape(-1)
            return float(-0.5 * np.sum((u - 0.5) ** 2))

        def prior(u):
            return np.asarray(u, dtype=float)

        with tempfile.TemporaryDirectory() as tmp:
            sampler = MultiNestSampler()
            sampler.set_config(
                {
                    "task_result_dir": tmp,
                    "Sampling": {
                        "Method": "MultiNest",
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
            ns = NestedSampler(
                loglikelihood=loglike,
                prior_transform=prior,
                ndim=2,
                nlive=10,
                rstate=np.random.default_rng(1),
            )
            ns.run_nested(maxiter=8, dlogz=10.0, print_progress=False)
            sampler._sampler = ns
            path = sampler.save_multinest_results_to_csv()
            expected = os.path.join(tmp, "DATABASE", "multinest_result.csv")
            self.assertEqual(path, expected)
            self.assertTrue(os.path.isfile(expected))
            # Must not write the Dynesty filename by default
            self.assertFalse(
                os.path.isfile(os.path.join(tmp, "DATABASE", "dynesty_result.csv"))
            )
            with open(expected, encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                fieldnames = list(reader.fieldnames or [])
                rows = list(reader)
            self.assertGreater(len(rows), 0)
            self.assertIn("log_Like", fieldnames)
            self.assertIn("log_Evidence", fieldnames)
            # Dynesty parity: physical Sampling.Variable names must be present.
            self.assertIn("x", fieldnames)
            self.assertIn("y", fieldnames)
            self.assertIn("samples_v[0]", fieldnames)

    def test_export_schema_matches_dynesty_plot_columns(self) -> None:
        n = 3
        results = {
            "logl": np.array([-2.0, -1.0, -0.5]),
            "logwt": np.array([-1.0, -0.5, -0.1]),
            "logvol": np.array([0.0, -0.5, -1.0]),
            "logz": np.array([-3.0, -2.5, -2.0]),
            "logzerr": np.full(n, 0.1),
            "samples_n": np.full(n, 15),
            "ncall": np.arange(1, n + 1),
            "samples_it": np.arange(n),
            "samples_id": np.arange(n),
            "information": np.linspace(0, 1, n),
            "samples": np.zeros((n, 1)),
            "samples_u": np.zeros((n, 1)),
            "samples_uid": ["a", "b", "c"],
            "nlive": 15,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = export_dynesty_results_csv(
                results, os.path.join(tmp, "multinest_result.csv")
            )
            with open(path, encoding="utf-8", newline="") as handle:
                fieldnames = list(csv.DictReader(handle).fieldnames or [])
        for col in (
            "uuid",
            "log_weight",
            "log_Like",
            "log_PriorVolume",
            "log_Evidence",
            "log_Evidence_err",
            "samples_nlive",
        ):
            self.assertIn(col, fieldnames)


class MultiNestParallelEvolveTests(unittest.TestCase):
    """Static NestedSampler uses the same Redis evolve path as Dynesty."""

    def test_constructor_sets_queue_size_from_batch_size(self) -> None:
        sampler = MultiNestSampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "MultiNest",
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
                    "Bounds": {"nlive": 16, "dlogz": 1.0, "seed": 3},
                },
                "EnvReqs": {"V2": {"workers": 4, "batch_size": 4}},
            }
        )
        self.assertFalse(sampler._use_dynamic)
        self.assertEqual(sampler._batch_size, 4)
        queue = make_fakeredis_queue()
        pool = RedisEvaluationPool(
            queue,
            build_sample=lambda payload, uuid: Sample(
                uuid=uuid, u_coords=np.asarray(payload, dtype=float)
            ),
            batch_size=sampler._batch_size,
            method="MultiNest",
        )
        kwargs = sampler._build_constructor_kwargs(
            pool=pool, rstate=np.random.default_rng(0)
        )
        self.assertEqual(kwargs["queue_size"], 4)
        self.assertIs(kwargs["pool"], pool)
        self.assertIsInstance(kwargs["loglikelihood"], NestedRedisLogL)

    def test_static_nested_sampler_evolve_runs_in_parallel(self) -> None:
        """NestedSampler._fill_queue must overlap Redis logL across queue_size."""
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty import NestedSampler

        queue = make_fakeredis_queue()
        logs: list[str] = []

        class _CapturingLogger:
            def info(self, msg: str, *args: Any, **kwargs: Any) -> None:
                logs.append(msg % args if args else str(msg))

            def warning(self, msg: str, *args: Any, **kwargs: Any) -> None:
                logs.append(msg % args if args else str(msg))

        pool = RedisEvaluationPool(
            queue,
            build_sample=lambda payload, uuid: Sample(
                uuid=uuid, u_coords=np.asarray(payload, dtype=float)
            ),
            batch_size=4,
            seed=5,
            timeout=30.0,
            method="MultiNest",
            logger=_CapturingLogger(),
        )
        stop = threading.Event()
        delay = 0.2

        def worker() -> None:
            while not stop.is_set():
                task = queue.pull_task(timeout=1)
                if task is None:
                    continue
                time.sleep(delay)
                u = np.asarray(task.get("u_coords") or [], dtype=float)
                queue.publish_feedback(
                    {"uuid": task["uuid"], "logL": float(-np.sum((u - 0.5) ** 2))}
                )

        workers = [threading.Thread(target=worker, daemon=True) for _ in range(4)]
        for thread in workers:
            thread.start()
        evolve_sizes: list[int] = []
        orig = pool._map_sampler_arguments

        def _spy(func, items):
            evolve_sizes.append(len(items))
            return orig(func, items)

        pool._map_sampler_arguments = _spy  # type: ignore[method-assign]
        try:
            t0 = time.monotonic()
            sampler = NestedSampler(
                loglikelihood=NestedRedisLogL(pool),
                prior_transform=_jarvis_prior_transform,
                ndim=2,
                nlive=8,
                pool=pool,
                queue_size=4,
                rstate=np.random.default_rng(0),
            )
            sampler.run_nested(
                maxiter=4, dlogz=1e9, print_progress=False, add_live=False
            )
            elapsed = time.monotonic() - t0
            self.assertGreaterEqual(int(getattr(sampler, "ncall", 0) or 0), 8)
            self.assertTrue(evolve_sizes, msg="NestedSampler never mapped SamplerArgument")
            self.assertGreaterEqual(max(evolve_sizes), 2)
            self.assertTrue(
                any("MultiNest evolve" in line for line in logs),
                msg=f"expected MultiNest evolve log, got {logs!r}",
            )
            # Serial live-init + 4 queue fills ≈ 8*0.2 + 4*4*0.2 = 4.8s.
            self.assertLess(
                elapsed, 2.5, msg=f"MultiNest NestedSampler looked serial: {elapsed:.3f}s"
            )
        finally:
            stop.set()
            for thread in workers:
                thread.join(timeout=2.0)
            pool.close()


if __name__ == "__main__":
    unittest.main()
