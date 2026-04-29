#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import threading
import time
import unittest

import numpy as np


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Sampling.Source.Dynesty.py.dynesty.pool import JarvisFactoryAsyncPool  # noqa: E402
from jarvishep.Sampling.Source.Dynesty.py.dynesty.sampler import Sampler  # noqa: E402


class _GuardedIterable:
    """Raise if consumer tries to over-read before gate is opened."""

    def __init__(self, total: int, prefetch_limit: int, gate: threading.Event):
        self.total = int(total)
        self.prefetch_limit = int(prefetch_limit)
        self.gate = gate
        self.idx = 0
        self.next_calls = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.idx >= self.total:
            raise StopIteration
        self.next_calls += 1
        if self.next_calls > self.prefetch_limit and not self.gate.is_set():
            raise RuntimeError("eager iterable consumption detected")
        value = self.idx
        self.idx += 1
        return value


class _AsyncQueueProbeSampler(Sampler):
    def __init__(self, pool, delays):
        self._delays = [float(x) for x in delays]
        self._next_idx = 0
        live_u = np.zeros((2, 1), dtype=float)
        live_v = np.zeros((2, 1), dtype=float)
        live_uid = np.array(["u0", "u1"], dtype=object)
        live_logl = np.array([1.0, 1.0], dtype=float)
        super().__init__(
            loglikelihood=lambda x: 0.0,
            prior_transform=lambda x: x,
            npdim=1,
            live_points=(live_u, live_v, live_uid, live_logl),
            update_interval=1,
            first_update={},
            rstate=np.random.default_rng(0),
            queue_size=2,
            pool=pool,
            use_pool={
                "prior_transform": False,
                "loglikelihood": False,
                "propose_point": True,
                "update_bound": False,
            },
            ncdim=1,
        )
        self.method = "rwalk"
        self.unit_cube_sampling = False
        self.scale = 1.0

    def save(self, fname):
        return None

    def propose_point(self, *args):
        idx = self._next_idx
        self._next_idx += 1
        return np.array([idx], dtype=float), np.identity(self.ncdim)

    def evolve_point(self, sampler_arg):
        idx = int(sampler_arg.u[0])
        time.sleep(self._delays[idx])
        return (
            sampler_arg.u,
            np.array([idx, f"uid-{idx}"], dtype=object),
            float(idx),
            1,
            None,
        )

    def update_proposal(self, *args, **kwargs):
        return None

    def update(self, subset=None):
        return None


class DynestyPoolAdapterTests(unittest.TestCase):
    def test_map_streams_iterable_without_eager_full_materialization(self):
        gate = threading.Event()
        pool = JarvisFactoryAsyncPool(njobs=2)
        data = _GuardedIterable(total=7, prefetch_limit=2, gate=gate)

        def _worker(x):
            time.sleep(0.01)
            gate.set()
            return x * x

        try:
            out = pool.map(_worker, data)
        finally:
            pool.shutdown(wait_for_tasks=True, cancel_futures=True)

        self.assertEqual(out, [0, 1, 4, 9, 16, 25, 36])
        self.assertEqual(data.next_calls, 7)

    def test_map_raises_after_shutdown(self):
        pool = JarvisFactoryAsyncPool(njobs=1)
        pool.shutdown(wait_for_tasks=True, cancel_futures=True)
        with self.assertRaises(RuntimeError):
            pool.map(lambda x: x, [1, 2, 3])

    def test_sampler_queue_returns_first_completed_without_batch_barrier(self):
        pool = JarvisFactoryAsyncPool(njobs=2)
        sampler = _AsyncQueueProbeSampler(pool=pool, delays=[0.25, 0.02])
        try:
            t0 = time.monotonic()
            _, vid, _, _, _ = sampler._get_point_value(loglstar=0.0)
            elapsed = time.monotonic() - t0
        finally:
            pool.shutdown(wait_for_tasks=True, cancel_futures=True)

        self.assertLess(elapsed, 0.2)
        self.assertEqual(int(vid[0]), 1)
        self.assertEqual(sampler.nqueue, 1)


if __name__ == "__main__":
    unittest.main()
