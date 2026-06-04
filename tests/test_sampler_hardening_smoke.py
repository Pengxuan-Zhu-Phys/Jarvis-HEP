#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import json
import csv
import os
import sys
import tempfile
import time
import threading
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import sympy as sp


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Sampling.bridson import Bridson  # noqa: E402
from jarvishep.Sampling.diver import Diver  # noqa: E402
from jarvishep.Sampling.dynesty import Dynesty  # noqa: E402
from jarvishep.Sampling.grid import Grid  # noqa: E402
from jarvishep.Sampling.ammcmc import AMMCMC  # noqa: E402
from jarvishep.Sampling.demcmc import DEMCMC  # noqa: E402
from jarvishep.Sampling.dream import DREAM  # noqa: E402
from jarvishep.Sampling.dream_lite import DREAMLite  # noqa: E402
from jarvishep.Sampling.dram import DRAM  # noqa: E402
from jarvishep.Sampling.ess import ESS  # noqa: E402
from jarvishep.Sampling.ensemblemcmc import EnsembleMCMC  # noqa: E402
from jarvishep.Sampling.mcmc_standard import MCMC  # noqa: E402
from jarvishep.Sampling.multinest import MultiNest  # noqa: E402
from jarvishep.Sampling.hmc import HMC  # noqa: E402
from jarvishep.Sampling.mala import MALA  # noqa: E402
from jarvishep.Sampling.nuts import NUTS  # noqa: E402
from jarvishep.Sampling.pt_ensemble import PTEnsemble  # noqa: E402
from jarvishep.Sampling.randoms import RandomS  # noqa: E402
from jarvishep.Sampling.robustam import RobustAM  # noqa: E402
from jarvishep.Sampling.slicemcmc import SliceMCMC  # noqa: E402
from jarvishep.Sampling.tpmcmc import PTMCMC  # noqa: E402
from jarvishep.distributor import Distributor  # noqa: E402
from jarvishep.moduleManager import ModuleManager  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _CaptureLogger(_NoopLogger):
    def __init__(self):
        self.records = []

    def info(self, *args, **kwargs):
        self.records.append(("INFO", args, kwargs))
        return None

    def warning(self, *args, **kwargs):
        self.records.append(("WARNING", args, kwargs))
        return None

    def error(self, *args, **kwargs):
        self.records.append(("ERROR", args, kwargs))
        return None


class _FakeSample:
    seq = 0
    close_calls = 0

    @classmethod
    def reset(cls):
        cls.seq = 0
        cls.close_calls = 0

    def __init__(self, params):
        type(self).seq += 1
        self.uuid = f"fake-sample-{type(self).seq:03d}"
        self.params = dict(params)
        self.info = {
            "uuid": self.uuid,
            "params": dict(self.params),
            "observables": {**self.params, "uuid": self.uuid},
            "status": "Accept",
        }

    def set_config(self, config):
        if config:
            self.info.update(dict(config))
        self.info["uuid"] = self.uuid
        self.info["params"] = dict(self.params)
        observables = dict(self.info.get("observables", {}))
        observables.update(self.params)
        observables["uuid"] = self.uuid
        self.info["observables"] = observables

    def update_uuid(self, uuid):
        self.uuid = str(uuid)
        self.info["uuid"] = self.uuid
        self.info["observables"]["uuid"] = self.uuid

    def evaluate_output(self, outputs):
        observables = self.info.get("observables", {})
        return {name: observables.get(name, 0.0) for name in outputs}

    def start(self):
        self.info["status"] = "Running"

    def record(self):
        return None

    def combine_nuisance_card(self):
        self.info["status"] = "Accept"

    def close(self):
        type(self).close_calls += 1


class _CaptureConfigSample(_FakeSample):
    save_dirs = []

    @classmethod
    def reset(cls):
        super().reset()
        cls.save_dirs = []

    def set_config(self, config):
        if isinstance(config, dict):
            sd = config.get("save_dir")
            if sd:
                type(self).save_dirs.append(str(sd))
        return super().set_config(config)


class _ImmediateFactory:
    def __init__(self):
        self.calls = 0

    def submit_task(self, sample_info):
        self.calls += 1
        future = concurrent.futures.Future()
        observables = sample_info.setdefault("observables", {})
        params = sample_info.get("params", {})
        value = float(self.calls)
        if isinstance(params, dict) and params:
            value = float(sum(float(v) for v in params.values()))
        observables.setdefault("obs1", value + 1.0)
        observables["LogL"] = value
        future.set_result(value)
        return future


class _BlockingFactory:
    def __init__(self, max_workers=16, delay=0.02):
        self._max_workers = int(max_workers)
        self._delay = float(delay)
        self.calls = 0
        self._active = 0
        self.max_active = 0
        self._lock = threading.Lock()
        self._executor = concurrent.futures.ThreadPoolExecutor(max_workers=self._max_workers)

    def submit_task(self, sample_info):
        self.calls += 1

        def _worker():
            with self._lock:
                self._active += 1
                self.max_active = max(self.max_active, self._active)
            try:
                time.sleep(self._delay)
                return float(sample_info.get("params", {}).get("x", 0.0))
            finally:
                with self._lock:
                    self._active -= 1

        return self._executor.submit(_worker)

    def shutdown(self):
        self._executor.shutdown(wait=True, cancel_futures=True)


class _FakeBucketAllocator:
    def __init__(self, out_dir):
        self.out_dir = out_dir

    def next_bucket_dir(self):
        return self.out_dir


class _FakeChain:
    def __init__(self):
        self.last_loglikelihood = 0.0

    def update(self, value):
        self.last_loglikelihood = float(value)


class _FakeVar:
    def __init__(self, name):
        self.name = name

    def map_standard_random_to_distribution(self, value):
        return float(value)


class _FakeBridsonVar(_FakeVar):
    def __init__(self, name, length=1.0):
        super().__init__(name)
        self._parameters = {"length": float(length)}


class _FakeResults(dict):
    def summary(self):
        return "fake-summary"


class _FakeDynamicNestedSampler:
    def __init__(self, loglikelihood, prior_transform, ndim, **_kwargs):
        self.loglikelihood = loglikelihood
        self.prior_transform = prior_transform
        self.ndim = int(ndim)
        self.results = _FakeResults()
        self.run_calls = 0
        self.last_resume = None

    def run_nested(self, resume=False, **_kwargs):
        self.run_calls += 1
        self.last_resume = bool(resume)
        uu = np.array(
            [
                [0.1, 0.2],
                [0.2, 0.3],
                [0.3, 0.4],
            ],
            dtype=float,
        )
        logl = []
        sample_uid = []
        for row in uu:
            point = self.prior_transform(row)
            sample_uid.append(point[-1])
            logl.append(float(self.loglikelihood(point)))
        n = len(logl)
        self.results = _FakeResults(
            {
                "samples_uid": np.array(sample_uid),
                "logwt": np.zeros(n),
                "logl": np.array(logl, dtype=float),
                "logvol": np.linspace(-1.0, -0.2, n),
                "logz": np.linspace(0.1, 0.3, n),
                "logzerr": np.full(n, 0.01),
                "samples_n": np.full(n, 10),
                "ncall": np.full(n, 1),
                "samples_it": np.arange(n),
                "samples_id": np.arange(n),
                "information": np.full(n, 0.2),
                "samples": uu.copy(),
                "samples_u": uu.copy(),
            }
        )


class _FakeDynamicNestedSamplerUsesPool:
    def __init__(self, loglikelihood, prior_transform, ndim, pool=None, **_kwargs):
        self.loglikelihood = loglikelihood
        self.prior_transform = prior_transform
        self.ndim = int(ndim)
        self.pool = pool
        self.results = _FakeResults()

    def run_nested(self, **_kwargs):
        uu = np.array(
            [
                [0.1, 0.2],
                [0.2, 0.3],
                [0.3, 0.4],
                [0.4, 0.5],
                [0.5, 0.6],
                [0.6, 0.7],
            ],
            dtype=float,
        )
        params = [self.prior_transform(row) for row in uu]
        if self.pool is None:
            logl = [float(self.loglikelihood(point)) for point in params]
        else:
            logl = [float(v) for v in self.pool.map(self.loglikelihood, params)]
        n = len(logl)
        self.results = _FakeResults(
            {
                "samples_uid": np.array([p[-1] for p in params]),
                "logwt": np.zeros(n),
                "logl": np.array(logl, dtype=float),
                "logvol": np.linspace(-1.0, -0.2, n),
                "logz": np.linspace(0.1, 0.3, n),
                "logzerr": np.full(n, 0.01),
                "samples_n": np.full(n, 10),
                "ncall": np.full(n, 1),
                "samples_it": np.arange(n),
                "samples_id": np.arange(n),
                "information": np.full(n, 0.2),
                "samples": uu.copy(),
                "samples_u": uu.copy(),
            }
        )


class _FakeNestedSampler:
    def __init__(self, loglikelihood, prior_transform, ndim, **_kwargs):
        self.loglikelihood = loglikelihood
        self.prior_transform = prior_transform
        self.ndim = int(ndim)
        self.results = _FakeResults()
        self.run_calls = 0
        self.last_resume = None

    def run_nested(self, resume=False, **_kwargs):
        self.run_calls += 1
        self.last_resume = bool(resume)
        uu = np.array(
            [
                [0.1, 0.2],
                [0.2, 0.3],
                [0.3, 0.4],
            ],
            dtype=float,
        )
        logl = []
        sample_uid = []
        for row in uu:
            point = self.prior_transform(row)
            sample_uid.append(point[-1])
            logl.append(float(self.loglikelihood(point)))
        n = len(logl)
        self.results = _FakeResults(
            {
                "samples_uid": np.array(sample_uid),
                "logwt": np.zeros(n),
                "logl": np.array(logl, dtype=float),
                "logvol": np.linspace(-1.0, -0.2, n),
                "logz": np.linspace(0.1, 0.3, n),
                "logzerr": np.full(n, 0.01),
                "samples_n": np.full(n, 10),
                "ncall": np.full(n, 1),
                "samples_it": np.arange(n),
                "samples_id": np.arange(n),
                "information": np.full(n, 0.2),
                "samples": uu.copy(),
                "samples_u": uu.copy(),
            }
        )


class _FakeNestedSamplerStrict:
    def __init__(self, loglikelihood, prior_transform, ndim, nlive=None, pool=None, rstate=None, queue_size=None):
        self.loglikelihood = loglikelihood
        self.prior_transform = prior_transform
        self.ndim = int(ndim)
        self.nlive = nlive
        self.pool = pool
        self.rstate = rstate
        self.queue_size = queue_size
        self.results = _FakeResults()

    def run_nested(self, **_kwargs):
        uu = np.array(
            [
                [0.1, 0.2],
                [0.2, 0.3],
                [0.3, 0.4],
            ],
            dtype=float,
        )
        logl = []
        sample_uid = []
        for row in uu:
            point = self.prior_transform(row)
            sample_uid.append(point[-1])
            logl.append(float(self.loglikelihood(point)))
        n = len(logl)
        self.results = _FakeResults(
            {
                "samples_uid": np.array(sample_uid),
                "logwt": np.zeros(n),
                "logl": np.array(logl, dtype=float),
                "logvol": np.linspace(-1.0, -0.2, n),
                "logz": np.linspace(0.1, 0.3, n),
                "logzerr": np.full(n, 0.01),
                "ncall": np.full(n, 1),
                "samples_n": np.full(n, 10),
                "samples_it": np.arange(n),
                "samples_id": np.arange(n),
                "information": np.full(n, 0.2),
                "samples": uu.copy(),
                "samples_u": uu.copy(),
            }
        )


class _FakeNestedSamplerWithInnerLogger(_FakeNestedSampler):
    last_inner_logger = None

    def __init__(self, loglikelihood, prior_transform, ndim, inner_logger=None, **_kwargs):
        super().__init__(loglikelihood, prior_transform, ndim, **_kwargs)
        type(self).last_inner_logger = inner_logger
        self._inner_logger = inner_logger

    def run_nested(self, **_kwargs):
        if self._inner_logger is not None:
            self._inner_logger.info("nested sampler internal log")
        return super().run_nested(**_kwargs)


class _InMemoryDb:
    def __init__(self):
        self.rows = []

    def add_data(self, row):
        self.rows.append(dict(row))


class _GoodPool:
    def execute(self, observables, _sample_info):
        return {"ok": observables.get("seed", 0) + 1}


class _FailPool:
    def execute(self, _observables, _sample_info):
        raise RuntimeError("forced module failure")


class SamplerHardeningSmokeTests(unittest.TestCase):
    def setUp(self):
        _FakeSample.reset()
        ModuleManager._instance = None
        self.tempdir = tempfile.TemporaryDirectory()
        os.environ["MPLCONFIGDIR"] = self.tempdir.name

    def tearDown(self):
        self.tempdir.cleanup()
        ModuleManager._instance = None

    def _sample_cfg(self):
        return {"sample": {"save_dir": self.tempdir.name, "task_result_dir": self.tempdir.name}}

    def _sample_cfg_with_bucket_root(self):
        sample_root = os.path.join(self.tempdir.name, "SAMPLE")
        os.makedirs(sample_root, exist_ok=True)
        return {
            "sample": {
                "save_dir": self.tempdir.name,
                "task_result_dir": self.tempdir.name,
                "sample_dirs": sample_root,
                "archive_samples": False,
            }
        }

    @staticmethod
    def _fake_torch_modules():
        fake_torch = types.ModuleType("torch")
        fake_nn = types.ModuleType("torch.nn")

        class _Layer:
            def __init__(self, *_args, **_kwargs):
                return

        class _NoGrad:
            def __enter__(self):
                return None

            def __exit__(self, *_exc):
                return False

        fake_nn.Linear = _Layer
        fake_nn.ReLU = _Layer
        fake_nn.Dropout = _Layer
        fake_nn.MSELoss = _Layer
        fake_nn.BCELoss = _Layer
        fake_nn.Sequential = lambda *_args, **_kwargs: object()

        fake_torch.nn = fake_nn
        fake_torch.backends = types.SimpleNamespace(mps=types.SimpleNamespace(is_available=lambda: False))
        fake_torch.cuda = types.SimpleNamespace(is_available=lambda: False)
        fake_torch.xpu = types.SimpleNamespace(is_available=lambda: False)
        fake_torch.xla = types.SimpleNamespace(is_available=lambda: False)
        fake_torch.hip = types.SimpleNamespace(is_available=lambda: False)
        fake_torch.optim = types.SimpleNamespace(Adam=lambda *_args, **_kwargs: object())
        fake_torch.device = lambda name: name
        fake_torch.float32 = np.float32
        fake_torch.tensor = lambda *_args, **_kwargs: np.array([])
        fake_torch.no_grad = _NoGrad
        fake_torch.save = lambda *_args, **_kwargs: None
        return fake_torch, fake_nn

    def _load_dnn_class(self):
        import importlib

        fake_torch, fake_nn = self._fake_torch_modules()
        fake_mpl = types.ModuleType("matplotlib")
        fake_pyplot = types.ModuleType("matplotlib.pyplot")
        fake_pyplot.figure = lambda *_args, **_kwargs: None
        fake_pyplot.savefig = lambda *_args, **_kwargs: None
        fake_pyplot.close = lambda *_args, **_kwargs: None
        fake_pyplot.colorbar = lambda *_args, **_kwargs: None
        fake_pyplot.draw = lambda *_args, **_kwargs: None
        fake_mpl.pyplot = fake_pyplot

        with patch.dict(
            sys.modules,
            {
                "torch": fake_torch,
                "torch.nn": fake_nn,
                "matplotlib": fake_mpl,
                "matplotlib.pyplot": fake_pyplot,
            },
            clear=False,
        ):
            if "jarvishep.dataconvert" in sys.modules:
                importlib.reload(sys.modules["jarvishep.dataconvert"])
            dnn_module = importlib.import_module("jarvishep.Sampling.dnn")
            dnn_module = importlib.reload(dnn_module)
            return dnn_module.DNN

    def test_random_sampler_smoke_closes_samples(self):
        sampler = RandomS()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.factory = _ImmediateFactory()

        points = iter([{"x": 0.1}, {"x": 0.2}, {"x": 0.3}])
        sampler.next_sample = lambda: next(points)

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.randoms.Sample", _FakeSample):
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_FakeSample.close_calls, 3)

    def test_grid_sampler_smoke_closes_samples(self):
        sampler = Grid()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.factory = _ImmediateFactory()
        sampler.tasks = []
        sampler.future_to_sample = {}

        points = iter([{"x": 0.1}, {"x": 0.2}, {"x": 0.3}])
        sampler.next_sample = lambda: next(points)

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.grid.Sample", _FakeSample):
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_FakeSample.close_calls, 3)

    def test_random_sampler_uses_bucketed_sample_save_dir_and_completion_tracking(self):
        sampler = RandomS()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg_with_bucket_root()
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Point number": 3,
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 200, "width": 6}},
            }
        )

        completed = []
        original_on_completed = sampler._on_sample_completed

        def _capture_completed(sample_info):
            completed.append(sample_info.get("uuid"))
            original_on_completed(sample_info)

        sampler._on_sample_completed = _capture_completed

        points = iter([{"x": 0.1}, {"x": 0.2}, {"x": 0.3}])
        sampler.next_sample = lambda: next(points)

        _CaptureConfigSample.reset()
        with patch("jarvishep.Sampling.randoms.Sample", _CaptureConfigSample):
            sampler.run_nested()

        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_CaptureConfigSample.close_calls, 3)
        self.assertEqual(len(completed), 3)
        self.assertGreaterEqual(len(_CaptureConfigSample.save_dirs), 1)
        sample_root = sampler.info["sample"]["sample_dirs"]
        for save_dir in _CaptureConfigSample.save_dirs:
            self.assertTrue(save_dir.startswith(sample_root))
            self.assertEqual(os.path.basename(save_dir), "000001")

    def test_grid_sampler_uses_bucketed_sample_save_dir_and_completion_tracking(self):
        sampler = Grid()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg_with_bucket_root()
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0, "num": 4}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 200, "width": 6}},
            }
        )

        completed = []
        original_on_completed = sampler._on_sample_completed

        def _capture_completed(sample_info):
            completed.append(sample_info.get("uuid"))
            original_on_completed(sample_info)

        sampler._on_sample_completed = _capture_completed

        points = iter([{"x": 0.1}, {"x": 0.2}, {"x": 0.3}])
        sampler.next_sample = lambda: next(points)

        _CaptureConfigSample.reset()
        with patch("jarvishep.Sampling.grid.Sample", _CaptureConfigSample):
            sampler.run_nested()

        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_CaptureConfigSample.close_calls, 3)
        self.assertEqual(len(completed), 3)
        self.assertGreaterEqual(len(_CaptureConfigSample.save_dirs), 1)
        sample_root = sampler.info["sample"]["sample_dirs"]
        for save_dir in _CaptureConfigSample.save_dirs:
            self.assertTrue(save_dir.startswith(sample_root))
            self.assertEqual(os.path.basename(save_dir), "000001")

    def test_mcmc_sampler_smoke_closes_samples(self):
        sampler = MCMC()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 2,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch(
            "jarvishep.Sampling.Source.MCMC.state_machine_base.Sample",
            _FakeSample,
        ), patch(
            "jarvishep.Sampling.Source.MCMC.state_machine_multistage_base.Sample",
            _FakeSample,
        ):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertGreaterEqual(sampler.factory.calls, 4)
        self.assertEqual(_FakeSample.close_calls, sampler.factory.calls)

    def test_tpmcmc_sampler_smoke_closes_samples(self):
        sampler = PTMCMC()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 2,
                        "num_iters": 2,
                        "exchange_interval": 1,
                        "proposal_scales": [0.1, 0.2],
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch(
            "jarvishep.Sampling.Source.MCMC.state_machine_base.Sample",
            _FakeSample,
        ), patch(
            "jarvishep.Sampling.Source.MCMC.state_machine_multistage_base.Sample",
            _FakeSample,
        ):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertGreaterEqual(sampler.factory.calls, 4)
        self.assertEqual(_FakeSample.close_calls, sampler.factory.calls)

    def test_ammcmc_sampler_smoke_closes_samples(self):
        sampler = AMMCMC()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 2,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "adapt_enabled": True,
                        "adapt_start_iter": 1,
                        "adapt_window": 1,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch(
            "jarvishep.Sampling.Source.MCMC.state_machine_base.Sample",
            _FakeSample,
        ), patch(
            "jarvishep.Sampling.Source.MCMC.state_machine_multistage_base.Sample",
            _FakeSample,
        ):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertGreaterEqual(sampler.factory.calls, 4)
        self.assertEqual(_FakeSample.close_calls, sampler.factory.calls)

    def test_demcmc_sampler_smoke_closes_samples(self):
        sampler = DEMCMC()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 3,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "de_gamma": 0.0,
                        "de_noise": 1.0e-3,
                        "de_crossover": 0.9,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 6)
        self.assertEqual(_FakeSample.close_calls, 6)

    def test_dream_lite_sampler_smoke_closes_samples(self):
        sampler = DREAMLite()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 3,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "de_gamma": 0.0,
                        "de_noise": 1.0e-3,
                        "de_crossover": 0.9,
                        "dream_snooker_prob": 0.2,
                        "dream_archive_size": 64,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 6)
        self.assertEqual(_FakeSample.close_calls, 6)

    def test_dream_sampler_smoke_closes_samples(self):
        sampler = DREAM()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 3,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "de_gamma": 0.0,
                        "de_noise": 1.0e-3,
                        "de_crossover": 0.9,
                        "dream_snooker_prob": 0.2,
                        "dream_archive_size": 128,
                        "dream_crossover_values": [0.2, 0.5, 0.9],
                        "dream_crossover_adapt_interval": 16,
                        "dream_scale_jitter": 0.1,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 6)
        self.assertEqual(_FakeSample.close_calls, 6)

    def test_ensemblemcmc_sampler_smoke_closes_samples(self):
        sampler = EnsembleMCMC()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 4,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "stretch_a": 2.0,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 8)
        self.assertEqual(_FakeSample.close_calls, 8)

    def test_ptensemble_sampler_smoke_closes_samples(self):
        sampler = PTEnsemble()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 4,
                        "num_iters": 2,
                        "exchange_interval": 1,
                        "proposal_scales": [0.1, 0.1, 0.1, 0.1],
                        "temperature_ladder": [1.0, 1.6, 2.4, 3.6],
                        "stretch_a": 2.0,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 8)
        self.assertEqual(_FakeSample.close_calls, 8)

    def test_slicemcmc_sampler_smoke_closes_samples(self):
        sampler = SliceMCMC()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 3,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "slice_mode": "random_direction",
                        "slice_width": 0.2,
                        "slice_max_steps_out": 16,
                        "slice_max_shrink": 32,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 6)
        self.assertEqual(_FakeSample.close_calls, 6)

    def test_ess_sampler_smoke_closes_samples(self):
        sampler = ESS()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 3,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "ess_prior_cov": 1.0,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 6)
        self.assertEqual(_FakeSample.close_calls, 6)

    def test_mala_sampler_smoke_closes_samples(self):
        sampler = MALA()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 3,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "mala_step_size": 0.1,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 6)
        self.assertEqual(_FakeSample.close_calls, 6)

    def test_hmc_sampler_smoke_closes_samples(self):
        sampler = HMC()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 3,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "hmc_step_size": 0.05,
                        "hmc_leapfrog_steps": 6,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 6)
        self.assertEqual(_FakeSample.close_calls, 6)

    def test_nuts_sampler_smoke_closes_samples(self):
        sampler = NUTS()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 3,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "nuts_step_size": 0.05,
                        "nuts_max_depth": 5,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 6)
        self.assertEqual(_FakeSample.close_calls, 6)

    def test_robustam_sampler_smoke_closes_samples(self):
        sampler = RobustAM()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 2,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "adapt_enabled": True,
                        "adapt_start_iter": 1,
                        "adapt_window": 1,
                        "global_jump_prob": 0.2,
                        "heavy_tail_df": 4.0,
                        "global_scale": 2.0,
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 4)
        self.assertEqual(_FakeSample.close_calls, 4)

    def test_dram_sampler_smoke_closes_samples(self):
        sampler = DRAM()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.info["sample"]["sample_dirs"] = self.tempdir.name
        sampler.info["sample"]["archive_samples"] = False
        sampler.factory = _ImmediateFactory()
        sampler.set_config(
            {
                "Sampling": {
                    "Bounds": {
                        "num_chains": 2,
                        "num_iters": 2,
                        "proposal_scale": 0.1,
                        "adapt_enabled": True,
                        "adapt_start_iter": 1,
                        "adapt_window": 1,
                        "dr_scale_factors": [1.0, 0.5],
                    },
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
                "Scan": {"sample_directory": {"limit": 20, "width": 4}},
            }
        )

        t0 = time.perf_counter()
        with patch(
            "jarvishep.Sampling.Source.MCMC.state_machine_base.Sample",
            _FakeSample,
        ), patch(
            "jarvishep.Sampling.Source.MCMC.state_machine_multistage_base.Sample",
            _FakeSample,
        ):
            sampler.initialize()
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertGreaterEqual(sampler.factory.calls, 4)
        self.assertEqual(_FakeSample.close_calls, sampler.factory.calls)

    def test_dnn_sampler_dataset_smoke_closes_samples(self):
        DNN = self._load_dnn_class()
        sampler = DNN()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.factory = _ImmediateFactory()
        sampler._outputs = ["obs1"]
        sampler._vcolum = ["x"]
        sampler._ocolum = ["obs1"]
        sampler._device = "cpu"
        sampler._iter = 0

        points = iter([{"x": 0.1}, {"x": 0.2}, {"x": 0.3}])
        sampler.next_sample = lambda: next(points)

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.dnn.Sample", _FakeSample):
            dataset = sampler.create_dataset(3)
        sampler.plotter.shutdown(wait=False, cancel_futures=True)
        self.assertLess(time.perf_counter() - t0, 1.0)
        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_FakeSample.close_calls, 3)
        self.assertEqual(dataset.data.shape[0], 3)

    def test_dnn_output_likelihood_symbol_match(self):
        DNN = self._load_dnn_class()
        sampler = DNN()
        sampler._outputs = ["z"]
        sampler._loglike = types.SimpleNamespace(variables={sp.Symbol("z")})
        self.assertTrue(sampler.check_output_in_likelihood())

    def test_dnn_run_info_json_handles_numpy_scalar(self):
        DNN = self._load_dnn_class()
        sampler = DNN()
        run_path = os.path.join(self.tempdir.name, "run_info.json")
        sampler.info = {
            "DNN_info": {
                "run": run_path,
                "iter_data": [],
                "loss": {"model_obs1": [np.float32(0.123)]},
                "model": {},
            }
        }

        sampler._write_run_info()
        sampler.plotter.shutdown(wait=False, cancel_futures=True)

        with open(run_path, "r") as f1:
            payload = json.load(f1)
        self.assertAlmostEqual(payload["loss"]["model_obs1"][0], 0.123, places=6)

    def test_dynesty_sampler_smoke_closes_samples(self):
        sampler = Dynesty()
        sampler.logger = _NoopLogger()
        sampler.info = {
            "sample": {"save_dir": self.tempdir.name, "task_result_dir": self.tempdir.name},
            "logfile": os.path.join(self.tempdir.name, "dynesty.log"),
        }
        sampler.factory = _ImmediateFactory()
        sampler._dimensions = 2
        sampler._nlive = 5
        sampler._rstate = np.random.default_rng(0)
        sampler._runnested = {}
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))
        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.DynamicNestedSampler = _FakeDynamicNestedSampler

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.dynesty.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        
        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_FakeSample.close_calls, 3)

    def test_multinest_sampler_smoke_closes_samples(self):
        sampler = MultiNest()
        sampler.logger = _NoopLogger()
        sampler.info = {
            "sample": {"save_dir": self.tempdir.name, "task_result_dir": self.tempdir.name},
            "logfile": os.path.join(self.tempdir.name, "multinest.log"),
        }
        sampler.factory = _ImmediateFactory()
        sampler._dimensions = 2
        sampler._nlive = 5
        sampler._rstate = np.random.default_rng(0)
        sampler._runnested = {}
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))
        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.NestedSampler = _FakeNestedSampler

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.multinest.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()
        self.assertLess(time.perf_counter() - t0, 1.0)
        
        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_FakeSample.close_calls, 3)
        sampler.finalize()
        summary_json = os.path.join(self.tempdir.name, "DATABASE", "multinest_summary.json")
        self.assertTrue(os.path.exists(summary_json))
        with open(summary_json, "r", encoding="utf-8") as f1:
            payload = json.load(f1)
        self.assertEqual(payload.get("method"), "MultiNest")
        self.assertIn("niter", payload)
        self.assertIn("ncall", payload)
        self.assertIn("logz", payload)

    def test_dynesty_sampler_native_checkpoint_resume(self):
        sampler = Dynesty()
        sampler.logger = _NoopLogger()
        sample_root = os.path.join(self.tempdir.name, "SAMPLE")
        sampler.info = {
            "sample": {
                "save_dir": self.tempdir.name,
                "task_result_dir": self.tempdir.name,
                "sample_dirs": sample_root,
                "archive_samples": False,
            },
            "logfile": os.path.join(self.tempdir.name, "dynesty.log"),
        }
        sampler.factory = _ImmediateFactory()
        config = {
            "Sampling": {
                "Method": "Dynesty",
                "Bounds": {"nlive": 5, "rseed": 1, "run_nested": {"maxiter": 1}},
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    },
                    {
                        "name": "y",
                        "description": "y",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    },
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }
        sampler.set_config(config)
        sampler._runnested = {"maxiter": 1}
        checkpoint_root = os.path.join(self.tempdir.name, "checkpoints", "dynesty")
        sampler.configure_runtime_checkpointing(checkpoint_root, logger=_NoopLogger())
        sampler.set_runtime_checkpoint_context(run_spec={"normalized_config": sampler.config}, factory_blueprint={})

        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.DynamicNestedSampler = _FakeDynamicNestedSampler
        with patch("jarvishep.Sampling.dynesty.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()
        self.assertEqual(sampler.sampler.run_calls, 1)
        self.assertFalse(sampler.sampler.last_resume)
        self.assertTrue(sampler.persist_runtime_checkpoint(force=True, reason="unit-test"))

        restored = Dynesty()
        restored.logger = _NoopLogger()
        restored.info = {
            "sample": {
                "save_dir": self.tempdir.name,
                "task_result_dir": self.tempdir.name,
                "sample_dirs": sample_root,
                "archive_samples": False,
            },
            "logfile": os.path.join(self.tempdir.name, "dynesty.log"),
        }
        restored.factory = _ImmediateFactory()
        restored.set_config(config)
        restored.configure_runtime_checkpointing(checkpoint_root, logger=_NoopLogger())
        self.assertTrue(restored.restore_runtime_checkpoint_if_available())
        with patch("jarvishep.Sampling.dynesty.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            restored.run_nested()
        self.assertEqual(restored.sampler.run_calls, 2)
        self.assertTrue(restored.sampler.last_resume)

    def test_multinest_sampler_native_checkpoint_resume(self):
        sampler = MultiNest()
        sampler.logger = _NoopLogger()
        sample_root = os.path.join(self.tempdir.name, "SAMPLE")
        sampler.info = {
            "sample": {
                "save_dir": self.tempdir.name,
                "task_result_dir": self.tempdir.name,
                "sample_dirs": sample_root,
                "archive_samples": False,
            },
            "logfile": os.path.join(self.tempdir.name, "multinest.log"),
        }
        sampler.factory = _ImmediateFactory()
        config = {
            "Sampling": {
                "Method": "MultiNest",
                "Bounds": {"nlive": 5, "rseed": 1, "run_nested": {"maxiter": 1}},
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    },
                    {
                        "name": "y",
                        "description": "y",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    },
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }
        sampler.set_config(config)
        sampler._runnested = {"maxiter": 1}
        checkpoint_root = os.path.join(self.tempdir.name, "checkpoints", "multinest")
        sampler.configure_runtime_checkpointing(checkpoint_root, logger=_NoopLogger())
        sampler.set_runtime_checkpoint_context(run_spec={"normalized_config": sampler.config}, factory_blueprint={})

        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.NestedSampler = _FakeNestedSampler
        with patch("jarvishep.Sampling.multinest.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()
        self.assertEqual(sampler.sampler.run_calls, 1)
        self.assertFalse(sampler.sampler.last_resume)
        self.assertTrue(sampler.persist_runtime_checkpoint(force=True, reason="unit-test"))

        restored = MultiNest()
        restored.logger = _NoopLogger()
        restored.info = {
            "sample": {
                "save_dir": self.tempdir.name,
                "task_result_dir": self.tempdir.name,
                "sample_dirs": sample_root,
                "archive_samples": False,
            },
            "logfile": os.path.join(self.tempdir.name, "multinest.log"),
        }
        restored.factory = _ImmediateFactory()
        restored.set_config(config)
        restored.configure_runtime_checkpointing(checkpoint_root, logger=_NoopLogger())
        self.assertTrue(restored.restore_runtime_checkpoint_if_available())
        with patch("jarvishep.Sampling.multinest.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            restored.run_nested()
        self.assertEqual(restored.sampler.run_calls, 2)
        self.assertTrue(restored.sampler.last_resume)

    def test_distributor_supports_multinest_method(self):
        sampler = Distributor.set_method("MultiNest")
        self.assertEqual(getattr(sampler, "method", None), "MultiNest")
        self.assertTrue(str(getattr(sampler, "schema", "")).endswith("MultiNest_schema.json"))

    def test_multinest_sampler_accepts_nestedsampler_without_log_file_path(self):
        sampler = MultiNest()
        sampler.logger = _NoopLogger()
        sampler.info = {
            "sample": {"save_dir": self.tempdir.name, "task_result_dir": self.tempdir.name},
            "logfile": os.path.join(self.tempdir.name, "multinest.log"),
        }
        sampler.factory = _ImmediateFactory()
        sampler._dimensions = 2
        sampler._nlive = 5
        sampler._rstate = np.random.default_rng(0)
        sampler._runnested = {}
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))
        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.NestedSampler = _FakeNestedSamplerStrict

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.multinest.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()
            sampler.finalize()
        self.assertLess(time.perf_counter() - t0, 1.0)       
        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_FakeSample.close_calls, 3)

    def test_multinest_forwards_inner_logger_to_nestedsampler(self):
        sampler = MultiNest()
        cap_logger = _CaptureLogger()
        sampler.logger = cap_logger
        sampler.info = {
            "sample": {"save_dir": self.tempdir.name, "task_result_dir": self.tempdir.name},
            "logfile": os.path.join(self.tempdir.name, "multinest.log"),
        }
        sampler.factory = _ImmediateFactory()
        sampler._dimensions = 2
        sampler._nlive = 5
        sampler._rstate = np.random.default_rng(0)
        sampler._runnested = {}
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))
        _FakeNestedSamplerWithInnerLogger.last_inner_logger = None
        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.NestedSampler = _FakeNestedSamplerWithInnerLogger

        with patch("jarvishep.Sampling.multinest.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()

        self.assertIs(_FakeNestedSamplerWithInnerLogger.last_inner_logger, cap_logger)
        self.assertTrue(any("nested sampler internal log" in str(args[0]) for _, args, _ in cap_logger.records if args))

    def test_multinest_finalize_is_noop_when_sampler_missing(self):
        sampler = MultiNest()
        sampler.logger = _NoopLogger()
        sampler.sampler = None
        sampler.finalize()

    def test_multinest_combine_data_fallback_uses_converted_csv(self):
        sampler = MultiNest()
        sampler.logger = _NoopLogger()
        db_dir = os.path.join(self.tempdir.name, "DATABASE")
        os.makedirs(db_dir, exist_ok=True)
        sampler.info = {"sample": {"task_result_dir": self.tempdir.name}}
        sampler.df = pd.DataFrame(
            {
                "uuid": ["u-001", "u-002"],
                "log_Like": [1.0, 2.0],
                "log_PriorVolume": [-1.0, -2.0],
            }
        )
        sampler.lnX_from_LogLike = None

        converted_csv = os.path.join(db_dir, "samples.0.csv")
        with open(converted_csv, "w", newline="") as f1:
            writer = csv.DictWriter(f1, fieldnames=["uuid", "LogL", "obs"])
            writer.writeheader()
            writer.writerow({"uuid": "u-001", "LogL": 1.0, "obs": 10.0})
            writer.writerow({"uuid": "u-002", "LogL": 2.0, "obs": 20.0})

        with open(os.path.join(db_dir, "running.json"), "w", encoding="utf-8") as f1:
            json.dump({"converted": [converted_csv], "pathroot": os.path.join(db_dir, "samples"), "activeNO": 0}, f1)

        sampler.combine_data(os.path.join(db_dir, "samples.hdf5"))
        merged = os.path.join(db_dir, "multinest_full.csv")
        self.assertTrue(os.path.exists(merged))
        out = pd.read_csv(merged)
        self.assertEqual(out.shape[0], 2)

    def test_multinest_interpolator_handles_duplicate_loglike_points(self):
        sampler = MultiNest()
        sampler.logger = _NoopLogger()
        sampler.df = pd.DataFrame(
            {
                "log_Like": [1.0, 1.0, 2.0, 3.0],
                "log_PriorVolume": [-0.5, -0.7, -1.2, -1.8],
            }
        )
        sampler._build_logl_to_logx_interpolator()
        self.assertIsNotNone(sampler.lnX_from_LogLike)
        val = float(sampler.lnX_from_LogLike(1.5))
        self.assertTrue(np.isfinite(val))

    def test_dynesty_samples_u_mapping_uses_unit_cube_source(self):
        sampler = Dynesty()
        sampler.logger = _NoopLogger()
        sampler.info = {
            "sample": {"save_dir": self.tempdir.name, "task_result_dir": self.tempdir.name},
            "logfile": os.path.join(self.tempdir.name, "dynesty.log"),
        }
        sampler.factory = _ImmediateFactory()
        sampler._dimensions = 2
        sampler._nlive = 5
        sampler._rstate = np.random.default_rng(0)
        sampler._runnested = {}
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))

        class _FakeDynestyMapping:
            def __init__(self, loglikelihood, prior_transform, ndim, **_kwargs):
                self.loglikelihood = loglikelihood
                self.prior_transform = prior_transform
                self.ndim = int(ndim)
                self.results = _FakeResults()

            def run_nested(self, **_kwargs):
                uu = np.array(
                    [
                        [0.1, 0.2],
                        [0.2, 0.3],
                        [0.3, 0.4],
                        [0.4, 0.5],
                    ],
                    dtype=float,
                )
                vv = np.array(
                    [
                        [10.0, 20.0],
                        [20.0, 30.0],
                        [30.0, 40.0],
                        [40.0, 50.0],
                    ],
                    dtype=float,
                )
                sample_uid = []
                logl = []
                for row in uu:
                    point = self.prior_transform(row)
                    sample_uid.append(point[-1])
                    logl.append(float(self.loglikelihood(point)))
                n = len(sample_uid)
                self.results = _FakeResults(
                    {
                        "samples_uid": np.array(sample_uid),
                        "logwt": np.zeros(n),
                        "logl": np.array(logl, dtype=float),
                        "logvol": np.linspace(-1.0, -0.2, n),
                        "logz": np.linspace(0.1, 0.4, n),
                        "logzerr": np.full(n, 0.01),
                        "samples_n": np.full(n, 10),
                        "ncall": np.full(n, 1),
                        "samples_it": np.arange(n),
                        "samples_id": np.arange(n),
                        "information": np.full(n, 0.2),
                        "samples": vv,
                        "samples_u": uu,
                    }
                )

        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.DynamicNestedSampler = _FakeDynestyMapping

        with patch("jarvishep.Sampling.dynesty.Sample", _FakeSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()
            os.makedirs(os.path.dirname(sampler.info["db"]["nested result"]), exist_ok=True)
            sampler.finalize()

        out_csv = sampler.info["db"]["nested result"]
        with open(out_csv, "r", newline="") as f1:
            rows = list(csv.DictReader(f1))
        self.assertGreaterEqual(len(rows), 1)
        self.assertEqual(float(rows[0]["samples_u[0]"]), 0.1)
        self.assertEqual(float(rows[0]["samples_u[1]"]), 0.2)
        self.assertEqual(float(rows[0]["samples_v[0]"]), 10.0)
        self.assertEqual(float(rows[0]["samples_v[1]"]), 20.0)

    def test_dynesty_uses_new_pool_adapter(self):
        sampler = Dynesty()
        sampler.logger = _NoopLogger()
        sampler.info = {
            "sample": {"save_dir": self.tempdir.name, "task_result_dir": self.tempdir.name},
            "logfile": os.path.join(self.tempdir.name, "dynesty.log"),
        }
        sampler.factory = _ImmediateFactory()
        sampler._dimensions = 2
        sampler._nlive = 5
        sampler._rstate = np.random.default_rng(0)
        sampler._runnested = {}
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))
        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.DynamicNestedSampler = _FakeDynamicNestedSampler
        import jarvishep.Sampling.dynesty as dynesty_module
        pool_build_counter = {"count": 0}
        original_pool_cls = dynesty_module.DynestyFactoryPool

        class _SpyPool(original_pool_cls):
            def __init__(self, *args, **kwargs):
                pool_build_counter["count"] += 1
                super().__init__(*args, **kwargs)

        with patch("jarvishep.Sampling.dynesty.DynestyFactoryPool", _SpyPool), patch(
            "jarvishep.Sampling.dynesty.Sample",
            _FakeSample,
        ), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()
        self.assertGreaterEqual(pool_build_counter["count"], 1)

        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_FakeSample.close_calls, 3)

    def test_dynesty_backpressure_limits_factory_concurrency(self):
        sampler = Dynesty()
        sampler.logger = _NoopLogger()
        sampler.info = {
            "sample": {"save_dir": self.tempdir.name, "task_result_dir": self.tempdir.name},
            "logfile": os.path.join(self.tempdir.name, "dynesty.log"),
        }
        blocking_factory = _BlockingFactory(max_workers=16, delay=0.02)
        sampler.factory = blocking_factory
        sampler.set_max_workers(6)
        sampler.set_execution_profile(dynesty_workers=6, max_pending_factory=2)
        sampler._dimensions = 2
        sampler._nlive = 12
        sampler._rstate = np.random.default_rng(0)
        sampler._runnested = {}
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))

        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.DynamicNestedSampler = _FakeDynamicNestedSamplerUsesPool

        try:
            with patch("jarvishep.Sampling.dynesty.Sample", _FakeSample), patch.dict(
                sys.modules,
                {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
                clear=False,
            ):
                sampler.run_nested()
            self.assertGreaterEqual(blocking_factory.calls, 6)
            self.assertLessEqual(blocking_factory.max_active, 2)
            self.assertEqual(_FakeSample.close_calls, blocking_factory.calls)
        finally:
            blocking_factory.shutdown()

    def test_dynesty_uses_bucketed_sample_save_dir(self):
        sampler = Dynesty()
        sampler.logger = _NoopLogger()
        sample_root = os.path.join(self.tempdir.name, "SAMPLE")
        sampler.info = {
            "sample": {
                "save_dir": self.tempdir.name,
                "task_result_dir": self.tempdir.name,
                "sample_dirs": sample_root,
            },
            "logfile": os.path.join(self.tempdir.name, "dynesty.log"),
        }
        sampler.factory = _ImmediateFactory()
        sampler._dimensions = 2
        sampler._nlive = 5
        sampler._rstate = np.random.default_rng(0)
        sampler._runnested = {}
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))

        fake_dynesty_mod = types.ModuleType("jarvishep.Sampling.Source.Dynesty.py.dynesty")
        fake_dynesty_mod.DynamicNestedSampler = _FakeDynamicNestedSampler

        _CaptureConfigSample.reset()
        with patch("jarvishep.Sampling.dynesty.Sample", _CaptureConfigSample), patch.dict(
            sys.modules,
            {"jarvishep.Sampling.Source.Dynesty.py.dynesty": fake_dynesty_mod},
            clear=False,
        ):
            sampler.run_nested()

        self.assertGreaterEqual(len(_CaptureConfigSample.save_dirs), 1)
        for save_dir in _CaptureConfigSample.save_dirs:
            self.assertTrue(save_dir.startswith(sample_root))
            self.assertEqual(os.path.basename(save_dir), "000001")

    def test_bridson_sampler_smoke_closes_samples(self):
        sampler = Bridson()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.factory = _ImmediateFactory()
        sampler.bucket_alloc = _FakeBucketAllocator(self.tempdir.name)
        sampler.total_core = 2
        sampler.tasks = {}

        points = iter([{"x": 0.1}, {"x": 0.2}, {"x": 0.3}])
        sampler.next_sample = lambda: next(points)

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.bridson.Sample", _FakeSample):
            sampler.run_wo_nuisance()
        self.assertLess(time.perf_counter() - t0, 1.0)
        
        self.assertEqual(sampler.factory.calls, 3)
        self.assertEqual(_FakeSample.close_calls, 3)

    def test_bridson_progress_logging_uses_info_for_permille_and_warning_for_percent(self):
        sampler = Bridson()
        sampler.logger = _CaptureLogger()
        sampler._P = np.zeros((1000, 2), dtype=float)
        sampler.barinfo = {}

        sampler.progress_bar()
        sampler._index = 1
        sampler.progress_bar()
        sampler._index = 10
        sampler.progress_bar()

        progress_records = [
            (level, args[0])
            for level, args, _kwargs in sampler.logger.records
            if args and "samples submited" in str(args[0])
        ]
        self.assertEqual(progress_records[0][0], "INFO")
        self.assertIn("0‰", progress_records[0][1])
        self.assertEqual(progress_records[1][0], "INFO")
        self.assertIn("1‰", progress_records[1][1])
        self.assertEqual(progress_records[2][0], "WARNING")
        self.assertIn("10‰", progress_records[2][1])

    def test_diver_sampler_smoke_closes_samples(self):
        sampler = Diver()
        sampler.logger = _NoopLogger()
        sampler.info = self._sample_cfg()
        sampler.factory = _ImmediateFactory()
        sampler.bucket_alloc = _FakeBucketAllocator(self.tempdir.name)
        sampler.vars = (_FakeVar("x"), _FakeVar("y"))

        population = np.array([[0.1, 0.2], [0.2, 0.3], [0.3, 0.4]], dtype=float)

        t0 = time.perf_counter()
        with patch("jarvishep.Sampling.diver.Sample", _FakeSample):
            values = sampler._evaluate_population_loglike(population)
        self.assertLess(time.perf_counter() - t0, 1.0)
        
        self.assertEqual(sampler.factory.calls, 3)
        
        self.assertEqual(_FakeSample.close_calls, 3)
        self.assertEqual(values.shape[0], 3)

    def test_selection_contract_restored(self):
        sampler = RandomS()
        self.assertTrue(sampler.evaluate_selection("x > 0.5", {"x": 0.7}))
        self.assertFalse(sampler.evaluate_selection("x > 0.5", {"x": 0.2}))

    def test_bridson_selection_filters_generated_points(self):
        sampler = Bridson()
        sampler.logger = _NoopLogger()
        sampler.vars = (
            _FakeBridsonVar("HiggsMass"),
            _FakeBridsonVar("BottomMass"),
        )
        sampler._P = np.array(
            [
                [0.20, 0.20],
                [0.90, 0.30],
                [0.30, 0.25],
            ],
            dtype=float,
        )
        sampler._selectionexp = "2.0 * BottomMass < HiggsMass"

        sample = sampler.next_sample()

        self.assertEqual(sample, {"HiggsMass": 0.9, "BottomMass": 0.3})
        self.assertEqual(sampler._index, 2)

    def test_module_failure_policy_fail_fast_override(self):
        manager = ModuleManager()
        manager.set_logger(_NoopLogger())
        manager.set_config({"Sampling": {"ModuleFailurePolicy": "fail-fast"}})
        manager.workflow = {2: ["bad"]}
        manager.module_pools = {"bad": _FailPool()}
        manager._database = _InMemoryDb()

        sample_info = {"uuid": "s-001", "observables": {"seed": 1}, "logger": _NoopLogger()}
        with self.assertRaises(RuntimeError):
            manager.execute_workflow(sample_info)

    def test_module_failure_policy_continue_default_skips_remaining_modules_for_sample(self):
        manager = ModuleManager()
        manager.set_logger(_NoopLogger())
        manager.set_config({"Sampling": {}})
        manager.workflow = {2: ["bad", "good"]}
        manager.module_pools = {"bad": _FailPool(), "good": _GoodPool()}
        manager._database = _InMemoryDb()

        sample_info = {"uuid": "s-002", "observables": {"seed": 1}, "logger": _NoopLogger()}
        result = manager.execute_workflow(sample_info)
        self.assertEqual(result, 1.0)
        self.assertNotIn("ok", sample_info["observables"])
        self.assertEqual(len(manager.database.rows), 1)
        self.assertEqual(manager.database.rows[0]["seed"], 1)


if __name__ == "__main__":
    unittest.main()
