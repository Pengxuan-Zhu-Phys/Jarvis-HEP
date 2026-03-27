#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from copy import deepcopy
import concurrent.futures
from unittest.mock import patch

import numpy as np
import pandas as pd


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Sampling.dnn import DNN  # noqa: E402
from jarvishep.Sampling.diver import Diver  # noqa: E402
from jarvishep.Sampling.dynesty import Dynesty  # noqa: E402
from jarvishep.Sampling.multinest import MultiNest  # noqa: E402
from jarvishep.Sampling.Source.Diver.de import DEConfig, DifferentialEvolution  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _ImmediateFactory:
    def submit_task(self, sample_info):
        future = concurrent.futures.Future()
        params = sample_info.get("params", {})
        logl = float(sum(float(value) for value in params.values())) if isinstance(params, dict) else 0.0
        sample_info.setdefault("observables", {})
        sample_info["observables"]["LogL"] = logl
        future.set_result(logl)
        return future


class _FakeDiverSample:
    seq = 0

    @classmethod
    def reset(cls):
        cls.seq = 0

    def __init__(self, params):
        type(self).seq += 1
        self.uuid = f"diver-fake-{type(self).seq:04d}"
        self.params = dict(params)
        self.info = {
            "uuid": self.uuid,
            "params": dict(self.params),
            "observables": {**self.params, "uuid": self.uuid},
        }

    def set_config(self, config):
        if isinstance(config, dict):
            self.info.update(dict(config))
        self.info["uuid"] = self.uuid
        self.info["params"] = dict(self.params)
        observables = dict(self.info.get("observables", {}))
        observables.update(self.params)
        observables["uuid"] = self.uuid
        self.info["observables"] = observables

    def close(self):
        return None


def _sanitize_config(config: dict) -> dict:
    config = deepcopy(config)
    config.setdefault("Scan", {})
    config.setdefault("Sampling", {})
    config.setdefault("EnvReqs", {"OS": [], "Python": {"version": ">=3.10", "Dependencies": []}})
    config.setdefault("Calculators", {"make_paraller": 1, "Modules": []})
    return config


class TestCheckpointPayloadRoundtrip(unittest.TestCase):
    def test_dnn_sampler_checkpoint_payload_roundtrip(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-dnn-") as tmp:
            sampler = DNN()
            sampler.logger = _NoopLogger()
            sampler.path["task_root"] = tmp
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_config(
                _sanitize_config(
                    {
                        "Sampling": {
                            "Method": "DNN",
                            "Variables": [
                                {
                                    "name": "x",
                                    "description": "x",
                                    "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                                }
                            ],
                            "Bounds": {
                                "Niters": 3,
                                "Hidden_layers": [4],
                                "Batch_size": 2,
                                "Ninit": 2,
                                "Nepoch": 1,
                                "Learning_rate": 0.001,
                                "Outputs": ["LogL"],
                                "Prop_new": 0.5,
                            },
                        }
                    }
                )
            )
            sampler.initialize()
            sampler.dataset = pd.DataFrame([{"x": 0.1, "LogL": 1.0}])
            sampler.df = pd.DataFrame([{"x": 0.2}])
            sampler._iter = 2
            sampler.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            sampler.set_runtime_checkpoint_context(run_spec={"normalized_config": sampler.config}, factory_blueprint={})
            payload = sampler.export_runtime_state()

            restored = DNN()
            restored.logger = _NoopLogger()
            restored.path["task_root"] = tmp
            restored.info["sample"] = sampler.info["sample"]
            restored.set_config(deepcopy(sampler.config))
            restored.initialize()
            restored.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            restored.import_runtime_state(payload)
            self.assertEqual(restored._iter, 2)
            self.assertEqual(int(restored.dataset.shape[0]), 1)
            self.assertEqual(int(restored.df.shape[0]), 1)

    def test_diver_engine_state_roundtrip(self):
        cfg = DEConfig(dim=2, pop_size=5, max_gen=3, max_civ=2, seed=7)
        engine = DifferentialEvolution(cfg)
        payload = engine.export_state()
        clone = DifferentialEvolution(cfg)
        clone.import_state(payload)
        self.assertEqual(clone.cfg.dim, 2)
        self.assertEqual(clone.cfg.pop_size, 5)
        self.assertEqual(clone.export_state()["mean_cost"], payload["mean_cost"])

    def test_diver_sampler_checkpoint_payload_roundtrip(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-diver-") as tmp:
            sampler = Diver()
            sampler.logger = _NoopLogger()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_config(
                _sanitize_config(
                    {
                        "Sampling": {
                            "Method": "Diver",
                            "Variables": [
                                {
                                    "name": "x",
                                    "description": "x",
                                    "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                                }
                            ],
                            "run": {"NP": 5, "maxgen": 3, "maxciv": 2, "seed": 7},
                        }
                    }
                )
            )
            sampler._best_params = {"x": 0.5}
            sampler._de_result = {"best_cost": 1.0, "best_vector": np.array([0.1, 0.2])}
            sampler.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            sampler.set_runtime_checkpoint_context(run_spec={"normalized_config": sampler.config}, factory_blueprint={})
            payload = sampler.export_runtime_state()

            restored = Diver()
            restored.logger = _NoopLogger()
            restored.info["sample"] = sampler.info["sample"]
            restored.set_config(deepcopy(sampler.config))
            restored.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            restored.import_runtime_state(payload)
            self.assertEqual(restored._best_params, {"x": 0.5})
            self.assertEqual(restored._de_result["best_cost"], 1.0)

    def test_diver_runtime_checkpoint_resume_continues_from_next_civilization(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-diver-resume-") as tmp:
            checkpoint_root = os.path.join(tmp, "checkpoints", "diver")

            sampler = Diver()
            sampler.logger = _NoopLogger()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_config(
                _sanitize_config(
                    {
                        "Sampling": {
                            "Method": "Diver",
                            "Variables": [
                                {
                                    "name": "x",
                                    "description": "x",
                                    "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                                }
                            ],
                            "run": {"NP": 4, "maxgen": 1, "maxciv": 2, "seed": 11},
                        }
                    }
                )
            )
            sampler.set_factory(_ImmediateFactory())
            sampler.configure_runtime_checkpointing(checkpoint_root, logger=_NoopLogger())

            _FakeDiverSample.reset()
            call_count = {"n": 0}
            original_persist = sampler.persist_runtime_checkpoint

            def _persist_once_then_cutoff(*args, **kwargs):
                call_count["n"] += 1
                saved = original_persist(*args, **kwargs)
                if call_count["n"] >= 2:
                    raise RuntimeError("intentional diver checkpoint cutoff")
                return saved

            sampler.persist_runtime_checkpoint = _persist_once_then_cutoff  # type: ignore[assignment]

            with patch("jarvishep.Sampling.diver.Sample", _FakeDiverSample):
                with self.assertRaisesRegex(RuntimeError, "intentional diver checkpoint cutoff"):
                    sampler.run_nested()

            checkpoint_file = os.path.join(checkpoint_root, "state.pkl")
            self.assertTrue(os.path.exists(checkpoint_file))
            self.assertIsNotNone(sampler._de_state)
            self.assertEqual(int(sampler._de_state["runtime_state"]["next_civ"]), 2)

            restored = Diver()
            restored.logger = _NoopLogger()
            restored.info["sample"] = sampler.info["sample"]
            restored.set_config(deepcopy(sampler.config))
            restored.set_factory(_ImmediateFactory())
            restored.configure_runtime_checkpointing(checkpoint_root, logger=_NoopLogger())
            self.assertTrue(restored.restore_runtime_checkpoint_if_available())
            self.assertIsNotNone(restored._de_state)
            self.assertEqual(int(restored._de_state["runtime_state"]["next_civ"]), 2)

            with patch("jarvishep.Sampling.diver.Sample", _FakeDiverSample):
                restored.run_nested()

            self.assertIsNotNone(restored._de_result)
            self.assertIsNotNone(restored._best_params)
            self.assertGreaterEqual(restored._de_result.best_generation, 0)

    def test_dynesty_sampler_checkpoint_payload_roundtrip(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-dynesty-") as tmp:
            sampler = Dynesty()
            sampler.logger = _NoopLogger()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_config(
                _sanitize_config(
                    {
                        "Sampling": {
                            "Method": "Dynesty",
                            "Bounds": {"nlive": 5, "rseed": 3, "run_nested": {"maxiter": 1}},
                            "Variables": [
                                {
                                    "name": "x",
                                    "description": "x",
                                    "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                                }
                            ],
                        }
                    }
                )
            )
            sampler.sampler = type(
                "DummyDynestySampler",
                (),
                {"results": {"samples_uid": ["u0"], "logl": [0.1], "logvol": [0.2]}},
            )()
            sampler.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            sampler.set_runtime_checkpoint_context(run_spec={"normalized_config": sampler.config}, factory_blueprint={})
            payload = sampler.export_runtime_state()
            restored = Dynesty()
            restored.logger = _NoopLogger()
            restored.info["sample"] = sampler.info["sample"]
            restored.set_config(deepcopy(sampler.config))
            restored.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            restored.import_runtime_state(payload)
            self.assertEqual(restored._nlive, 5)
            self.assertEqual(restored._execution_profile, {})
            self.assertEqual(restored._sampler_results_snapshot["samples_uid"], ["u0"])

    def test_multinest_sampler_checkpoint_payload_roundtrip(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-multinest-") as tmp:
            sampler = MultiNest()
            sampler.logger = _NoopLogger()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_config(
                _sanitize_config(
                    {
                        "Sampling": {
                            "Method": "MultiNest",
                            "Bounds": {"nlive": 5, "rseed": 3, "run_nested": {"maxiter": 1}},
                            "Variables": [
                                {
                                    "name": "x",
                                    "description": "x",
                                    "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                                }
                            ],
                        }
                    }
                )
            )
            sampler.sampler = type(
                "DummyMultiNestSampler",
                (),
                {"results": {"samples_uid": ["u0"], "logl": [0.1], "logvol": [0.2], "samples": np.zeros((1, 1)), "samples_u": np.zeros((1, 1)), "logwt": np.array([0.0]), "logz": np.array([0.0]), "logzerr": np.array([0.0]), "samples_n": np.array([1]), "ncall": np.array([1]), "samples_it": np.array([0]), "samples_id": np.array([0]), "information": np.array([0.0])}},
            )()
            sampler.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            sampler.set_runtime_checkpoint_context(run_spec={"normalized_config": sampler.config}, factory_blueprint={})
            payload = sampler.export_runtime_state()
            restored = MultiNest()
            restored.logger = _NoopLogger()
            restored.info["sample"] = sampler.info["sample"]
            restored.set_config(deepcopy(sampler.config))
            restored.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            restored.import_runtime_state(payload)
            self.assertEqual(restored._nlive, 5)
            self.assertEqual(restored._execution_profile, {})
            self.assertEqual(restored._sampler_results_snapshot["samples_uid"], ["u0"])


if __name__ == "__main__":
    unittest.main()
