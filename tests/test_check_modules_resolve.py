#!/usr/bin/env python3
"""check-modules CSV resolution + sampler fallback (D13.9 follow-up)."""

from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock

import numpy as np

from jarvishep2.core import Jarvis2Core
from jarvishep2.sample import Sample
from jarvishep2.task_config import (
    check_modules_n_samples,
    check_modules_timeout_sec,
    get_check_modules_settings,
    resolve_check_modules_csv,
)


class CheckModulesResolveTests(unittest.TestCase):
    def test_settings_use_envreqs_v2_data_only(self) -> None:
        cfg = {
            "Sampling": {"data": "&J/data/mine.csv"},
            "EnvReqs": {"V2": {"check_modules": {"data": "&J/data/default.csv", "n_samples": 7}}},
        }
        settings = get_check_modules_settings(cfg)
        self.assertEqual(settings["data"], "&J/data/default.csv")
        self.assertEqual(settings["n_samples"], 7)

    def test_settings_default_n_samples(self) -> None:
        settings = get_check_modules_settings({})
        self.assertEqual(settings["n_samples"], 10)
        self.assertEqual(settings["timeout_sec"], 120.0)
        self.assertIn("check_modules_points.csv", settings["data"])

    def test_timeout_sec_from_envreqs(self) -> None:
        cfg = {"EnvReqs": {"V2": {"check_modules": {"timeout_sec": 30}}}}
        self.assertEqual(get_check_modules_settings(cfg)["timeout_sec"], 30.0)
        self.assertEqual(check_modules_timeout_sec(cfg), 30.0)

    def test_timeout_alias_is_removed_and_invalid_fallback(self) -> None:
        retired = {"EnvReqs": {"V2": {"check_modules": {"timeout": 45.5}}}}
        self.assertEqual(check_modules_timeout_sec(retired), 120.0)
        bad = {"EnvReqs": {"V2": {"check_modules": {"timeout_sec": -1}}}}
        self.assertEqual(check_modules_timeout_sec(bad), 120.0)
        self.assertEqual(check_modules_timeout_sec({}), 120.0)

    def test_resolve_csv_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "pts.csv")
            with open(csv_path, "w", encoding="utf-8") as handle:
                handle.write("uuid,x,y\na,0.1,0.2\n")
            cfg = {
                "project_root": tmp,
                "task_root": tmp,
                "task_yaml": os.path.join(tmp, "task.yaml"),
                "EnvReqs": {"V2": {"check_modules": {"data": "pts.csv"}}},
            }
            path, raw = resolve_check_modules_csv(cfg)
            self.assertEqual(raw, "pts.csv")
            self.assertEqual(path, csv_path)

    def test_resolve_csv_missing_returns_none(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cfg = {
                "project_root": tmp,
                "task_root": tmp,
                "task_yaml": os.path.join(tmp, "task.yaml"),
                "EnvReqs": {"V2": {"check_modules": {"data": "missing.csv"}}},
            }
            path, raw = resolve_check_modules_csv(cfg)
            self.assertIsNone(path)
            self.assertEqual(raw, "missing.csv")

    def test_n_samples_from_envreqs(self) -> None:
        cfg = {"EnvReqs": {"V2": {"check_modules": {"n_samples": 3}}}}
        self.assertEqual(check_modules_n_samples(cfg), 3)

    def test_run_check_modules_uses_config_timeout(self) -> None:
        core = Jarvis2Core(
            {"EnvReqs": {"V2": {"check_modules": {"timeout_sec": 12.5}}}}
        )
        core._resolve_sample_root = mock.Mock(return_value="/tmp/sample")  # type: ignore[method-assign]
        core._build_check_module_samples = mock.Mock(return_value=[Sample(uuid="u1")])  # type: ignore[method-assign]
        core.submit_samples = mock.Mock()  # type: ignore[method-assign]
        core.wait_for_results = mock.Mock()  # type: ignore[method-assign]
        core._finalize_sample_buckets = mock.Mock()  # type: ignore[method-assign]
        count = core.run_check_modules()
        self.assertEqual(count, 1)
        core.wait_for_results.assert_called_once_with(
            1,
            timeout=12.5,
            progress_total=1,
            progress_base=0,
            require_worker_completion=True,
        )

    def test_run_check_modules_cli_timeout_overrides_config(self) -> None:
        core = Jarvis2Core(
            {"EnvReqs": {"V2": {"check_modules": {"timeout_sec": 12.5}}}}
        )
        core._resolve_sample_root = mock.Mock(return_value="/tmp/sample")  # type: ignore[method-assign]
        core._build_check_module_samples = mock.Mock(return_value=[Sample(uuid="u1")])  # type: ignore[method-assign]
        core.submit_samples = mock.Mock()  # type: ignore[method-assign]
        core.wait_for_results = mock.Mock()  # type: ignore[method-assign]
        core._finalize_sample_buckets = mock.Mock()  # type: ignore[method-assign]
        core.run_check_modules(timeout=3.0)
        core.wait_for_results.assert_called_once_with(
            1,
            timeout=3.0,
            progress_total=1,
            progress_base=0,
            require_worker_completion=True,
        )


class CheckModulesBuildSamplesTests(unittest.TestCase):
    def test_csv_path_used_when_present(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "pts.csv")
            with open(csv_path, "w", encoding="utf-8") as handle:
                handle.write("uuid,x,y\np1,0.1,0.2\np2,0.3,0.4\n")
            core = Jarvis2Core(
                {
                    "project_root": tmp,
                    "task_root": tmp,
                    "task_yaml": os.path.join(tmp, "task.yaml"),
                    "Sampling": {
                        "Method": "Random",
                        "Bounds": {"point_number": 100},
                        "Variables": [
                            {
                                "name": "x",
                                "description": "x",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0.0, "max": 1.0},
                                },
                            },
                            {
                                "name": "y",
                                "description": "y",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0.0, "max": 1.0},
                                },
                            },
                        ],
                    },
                    "EnvReqs": {"V2": {"check_modules": {"data": "pts.csv"}}},
                }
            )
            from jarvishep2.Sampling.randoms import RandomS

            sampler = RandomS()
            sampler.set_config(core.config)
            core.sampler = sampler
            core._logger = mock.Mock()
            samples = core._build_check_module_samples()
            self.assertEqual(len(samples), 2)
            self.assertEqual(samples[0].uuid, "p1")
            core._logger.warning.assert_any_call(
                mock.ANY,  # first format string for CSV path
                mock.ANY,
                mock.ANY,
            )

    def test_sampler_fallback_when_csv_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            core = Jarvis2Core(
                {
                    "project_root": tmp,
                    "task_root": tmp,
                    "task_yaml": os.path.join(tmp, "task.yaml"),
                    "Sampling": {
                        "Method": "Random",
                        "Bounds": {"point_number": 100, "seed": 1},
                        "Variables": [
                            {
                                "name": "x",
                                "description": "x",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0.0, "max": 1.0},
                                },
                            },
                            {
                                "name": "y",
                                "description": "y",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0.0, "max": 1.0},
                                },
                            },
                        ],
                    },
                    "EnvReqs": {"V2": {"check_modules": {"data": "no_such.csv", "n_samples": 3}}},
                }
            )
            from jarvishep2.Sampling.randoms import RandomS

            sampler = RandomS()
            sampler.set_config(core.config)
            core.sampler = sampler
            core._logger = mock.Mock()
            samples = core._build_check_module_samples()
            self.assertEqual(len(samples), 3)
            for sample in samples:
                self.assertIsInstance(sample, Sample)
                self.assertEqual(len(sample.u_coords), 2)

    def test_unit_cube_fallback_without_propose_next(self) -> None:
        """Dynesty/MultiNest-like: no propose_next → random unit-cube draws."""
        with tempfile.TemporaryDirectory() as tmp:
            core = Jarvis2Core(
                {
                    "project_root": tmp,
                    "task_root": tmp,
                    "task_yaml": os.path.join(tmp, "task.yaml"),
                    "Sampling": {
                        "Method": "Dynesty",
                        "Bounds": {"seed": 7},
                        "Variables": [
                            {
                                "name": "xx",
                                "description": "x",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0.0, "max": 5.0},
                                },
                            },
                            {
                                "name": "yy",
                                "description": "y",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0.0, "max": 5.0},
                                },
                            },
                        ],
                    },
                    "EnvReqs": {"V2": {"check_modules": {"n_samples": 4}}},
                }
            )
            from jarvishep2.Sampling.sampler import SamplingVirtial

            sampler = SamplingVirtial()
            sampler.set_config(core.config)
            core.sampler = sampler
            core._logger = mock.Mock()
            samples = core._build_check_module_samples()
            self.assertEqual(len(samples), 4)
            for sample in samples:
                self.assertEqual(len(np.asarray(sample.u_coords)), 2)


if __name__ == "__main__":
    unittest.main()
