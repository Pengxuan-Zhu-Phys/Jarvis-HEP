#!/usr/bin/env python3
"""Check-modules contracts: CSV resolution, coordinate columns, fallback and pool policy."""

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
from jarvishep2.testing.check_modules import (
    build_check_module_u_coords,
    coordinate_columns_for_check_modules,
)
from jarvishep2.calculator_pools import resolve_calculator_pools
from jarvishep2.worker_config import build_worker_config


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


class CheckModulesIdentityContinuationTests(unittest.TestCase):
    def test_colliding_indices_are_shifted_to_append(self) -> None:
        core = Jarvis2Core()
        core._persisted_uuids = set()
        samples = [
            Sample(uuid="bridson-0", sample_index=0, u_coords=np.array([0.1, 0.2])),
            Sample(uuid="bridson-1", sample_index=1, u_coords=np.array([0.3, 0.4])),
        ]
        out = core._get_scan()._ensure_unique_check_identities(samples, 10)
        self.assertEqual([item.sample_index for item in out], [10, 11])
        self.assertNotEqual(out[0].uuid, "bridson-0")
        self.assertNotEqual(out[1].uuid, "bridson-1")

    def test_already_continued_indices_are_kept(self) -> None:
        core = Jarvis2Core()
        core._persisted_uuids = set()
        samples = [
            Sample(uuid="bridson-10", sample_index=10, u_coords=np.array([0.1, 0.2])),
            Sample(uuid="bridson-11", sample_index=11, u_coords=np.array([0.3, 0.4])),
        ]
        out = core._get_scan()._ensure_unique_check_identities(samples, 10)
        self.assertEqual([item.sample_index for item in out], [10, 11])
        self.assertEqual(out[0].uuid, "bridson-10")
        self.assertEqual(out[1].uuid, "bridson-11")

    def test_continue_fast_forwards_sampler_from_database_prefix(self) -> None:
        core = Jarvis2Core()
        core._persisted_index_prefix = 10
        core._persisted_records_count = 10
        core._logger = mock.Mock()
        advanced: list[int] = []

        class _Sampler:
            def advance_to_persisted_prefix(self, prefix: int) -> int:
                advanced.append(int(prefix))
                return int(prefix)

        core.sampler = _Sampler()  # type: ignore[assignment]
        start = core._get_scan()._continue_check_sampler_from_database()
        self.assertEqual(start, 10)
        self.assertEqual(advanced, [10])


class CheckModulesColumnTests(unittest.TestCase):
    def test_columns_from_sampling_variables(self) -> None:
        config = {
            "Sampling": {
                "Variables": [
                    {"name": "m0", "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}}},
                    {"name": "m12", "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}}},
                ]
            }
        }
        cols = coordinate_columns_for_check_modules(config, ["uuid", "m0", "m12", "extra"])
        self.assertEqual(cols, ["m0", "m12"])

    def test_columns_from_csv_headers_when_no_variables(self) -> None:
        cols = coordinate_columns_for_check_modules({}, ["uuid", "a", "b"])
        self.assertEqual(cols, ["a", "b"])

    def test_build_u_coords_order(self) -> None:
        row = {"m0": "0.1", "m12": "0.2", "uuid": "u1"}
        coords = build_check_module_u_coords(row, ["m0", "m12"])
        self.assertEqual(coords, [0.1, 0.2])


class CheckModulesPoolPolicyTests(unittest.TestCase):
    def test_check_policy_pins_make_parallel_and_pools_to_one(self) -> None:
        core = Jarvis2Core.__new__(Jarvis2Core)
        core.config = {
            "Runtime": {"workers": 8},
            "EnvReqs": {
                "V2": {
                    "workers": 8,
                    "sample_directory": {"pack": True, "enabled": True},
                    "archiver": {"pack_buckets": True},
                }
            },
            "Scan": {"name": "smoke"},
            "Calculators": {
                "make_parallel": 8,
                "Pools": {"MicroOMEGAs_Vector": 8},
                "Modules": [
                    {
                        "name": "MicroOMEGAs_Vector",
                        "make_parallel": 8,
                        "path": "&J/calculators/microOMEGAs/@PackID",
                    }
                ],
            },
            "Operas": {"Modules": []},
        }
        core.runtime = {}
        core.info = {}
        with tempfile.TemporaryDirectory() as tmpdir:
            core.info["task_result_dir"] = tmpdir
            core.config["task_result_dir"] = tmpdir
            core._apply_check_modules_runtime_policy()

            calc = core.config["Calculators"]
            self.assertEqual(calc["make_parallel"], 1)
            self.assertEqual(calc["Modules"][0]["make_parallel"], 1)
            self.assertEqual(calc["Pools"]["MicroOMEGAs_Vector"], 1)
            self.assertNotIn("Archiver", calc)
            self.assertEqual(core.config["Runtime"]["workers"], 1)
            self.assertEqual(core.config["EnvReqs"]["V2"]["workers"], 1)
            self.assertFalse(core.config["EnvReqs"]["V2"]["archiver"]["pack_buckets"])

            worker = build_worker_config(core.config, task_result_dir=tmpdir)
            pools = resolve_calculator_pools(worker)
            self.assertEqual(pools.get("MicroOMEGAs_Vector"), 1)
            self.assertTrue(
                os.path.isdir(os.path.join(tmpdir, "SAMPLE", "test"))
                or core.info.get("sample_root", "").endswith(
                    os.path.join("SAMPLE", "test")
                )
            )
            self.assertTrue(os.path.isdir(os.path.join(tmpdir, "DATABASE", "test")))
            self.assertEqual(
                core.info.get("database_dir"),
                os.path.join(tmpdir, "DATABASE", "test"),
            )

    def test_check_policy_keeps_smoke_workspace_and_full_scan_database(self) -> None:
        core = Jarvis2Core.__new__(Jarvis2Core)
        core.config = {"Runtime": {}, "EnvReqs": {"V2": {}}, "Scan": {"name": "smoke"}}
        core.runtime = {}
        core.info = {}
        with tempfile.TemporaryDirectory() as tmpdir:
            core.info["task_result_dir"] = tmpdir
            core.config["task_result_dir"] = tmpdir
            full_db = os.path.join(tmpdir, "DATABASE")
            os.makedirs(full_db)
            full_hdf5 = os.path.join(full_db, "samples.hdf5")
            with open(full_hdf5, "w", encoding="utf-8") as handle:
                handle.write("full-scan")
            smoke_hdf5 = os.path.join(full_db, "test", "samples.hdf5")
            os.makedirs(os.path.dirname(smoke_hdf5))
            with open(smoke_hdf5, "w", encoding="utf-8") as handle:
                handle.write("previous-check")
            leftover_sample = os.path.join(tmpdir, "SAMPLE", "test", "old-uuid")
            os.makedirs(leftover_sample)
            with open(os.path.join(leftover_sample, "output.json"), "w", encoding="utf-8") as handle:
                handle.write("{}")

            core._apply_check_modules_runtime_policy()

            self.assertTrue(os.path.isfile(full_hdf5))
            with open(full_hdf5, encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "full-scan")
            with open(smoke_hdf5, encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "previous-check")
            self.assertTrue(os.path.isfile(os.path.join(leftover_sample, "output.json")))
            self.assertEqual(
                core.info.get("database_dir"),
                os.path.join(tmpdir, "DATABASE", "test"),
            )


if __name__ == "__main__":
    unittest.main()
