#!/usr/bin/env python3
"""WP-D3.1 registered_executables + two-phase CommandParser tests."""

from __future__ import annotations

import os
import pickle
import tempfile
import unittest

from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.command_parser import CommandParser, prepare_calculator_modules
from jarvishep2.core import Jarvis2Core
from jarvishep2.factory import TaskFactory
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.sample import Sample
from jarvishep2.worker_config import build_worker_config

from test_worker_calculator import (
    EGGBOX_CALC_MODULE,
    EGGBOX_SCRIPT,
    FIXTURES,
    LIKELIHOOD_EXPRESSIONS,
    _load_csv_points,
    _normalize_database_records,
    _start_tcp_fakeredis,
    _worker_config,
)

TESTS_ROOT = os.path.dirname(__file__)
SAFE_SOURCE = os.path.join(TESTS_ROOT, "fixtures", "safe_calc", "source")


def _collect_command_strings(config: dict) -> list[str]:
    strings: list[str] = []

    def _walk(value) -> None:
        if isinstance(value, str):
            strings.append(value)
        elif isinstance(value, dict):
            for item in value.values():
                _walk(item)
        elif isinstance(value, list):
            for item in value:
                _walk(item)

    _walk(config)
    return strings


class CommandParserUnitTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.mkdtemp()
        self.project_root = self.tmpdir
        self.libdeps_root = os.path.join(self.tmpdir, "deps")
        os.makedirs(self.libdeps_root, exist_ok=True)
        self.config = {
            "task_result_dir": self.tmpdir,
            "Scan": {"name": "eggbox-scan"},
            "LibDeps": {
                "path": self.libdeps_root,
                "Modules": [
                    {
                        "name": "EggBoxSafe",
                        "path": os.path.dirname(EGGBOX_SCRIPT),
                    }
                ],
                "registered_executables": [
                    {
                        "name": "eggboxlk",
                        "source": EGGBOX_SCRIPT,
                        "resolution": "direct_path",
                    },
                    {
                        "name": "safetool",
                        "source": SAFE_SOURCE,
                        "resolution": "symlink",
                    },
                ],
            },
        }
        self.parser = CommandParser.from_config(self.config, project_root=self.project_root)

    def test_phase1_removes_static_tokens_from_calculator_config(self) -> None:
        raw_module = {
            "name": "EggBox",
            "path": "&J/calculators/runtime",
            "source": "${LibDeps:EggBoxSafe}",
            "execution": {
                "path": "&J/calculators/runtime",
                "commands": [{"cmd": "python3 eggboxlk @Sdir/input.json", "cwd": "@Sdir"}],
                "input": [{"path": "@Sdir/input.json", "type": "JSON"}],
                "output": [{"path": "@Sdir/output.json", "type": "JSON"}],
            },
        }
        resolved = self.parser.resolve_static_config(raw_module)
        for text in _collect_command_strings(resolved):
            self.assertFalse(self.parser.has_static_tokens(text))
            if "@Sdir" in text or "@SampleID" in text or "@PackID" in text:
                continue
            self.assertNotIn("&J", text)
            self.assertNotIn("${LibDeps:", text)
        self.assertIn(EGGBOX_SCRIPT, resolved["execution"]["commands"][0]["cmd"])
        self.assertIn("@Sdir", resolved["execution"]["commands"][0]["cmd"])

    def test_registered_executable_direct_path_and_symlink(self) -> None:
        direct = self.parser.registered["eggboxlk"]
        self.assertEqual(direct.resolution, "direct_path")
        self.assertEqual(os.path.realpath(direct.path), os.path.realpath(EGGBOX_SCRIPT))

        symlink = self.parser.registered["safetool"]
        self.assertEqual(symlink.resolution, "symlink")
        self.assertTrue(os.path.islink(symlink.path))
        self.assertEqual(os.path.realpath(symlink.path), os.path.realpath(SAFE_SOURCE))

    def test_phase2_resolves_sample_tokens(self) -> None:
        sample = Sample.from_params({"x": 0.1, "y": 0.2, "uuid": "phase2-sample"})
        sample.set_config(
            {
                "sample_dirs": self.tmpdir,
                "task_result_dir": self.tmpdir,
                "sample_artifacts": "always",
                "workflow_has_calculator": True,
                "workflow_references_sdir": True,
            }
        )
        sample.materialize()
        resolved = self.parser.resolve_sample(
            "@Sdir/input.json @SampleID @PackID",
            sample_info=sample.info,
            pack_id="pack-123",
        )
        self.assertIn(str(sample.save_dir), resolved)
        self.assertIn(str(sample.uuid), resolved)
        self.assertIn("pack-123", resolved)
        self.assertNotIn("@Sdir", resolved)

    def test_phase2_rejects_leftover_static_tokens(self) -> None:
        sample = Sample.from_params({"uuid": "bad-sample"})
        with self.assertRaises(ValueError):
            self.parser.resolve_sample(
                "python3 &J/deps/tool.py",
                sample_info=sample.info,
            )

    def test_phase2_matches_calculator_runtime_resolver(self) -> None:
        sample = Sample.from_params({"x": 0.25, "y": 0.25, "uuid": "parity-sample"})
        sample.set_config(
            {
                "sample_dirs": self.tmpdir,
                "task_result_dir": self.tmpdir,
                "sample_artifacts": "always",
                "workflow_has_calculator": True,
                "workflow_references_sdir": True,
            }
        )
        sample.materialize()
        module = CalculatorModule("EggBox", EGGBOX_CALC_MODULE)
        module.sample_info = dict(sample.info)
        module.acquire_pack_id("pack-parity")
        legacy = module._resolve_runtime_tokens(
            "@Sdir/input.json @SampleID",
            stage="execution",
            field="path",
        )
        module.attach_command_parser(self.parser)
        module.sample_info = dict(sample.info)
        modern = module._resolve_runtime_tokens(
            "@Sdir/input.json @SampleID",
            stage="execution",
            field="path",
        )
        self.assertEqual(modern, legacy)

    def test_phase1_config_pickles_under_spawn(self) -> None:
        resolved = prepare_calculator_modules([EGGBOX_CALC_MODULE], self.parser)
        blob = pickle.dumps(resolved, protocol=pickle.HIGHEST_PROTOCOL)
        restored = pickle.loads(blob)
        self.assertEqual(restored, resolved)
        self.assertEqual(get_spawn_context().get_start_method(), "spawn")


class CommandParserIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_registered_executable_eggbox_database_parity(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(
                    os.path.join(FIXTURES, "expected_calculator_records.json"),
                    encoding="utf-8",
                ) as handle:
                    import json

                    expected_records = json.load(handle)

                registered_module = {
                    "name": "EggBox",
                    "required_modules": [],
                    "clone_shadow": False,
                    "path": os.path.dirname(EGGBOX_SCRIPT),
                    "execution": {
                        "path": os.path.dirname(EGGBOX_SCRIPT),
                        "commands": [{"cmd": f"python3 eggboxlk", "cwd": "@Sdir"}],
                        "input": EGGBOX_CALC_MODULE["execution"]["input"],
                        "output": EGGBOX_CALC_MODULE["execution"]["output"],
                    },
                }
                task_config = {
                    "Runtime": {"mode": "redis", "workers": 1, "redis": redis_config},
                    "task_result_dir": tmpdir,
                    "Scan": {"name": "registered-eggbox"},
                    "LibDeps": {
                        "registered_executables": [
                            {
                                "name": "eggboxlk",
                                "source": EGGBOX_SCRIPT,
                                "resolution": "direct_path",
                            }
                        ]
                    },
                    "Likelihood": {"expressions": LIKELIHOOD_EXPRESSIONS},
                }
                core = Jarvis2Core(task_config)
                core.init_redis()
                core.init_command_parser()
                worker_config = core.build_worker_config(
                    calculator_modules=[registered_module],
                    calculator_pools={"EggBox": 1},
                )
                core.init_factory(worker_config)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=worker_config["calculator_modules"],
                    include_likelihood=True,
                )
                core.set_sampler(sampler)

                samples = []
                for row in _load_csv_points():
                    import numpy as np

                    sample = sampler._build_sample(
                        np.array([float(row["x"]), float(row["y"])], dtype=np.float64)
                    )
                    sample.uuid = str(row["uuid"])
                    samples.append(sample)

                core.submit_samples(samples)
                core.wait_for_results(len(samples), timeout=90.0)
                core.shutdown()

                from jarvishep2.database import SimpleHDF5Writer

                records = _normalize_database_records(
                    SimpleHDF5Writer(os.path.join(tmpdir, "DATABASE", "samples.hdf5")).read_records()
                )
                self.assertEqual(records, _normalize_database_records(expected_records))
        finally:
            server.shutdown()
            server.server_close()

    def test_build_worker_config_preserves_legacy_explicit_modules(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            worker_config = build_worker_config(
                {"task_result_dir": tmpdir},
                task_result_dir=tmpdir,
                calculator_modules=[EGGBOX_CALC_MODULE],
                likelihood_expressions=LIKELIHOOD_EXPRESSIONS,
                extra={"calculator_pools": {"EggBox": 1}},
            )
            cmd = worker_config["calculator_modules"][0]["execution"]["commands"][0]["cmd"]
            self.assertIn("eggbox.py", cmd)
            self.assertIn("command_parser", worker_config)
            self.assertEqual(worker_config["calculator_pools"], {"EggBox": 1})


if __name__ == "__main__":
    unittest.main()