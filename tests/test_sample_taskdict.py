#!/usr/bin/env python3
"""Unit tests for jarvishep2.sample (WP-D0.1)."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from uuid import uuid4

import numpy as np

from jarvishep2.runtime_config import should_eager_materialize
from jarvishep2.sample import (
    ExecutionStep,
    Sample,
    ensure_sample_materialized,
    materialize_failure_artifacts,
)
from jarvishep2.sample_logger import BufferedSampleLogger


class SampleTaskDictTests(unittest.TestCase):
    def test_round_trip_preserves_light_fields(self):
        original = Sample(
            uuid=str(uuid4()),
            u_coords=np.array([0.1, 0.2, 0.3], dtype=np.float64),
            execution_plan=[
                ExecutionStep(type="calculator", name="DemoCalc", layer=0),
                ExecutionStep(type="opera", name="LogL", layer=1, params={"x": 1}),
            ],
            opera_params={"shift": 0.5},
            sample_artifacts="auto",
            priority=2,
            params={"hidden": 99.0},
            info={"secret": True},
            observables={"z": 1.0},
        )
        wire = original.to_task_dict()
        rebuilt = Sample.from_task_dict(wire)

        self.assertEqual(rebuilt.uuid, original.uuid)
        self.assertTrue(np.allclose(rebuilt.u_coords, original.u_coords))
        self.assertEqual(len(rebuilt.execution_plan), 2)
        self.assertEqual(rebuilt.execution_plan[0].name, "DemoCalc")
        self.assertEqual(rebuilt.execution_plan[1].params, {"x": 1})
        self.assertEqual(rebuilt.opera_params, {"shift": 0.5})
        self.assertEqual(rebuilt.sample_artifacts, "auto")
        self.assertEqual(rebuilt.priority, 2)

        self.assertEqual(rebuilt.params, {})
        self.assertEqual(rebuilt.info, {})
        self.assertEqual(rebuilt.observables, {})
        self.assertFalse(rebuilt._materialized)
        self.assertIsNone(rebuilt._logger)

    def test_to_task_dict_is_json_serializable_without_logger(self):
        sample = Sample.from_params({"x": 1.0, "y": 2.0})
        with tempfile.TemporaryDirectory() as tmpdir:
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            self.assertIsInstance(sample.info.get("logger"), BufferedSampleLogger)

            wire = sample.to_task_dict()
            self.assertNotIn("logger", wire)
            json.dumps(wire)

    def test_lazy_operas_only_defers_directory_creation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            self.assertFalse(sample.info["_materialized"])
            self.assertIsNone(sample.info["save_dir"])
            self.assertFalse(os.path.exists(os.path.join(tmpdir, sample.uuid)))

    def test_always_mode_materializes_immediately(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                }
            )
            self.assertTrue(sample.info["_materialized"])
            self.assertTrue(os.path.isdir(sample.info["save_dir"]))
            self.assertTrue(os.path.exists(sample.info["run_log"]))

    def test_never_mode_skips_failure_materialization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_info = {
                "uuid": "dead-beef",
                "params": {"x": 1.0},
                "sample_dirs": tmpdir,
                "task_result_dir": tmpdir,
                "sample_artifacts": "never",
                "_materialized": False,
            }
            result = materialize_failure_artifacts(sample_info, error="boom")
            self.assertIsNone(result)
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "dead-beef")))

    def test_failure_materialization_writes_minimal_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_info = {
                "uuid": "dead-beef",
                "params": {"x": 1.0},
                "sample_dirs": tmpdir,
                "task_result_dir": tmpdir,
                "sample_artifacts": "auto",
                "workflow_has_calculator": False,
                "workflow_references_sdir": False,
                "_materialized": False,
            }
            save_dir = materialize_failure_artifacts(sample_info, error="boom")
            self.assertIsNotNone(save_dir)
            log_path = os.path.join(save_dir, "Sample_running.log")
            self.assertTrue(os.path.exists(log_path))
            with open(log_path, "r", encoding="utf-8") as handle:
                self.assertIn("Sample failed -> boom", handle.read())

    def test_sdir_resolution_triggers_materialization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            resolved = sample.resolve_token("echo @Sdir/out.txt", stage="runtime", field="cmd")
            self.assertTrue(resolved.endswith("/out.txt"))
            self.assertTrue(os.path.isdir(sample.info["save_dir"]))

    def test_resolve_pack_and_sample_tokens(self):
        sample = Sample(uuid="abc-123", u_coords=np.array([0.5]))
        sample.info = {"uuid": "abc-123", "pack_id": "pack-9"}
        resolved = sample.resolve_token("@SampleID/@PackID/log.txt")
        self.assertEqual(resolved, "abc-123/pack-9/log.txt")

    def test_status_transitions(self):
        sample = Sample.from_params({"x": 1.0})
        self.assertEqual(sample.status, "Created")
        sample.start()
        self.assertEqual(sample.status, "Running")
        self.assertEqual(sample.info["status"], "Running")
        sample.status = "Completed"
        self.assertEqual(sample.to_info_dict()["status"], "Completed")

    def test_to_info_dict_excludes_logger(self):
        sample = Sample.from_params({"x": 1.0})
        sample.info = {"logger": BufferedSampleLogger(extra={"module": "x"}), "uuid": sample.uuid}
        info = sample.to_info_dict()
        self.assertNotIn("logger", info)

    def test_calculator_workflow_is_eager_in_auto_mode(self):
        self.assertTrue(
            should_eager_materialize(
                {
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": True,
                    "workflow_references_sdir": False,
                }
            )
        )

    def test_ensure_sample_materialized_from_info_dict(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_info = {
                "uuid": "sdir-case",
                "params": {"x": 1.0},
                "sample_dirs": tmpdir,
                "task_result_dir": tmpdir,
                "sample_artifacts": "auto",
                "_materialized": False,
            }
            save_dir = ensure_sample_materialized(sample_info)
            self.assertIsNotNone(save_dir)
            self.assertTrue(os.path.isdir(save_dir))


if __name__ == "__main__":
    unittest.main()