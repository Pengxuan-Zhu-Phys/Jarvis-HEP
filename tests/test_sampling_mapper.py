#!/usr/bin/env python3
"""D22 Sampling.Mapper — list of {name, expression}, closed namespace, bilateral map."""

from __future__ import annotations

import math
import pickle
import unittest
from unittest import mock

import numpy as np

from jarvishep2.mapper import (
    MapperError,
    MapperPipeline,
    build_mapper,
    build_mapper_spec_from_config,
    mapper_block_fingerprint,
)
from jarvishep2.plot_scene import resolve_plot_axis_keys
from jarvishep2.Sampling.runtime_checkpoint import (
    build_payload,
    build_run_spec,
    check_mapper_fingerprint,
    check_operas_constants_fingerprint,
    operas_constants_fingerprint,
    stable_json_hash,
)
from jarvishep2.operas_functions import operas_expression_references_required
from jarvishep2.Sampling.sampling_utils import map_u_to_physical, physical_from_u
from jarvishep2.Sampling.variables import load_variables
from jarvishep2.task_schema import _schema_catalog, task_card_validator
from jarvishep2.task_validation import validate_task_config
from jarvishep2.worker_config import build_worker_config


def _clear_schema_cache() -> None:
    _schema_catalog.cache_clear()
    task_card_validator.cache_clear()


def _mapper_list(*pairs: tuple[str, str]) -> list[dict[str, str]]:
    return [{"name": name, "expression": expr} for name, expr in pairs]


def _base_card(*, method: str = "Bridson", with_mapper: bool = False) -> dict:
    sampling: dict = {
        "Method": method,
        "Bounds": {"radius": 0.05, "max_attempt": 10, "seed": 0},
        "Variables": [
            {
                "name": "t",
                "distribution": {
                    "type": "Flat",
                    "parameters": {"min": 0.0, "max": 6.283185307, "length": 1.0},
                },
            }
        ],
    }
    if method == "Random":
        sampling["Bounds"] = {"point_number": 20, "seed": 0}
    if method == "Dynesty":
        sampling["Bounds"] = {"nlive": 20, "dlogz": 0.5, "seed": 0}
    if with_mapper:
        sampling["Mapper"] = _mapper_list(("x", "cos(t)"), ("y", "sin(t)"))
        sampling["selection"] = "y > 0"
    return {
        "Scan": {"name": "mapper-unit"},
        "Sampling": sampling,
        "EnvReqs": {"V2": {"workers": 1}},
        "Calculators": {},
    }


class MapperPipelineUnitTests(unittest.TestCase):
    def test_variables_only_matches_legacy_distribution_map(self) -> None:
        card = _base_card(method="Random", with_mapper=False)
        pipeline = MapperPipeline.from_config(card)
        vars_ = load_variables(card)
        u = np.array([0.37])
        self.assertEqual(pipeline.map(u), map_u_to_physical(u, vars_))
        self.assertEqual(pipeline.derive_names, ())
        self.assertEqual(pipeline.output_names, ("t",))

    def test_mapper_unit_circle(self) -> None:
        card = _base_card(with_mapper=True)
        pipeline = MapperPipeline.from_config(card)
        # u=0.25 → t ≈ π/2 → (x,y) ≈ (0, 1)  (card max is truncated 2π)
        mapped = pipeline.map(np.array([0.25]))
        self.assertAlmostEqual(mapped["t"], math.pi / 2, places=9)
        self.assertAlmostEqual(mapped["x"], 0.0, places=9)
        self.assertAlmostEqual(mapped["y"], 1.0, places=9)
        self.assertAlmostEqual(mapped["x"] ** 2 + mapped["y"] ** 2, 1.0, places=12)
        self.assertEqual(pipeline.output_names, ("t", "x", "y"))

    def test_mapper_dag_order_independent_of_yaml_order(self) -> None:
        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = _mapper_list(
            ("z", "x + y"),
            ("x", "cos(t)"),
            ("y", "sin(t)"),
        )
        pipeline = MapperPipeline.from_config(card)
        mapped = pipeline.map(np.array([0.0]))  # t=0 → x=1, y=0, z=1
        self.assertAlmostEqual(mapped["x"], 1.0, places=12)
        self.assertAlmostEqual(mapped["y"], 0.0, places=12)
        self.assertAlmostEqual(mapped["z"], 1.0, places=12)
        self.assertEqual(set(pipeline.derive_names), {"x", "y", "z"})

    def test_closed_namespace_rejects_observable_symbol(self) -> None:
        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = _mapper_list(("x", "sin(t) + LogL"))
        with self.assertRaises(MapperError) as ctx:
            build_mapper_spec_from_config(card)
        self.assertEqual(ctx.exception.code, "JV2-MAP-002")

    def test_legacy_name_map_rejected(self) -> None:
        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = {"x": "cos(t)", "y": "sin(t)"}
        with self.assertRaises(MapperError) as ctx:
            build_mapper_spec_from_config(card)
        self.assertEqual(ctx.exception.code, "JV2-MAP-001")
        self.assertIn("list", ctx.exception.message)

    def test_m5_rejects_wrong_u_length(self) -> None:
        card = _base_card(method="Random")
        pipeline = MapperPipeline.from_config(card)
        with self.assertRaises(ValueError):
            pipeline.map(np.array([0.1, 0.2]))

    def test_control_and_worker_bilateral_identity(self) -> None:
        card = _base_card(with_mapper=True)
        control = MapperPipeline.from_config(card)
        worker_cfg = build_worker_config(card, task_result_dir="/tmp/mapper-test")
        worker = build_mapper(worker_cfg["mapper"])
        assert worker is not None
        for u in (0.0, 0.1, 0.37, 0.5, 0.99):
            coords = np.array([u])
            self.assertEqual(control.map(coords), worker.map(coords))

    def test_pipeline_spec_is_picklable(self) -> None:
        card = _base_card(with_mapper=True)
        spec = build_mapper_spec_from_config(card)
        restored = pickle.loads(pickle.dumps(spec))
        a = MapperPipeline.from_spec(spec).map(np.array([0.42]))
        b = MapperPipeline.from_spec(restored).map(np.array([0.42]))
        self.assertEqual(a, b)


class MapperValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        _clear_schema_cache()

    def test_optional_mapper_cards_validate(self) -> None:
        bare = _base_card(method="Random", with_mapper=False)
        with_mapper = _base_card(method="Random", with_mapper=True)
        self.assertTrue(validate_task_config(bare).ok)
        report = validate_task_config(with_mapper)
        self.assertTrue(report.ok, format_issues(report))

    def test_schema_rejects_legacy_name_map(self) -> None:
        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = {"x": "cos(t)"}
        report = validate_task_config(card)
        self.assertFalse(report.ok)
        # Schema type error and/or JV2-MAP-001 depending on branch order.
        codes = {i.code for i in report.errors()}
        self.assertTrue(
            "JV2-MAP-001" in codes or any(c.startswith("JV2-SCH") for c in codes),
            msg=format_issues(report),
        )

    def test_closed_namespace_at_validate(self) -> None:
        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = _mapper_list(("x", "sin(t) + LogL"))
        report = validate_task_config(card)
        self.assertFalse(report.ok)
        codes = {i.code for i in report.errors()}
        self.assertIn("JV2-MAP-002", codes)

    def test_name_conflict_with_variable(self) -> None:
        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = _mapper_list(("t", "t + 1"))
        report = validate_task_config(card)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-MAP-003" for i in report.errors()))

    def test_csv_rejects_mapper(self) -> None:
        card = {
            "Scan": {"name": "csv-mapper"},
            "Sampling": {
                "Method": "CSV",
                "Bounds": {"path": "points.csv"},
                "Mapper": _mapper_list(("x", "1")),
            },
            "EnvReqs": {"V2": {"workers": 1}},
        }
        report = validate_task_config(card)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-MAP-010" for i in report.errors()))

    def test_dynesty_dimension_expansion_warning(self) -> None:
        card = _base_card(method="Dynesty", with_mapper=True)
        report = validate_task_config(card)
        self.assertTrue(report.ok, format_issues(report))
        warnings = [i for i in report.warnings() if i.code == "JV2-MAP-050"]
        self.assertEqual(len(warnings), 1)

    def test_cycle_is_error(self) -> None:
        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = _mapper_list(
            ("a", "b + 1"),
            ("b", "a + 1"),
        )
        report = validate_task_config(card)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-MAP-004" for i in report.errors()))


class MapperCheckpointAndPlotTests(unittest.TestCase):
    def test_mapper_hash_in_checkpoint_payload(self) -> None:
        card = _base_card(with_mapper=True)
        run_spec = build_run_spec(
            config=card,
            scan_name="mapper-unit",
            task_root="/tmp",
            task_result_dir="/tmp/out",
            sampler_name="Bridson",
        )
        payload = build_payload(run_spec=run_spec, sampler_state={"index": 0})
        integrity = payload["integrity"]
        self.assertIn("mapper_hash", integrity)
        expected = stable_json_hash(mapper_block_fingerprint(card))
        self.assertEqual(integrity["mapper_hash"], expected)

    def test_resume_refuses_changed_mapper(self) -> None:
        card = _base_card(with_mapper=True)
        run_spec = build_run_spec(
            config=card,
            scan_name="mapper-unit",
            task_root="/tmp",
            task_result_dir="/tmp/out",
            sampler_name="Bridson",
        )
        payload = build_payload(run_spec=run_spec, sampler_state={"index": 0})
        ok, _ = check_mapper_fingerprint(payload, card)
        self.assertTrue(ok)
        changed = _base_card(with_mapper=True)
        changed["Sampling"]["Mapper"] = _mapper_list(("x", "sin(t)"), ("y", "sin(t)"))
        ok2, reason = check_mapper_fingerprint(payload, changed)
        self.assertFalse(ok2)
        self.assertIn("Mapper", reason)

    def test_plot_axes_prefer_mapper_names(self) -> None:
        card = _base_card(with_mapper=True)
        records = [
            {"t": 0.1, "x": 0.9, "y": 0.1, "LogL": -1.0},
            {"t": 0.2, "x": 0.8, "y": 0.2, "LogL": -2.0},
        ]
        x, y, color = resolve_plot_axis_keys(records=records, config=card)
        self.assertEqual(x, "x")
        self.assertEqual(y, "y")
        self.assertEqual(color, "LogL")


class PhysicalFromUHelperTests(unittest.TestCase):
    def test_physical_from_u_uses_pipeline(self) -> None:
        card = _base_card(with_mapper=True)
        pipeline = MapperPipeline.from_config(card)
        vars_ = load_variables(card)
        u = np.array([0.25])
        with_pipe = physical_from_u(u, vars_, pipeline)
        without = physical_from_u(u, vars_, None)
        self.assertIn("x", with_pipe)
        self.assertNotIn("x", without)


class MapperOperasIntegrationTests(unittest.TestCase):
    """D23.9: Mapper expressions may reference Operas constants/functions."""

    def test_gate_sees_mapper_expression_values(self) -> None:
        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = _mapper_list(
            ("mz", "pdg.mZ * cos(t)"),
        )
        self.assertTrue(operas_expression_references_required(card))
        # Legacy free-form map still scanned for discovery.
        legacy = _base_card(method="Random")
        legacy["Sampling"]["Mapper"] = {"mz": "pdg.mZ * cos(t)"}
        self.assertTrue(operas_expression_references_required(legacy))

    def test_mapper_can_use_pdg_constant(self) -> None:
        try:
            from jarvis_operas import build_constant_dicts
        except ImportError:
            self.skipTest("Jarvis-Operas not installed")
        if "pdg.mZ" not in build_constant_dicts():
            self.skipTest("pdg.mZ not registered")

        card = _base_card(method="Random")
        card["Sampling"]["Mapper"] = _mapper_list(
            ("mz", "pdg.mZ * cos(t)"),
        )
        report = validate_task_config(card)
        self.assertTrue(report.ok, format_issues(report))
        pipeline = MapperPipeline.from_config(card)
        mapped = pipeline.map(np.array([0.0]))  # t=0 → cos=1
        self.assertAlmostEqual(mapped["mz"], 91.1876, places=9)
        # Constant folded out of free symbols — only t is a sample var.
        self.assertEqual(set(pipeline.derive_names), {"mz"})

    def test_operas_constants_hash_in_checkpoint(self) -> None:
        card = _base_card(with_mapper=True)
        run_spec = build_run_spec(
            config=card,
            scan_name="mapper-unit",
            task_root="/tmp",
            task_result_dir="/tmp/out",
            sampler_name="Bridson",
        )
        payload = build_payload(run_spec=run_spec, sampler_state={"index": 0})
        integrity = payload["integrity"]
        self.assertIn("operas_constants_hash", integrity)
        fp = operas_constants_fingerprint()
        if fp is None:
            self.assertIsNone(integrity["operas_constants_hash"])
        else:
            expected = stable_json_hash(fp)
            self.assertEqual(integrity["operas_constants_hash"], expected)
        ok, reason = check_operas_constants_fingerprint(payload)
        self.assertTrue(ok)
        self.assertFalse(reason.startswith("skip:"))

    def test_operas_constants_hash_refuses_drift(self) -> None:
        card = _base_card(with_mapper=True)
        run_spec = build_run_spec(
            config=card,
            scan_name="mapper-unit",
            task_root="/tmp",
            task_result_dir="/tmp/out",
            sampler_name="Bridson",
        )
        payload = build_payload(run_spec=run_spec, sampler_state={"index": 0})
        # Simulate Operas drift by rewriting stored hash.
        payload["integrity"]["operas_constants_hash"] = "0" * 64
        ok, reason = check_operas_constants_fingerprint(payload)
        self.assertFalse(ok)
        self.assertIn("Operas", reason)
        self.assertNotIn("Mapper", reason)

    def test_operas_unavailable_skips_fingerprint_not_refuse(self) -> None:
        """D23.13: missing Operas must not look like constant drift."""
        card = _base_card(with_mapper=True)
        run_spec = build_run_spec(
            config=card,
            scan_name="mapper-unit",
            task_root="/tmp",
            task_result_dir="/tmp/out",
            sampler_name="Bridson",
        )
        payload = build_payload(run_spec=run_spec, sampler_state={"index": 0})
        # Real hash present; simulate resume-time Operas outage.
        if payload["integrity"].get("operas_constants_hash") is None:
            payload["integrity"]["operas_constants_hash"] = "ab" * 32
        with mock.patch(
            "jarvishep2.Sampling.runtime_checkpoint.operas_constants_fingerprint",
            return_value=None,
        ):
            ok, reason = check_operas_constants_fingerprint(payload)
        self.assertTrue(ok, msg=reason)
        self.assertTrue(reason.startswith("skip:"), msg=reason)
        self.assertNotIn("changed", reason)


def format_issues(report) -> str:
    return "\n".join(i.format_line() for i in report.issues)


if __name__ == "__main__":
    unittest.main()
