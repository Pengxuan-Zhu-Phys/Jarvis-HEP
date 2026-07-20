#!/usr/bin/env python3
"""D13.8 flat feedback return contract tests."""

from __future__ import annotations

import math
import unittest

import numpy as np

from jarvishep2.feedback_return import (
    WIRE_LOGL_KEY,
    build_feedback_record,
    decode_feedback_float,
    encode_feedback_float,
    feedback_logl,
    is_unusable_logl,
    normalize_feedback_record,
    resolve_feedback_return,
)
from jarvishep2.likelihood import LogLikelihoodEvaluator
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.sample import Sample


class CodecTests(unittest.TestCase):
    def test_encode_decode_neginf(self) -> None:
        enc = encode_feedback_float(float("-inf"))
        self.assertEqual(enc, "-inf")
        self.assertTrue(math.isinf(decode_feedback_float(enc)))
        self.assertLess(decode_feedback_float(enc), 0)


class BuildRecordTests(unittest.TestCase):
    def test_minimal_flat_uuid_logl(self) -> None:
        sample = Sample(uuid="abc", u_coords=np.array([0.1]))
        sample.observables = {"LogL": -3.5, "z": 99.0, "junk": [1, 2]}
        rec = build_feedback_record(
            sample, spec={"mode": "minimal", "include_logl": True, "fields": []}
        )
        self.assertEqual(set(rec.keys()), {"uuid", WIRE_LOGL_KEY})
        self.assertEqual(rec["uuid"], "abc")
        self.assertEqual(rec[WIRE_LOGL_KEY], -3.5)
        self.assertNotIn("status", rec)
        self.assertNotIn("observables", rec)

    def test_missing_logl_becomes_neginf(self) -> None:
        sample = Sample(uuid="x", u_coords=np.array([]))
        sample.observables = {"z": 1.0}
        rec = build_feedback_record(
            sample, spec={"mode": "minimal", "include_logl": True}
        )
        self.assertEqual(rec[WIRE_LOGL_KEY], "-inf")
        self.assertTrue(is_unusable_logl(feedback_logl(rec)))

    def test_fields_are_flat_siblings(self) -> None:
        sample = Sample(uuid="f", u_coords=np.array([0.0]))
        sample.observables = {"LogL": -1.0, "delta_chi2": 3.84, "m0": 100.0}
        rec = build_feedback_record(
            sample,
            spec={
                "mode": "fields",
                "include_logl": True,
                "fields": ["delta_chi2", "m0"],
            },
        )
        self.assertEqual(rec["uuid"], "f")
        self.assertEqual(rec[WIRE_LOGL_KEY], -1.0)
        self.assertEqual(rec["delta_chi2"], 3.84)
        self.assertEqual(rec["m0"], 100.0)
        self.assertNotIn("observables", rec)


class ResolveSpecTests(unittest.TestCase):
    def test_default_minimal(self) -> None:
        spec = resolve_feedback_return({"Sampling": {"Method": "Dynesty"}})
        self.assertEqual(spec["mode"], "minimal")
        self.assertTrue(spec["include_logl"])

    def test_als_fields_from_target_expression(self) -> None:
        spec = resolve_feedback_return(
            {
                "Sampling": {
                    "Method": "AdaptiveLevelSet",
                    "AdaptiveLevelSet": {"target_expression": "delta_chi2 + LogL"},
                }
            }
        )
        self.assertEqual(spec["mode"], "fields")
        self.assertIn("delta_chi2", spec["fields"])
        self.assertIn("LogL", spec["fields"])

    def test_yaml_override(self) -> None:
        spec = resolve_feedback_return(
            {
                "Sampling": {
                    "Method": "MCMC",
                    "FeedbackReturn": {"mode": "all_flat"},
                }
            }
        )
        self.assertEqual(spec["mode"], "all_flat")


class RedisRoundTripTests(unittest.TestCase):
    def test_flat_publish_pull(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_feedback({"uuid": "u1", "logL": -2.5})
        rec = queue.pull_feedback(timeout=1)
        assert rec is not None
        self.assertEqual(rec["uuid"], "u1")
        self.assertAlmostEqual(float(rec["logL"]), -2.5)
        self.assertNotIn("status", rec)
        self.assertNotIn("observables", rec)

    def test_neginf_roundtrip(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_feedback({"uuid": "bad", "logL": float("-inf")})
        rec = queue.pull_feedback(timeout=1)
        assert rec is not None
        self.assertTrue(math.isinf(float(rec["logL"])))
        self.assertLess(float(rec["logL"]), 0)

    def test_legacy_nested_still_accepted(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_feedback(
            {
                "uuid": "legacy",
                "status": "Completed",
                "observables": {"LogL": -7.0, "z": 1},
            }
        )
        rec = queue.pull_feedback(timeout=1)
        assert rec is not None
        self.assertAlmostEqual(float(rec["logL"]), -7.0)
        self.assertNotIn("status", rec)


class LikelihoodNegInfTests(unittest.TestCase):
    def test_empty_expressions_logl_neginf(self) -> None:
        ev = LogLikelihoodEvaluator([])
        out = ev.evaluate({"x": 1.0})
        self.assertTrue(math.isinf(out["LogL"]))
        self.assertLess(out["LogL"], 0)

    def test_nonfinite_term_becomes_neginf(self) -> None:
        # log(0) style: expression that yields -inf stays -inf
        ev = LogLikelihoodEvaluator(
            [{"name": "LogL", "expression": "log(x)"}]
        )
        out = ev.evaluate({"x": 0.0})
        self.assertTrue(math.isinf(out["LogL"]) or out["LogL"] < -1e100)

    def test_calculate_writes_neginf_on_missing_symbol(self) -> None:
        ev = LogLikelihoodEvaluator(
            [{"name": "LogL", "expression": "missing_var + 1"}]
        )
        info: dict = {"observables": {"x": 1.0}, "uuid": "t"}
        with self.assertRaises(KeyError):
            ev.calculate(info)
        self.assertTrue(math.isinf(float(info["observables"]["LogL"])))
        self.assertLess(float(info["observables"]["LogL"]), 0)


class NormalizeLegacyTests(unittest.TestCase):
    def test_normalize_nested(self) -> None:
        rec = normalize_feedback_record(
            {"uuid": "a", "status": "Failed", "observables": {"LogL": -1.0}}
        )
        self.assertEqual(rec["uuid"], "a")
        self.assertAlmostEqual(float(rec["logL"]), -1.0)
        self.assertNotIn("status", rec)


if __name__ == "__main__":
    unittest.main()
