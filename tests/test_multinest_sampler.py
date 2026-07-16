#!/usr/bin/env python3
"""MultiNest (static NestedSampler) + plot CSV path tests."""

from __future__ import annotations

import csv
import os
import tempfile
import unittest

import numpy as np

from jarvishep2.Sampling.dynesty_sampler import export_dynesty_results_csv
from jarvishep2.Sampling.multinest_sampler import MultiNestSampler, create_multinest
from jarvishep2.distributor import STATELESS_METHODS, Distributor


class MultiNestRegistrationTests(unittest.TestCase):
    def test_registered_as_feedback(self) -> None:
        self.assertNotIn("MultiNest", STATELESS_METHODS)
        sampler = Distributor.set_method("MultiNest")
        self.assertIsInstance(sampler, MultiNestSampler)
        self.assertEqual(sampler.method, "MultiNest")

    def test_create_factory(self) -> None:
        self.assertIsInstance(create_multinest(), MultiNestSampler)


class MultiNestConfigTests(unittest.TestCase):
    def test_always_static_even_if_dynamic_true(self) -> None:
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
                        "dynamic": True,  # Dynesty-style key must be ignored
                        "rseed": 3,
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
                rows = list(reader)
            self.assertGreater(len(rows), 0)
            self.assertIn("log_Like", reader.fieldnames or [])
            self.assertIn("log_Evidence", reader.fieldnames or [])

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


if __name__ == "__main__":
    unittest.main()
