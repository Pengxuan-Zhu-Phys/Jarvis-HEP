#!/usr/bin/env python3
"""D18.7 coverage for the five V1 distribution types missing from V2."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.slow

import unittest
import os
import tempfile

import numpy as np
import yaml
from scipy import special, stats

from jarvishep2.Sampling.variables import Variable, load_variables
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.factory import TaskFactory
from jarvishep2.mapper import DistributionUMapper
from jarvishep2.task_validation import validate_task_config
from test_worker_mvp import _start_tcp_fakeredis


_V1_DISTRIBUTIONS = [
    ("Binomial", {"n": 12, "p": 0.3}),
    ("Poisson", {"lambda": 4.0}),
    ("Beta", {"alpha": 2.0, "beta": 5.0}),
    ("Exponential", {"rate": 2.0}),
    ("Gamma", {"shape": 2.0, "scale": 3.0}),
]


def _random_config(variables: list[dict]) -> dict:
    return {
        "Scan": {"name": "v1-distributions"},
        "Sampling": {"Method": "Random", "Bounds": {"point_number": 8}, "Variables": variables},
        "EnvReqs": {"V2": {"workers": 1}},
    }


class V1DistributionMappingTests(unittest.TestCase):
    def test_v1_distribution_cards_validate_and_load(self) -> None:
        variables = [
            {
                "name": distribution.lower(),
                "distribution": {"type": distribution, "parameters": parameters},
            }
            for distribution, parameters in _V1_DISTRIBUTIONS
        ]
        self.assertTrue(validate_task_config(_random_config(variables)).ok)
        loaded = load_variables(_random_config(variables))
        self.assertEqual(
            [item.distribution for item in loaded],
            [item[0] for item in _V1_DISTRIBUTIONS],
        )

    def test_inverse_cdf_mapping_has_the_expected_distribution_semantics(self) -> None:
        u = 0.37
        binomial = Variable("x", "", "Binomial", {"n": 12, "p": 0.3})
        binomial_value = binomial.map_standard_random_to_distribution(u)
        self.assertIsInstance(binomial_value, int)
        self.assertGreaterEqual(stats.binom.cdf(binomial_value, 12, 0.3), u)
        self.assertLess(stats.binom.cdf(binomial_value - 1, 12, 0.3), u)

        poisson = Variable("x", "", "Poisson", {"lambda": 4.0})
        poisson_value = poisson.map_standard_random_to_distribution(u)
        self.assertIsInstance(poisson_value, int)
        self.assertGreaterEqual(stats.poisson.cdf(poisson_value, 4.0), u)
        self.assertLess(stats.poisson.cdf(poisson_value - 1, 4.0), u)

        beta = Variable("x", "", "Beta", {"alpha": 2.0, "beta": 5.0})
        self.assertAlmostEqual(
            special.betainc(2.0, 5.0, beta.map_standard_random_to_distribution(u)),
            u,
            places=12,
        )

        exponential = Variable("x", "", "Exponential", {"rate": 2.0})
        self.assertAlmostEqual(
            exponential.map_standard_random_to_distribution(u), -np.log1p(-u) / 2.0
        )

        gamma = Variable("x", "", "Gamma", {"shape": 2.0, "scale": 3.0})
        gamma_value = gamma.map_standard_random_to_distribution(u)
        self.assertAlmostEqual(special.gammainc(2.0, gamma_value / 3.0), u, places=12)

    def test_worker_mapper_maps_all_five_types(self) -> None:
        variables = [
            {
                "name": distribution.lower(),
                "distribution": {"type": distribution, "parameters": parameters},
            }
            for distribution, parameters in _V1_DISTRIBUTIONS
        ]
        mapped = DistributionUMapper(variables).map(np.full(5, 0.37))
        self.assertEqual(
            set(mapped), {distribution.lower() for distribution, _ in _V1_DISTRIBUTIONS}
        )
        self.assertIsInstance(mapped["binomial"], int)
        self.assertIsInstance(mapped["poisson"], int)
        self.assertGreaterEqual(mapped["beta"], 0.0)
        self.assertLessEqual(mapped["beta"], 1.0)
        self.assertGreaterEqual(mapped["exponential"], 0.0)
        self.assertGreaterEqual(mapped["gamma"], 0.0)

    def test_parameter_domains_are_rejected_before_a_worker_starts(self) -> None:
        variables = [
            {
                "name": "binomial",
                "distribution": {"type": "Binomial", "parameters": {"n": 1.5, "p": 1.2}},
            },
            {
                "name": "poisson",
                "distribution": {"type": "Poisson", "parameters": {"lambda": -1}},
            },
            {
                "name": "beta",
                "distribution": {"type": "Beta", "parameters": {"alpha": 0, "beta": 1}},
            },
            {
                "name": "exponential",
                "distribution": {"type": "Exponential", "parameters": {"rate": 0}},
            },
            {
                "name": "gamma",
                "distribution": {"type": "Gamma", "parameters": {"shape": 1, "scale": -1}},
            },
        ]
        report = validate_task_config(_random_config(variables))
        self.assertFalse(report.ok)
        codes = {entry.code for entry in report.errors()}
        self.assertTrue(
            {"JV2-VAR-038", "JV2-VAR-039", "JV2-VAR-040", "JV2-VAR-041"}.issubset(codes)
        )

    def test_grid_accepts_its_num_parameter_for_a_v1_distribution(self) -> None:
        config = _random_config(
            [
                {
                    "name": "beta",
                    "distribution": {
                        "type": "Beta",
                        "parameters": {"alpha": 2, "beta": 5, "num": 3},
                    },
                }
            ]
        )
        config["Sampling"] = {
            "Method": "Grid",
            "Variables": config["Sampling"]["Variables"],
        }
        self.assertTrue(validate_task_config(config).ok)

        config["Sampling"]["Variables"][0]["distribution"]["parameters"]["num"] = 1.5
        report = validate_task_config(config)
        self.assertTrue(any(issue.code == "JV2-VAR-036" for issue in report.errors()))

    def test_numeric_yaml_strings_map_like_native_numbers(self) -> None:
        variable = Variable("poisson", "", "Poisson", {"lambda": "4.0"})
        self.assertEqual(
            variable.map_standard_random_to_distribution(0.37),
            Variable("poisson", "", "Poisson", {"lambda": 4.0}).map_standard_random_to_distribution(0.37),
        )


class V1DistributionWorkerIntegrationTests(unittest.TestCase):
    def tearDown(self) -> None:
        factory = TaskFactory._instance
        if factory is not None:
            factory.shutdown(wait=False)
        TaskFactory.reset_instance()

    def test_random_worker_archives_all_five_v1_distribution_types(self) -> None:
        variables = [
            {
                "name": distribution.lower(),
                "distribution": {"type": distribution, "parameters": parameters},
            }
            for distribution, parameters in _V1_DISTRIBUTIONS
        ]
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as root:
                task_yaml = os.path.join(root, "v1-distributions.yaml")
                with open(task_yaml, "w", encoding="utf-8") as handle:
                    yaml.safe_dump(_random_config(variables), handle, sort_keys=False)

                core = Jarvis2Core()
                core.load_task_yaml(task_yaml)
                core.config["task_root"] = root
                core.config["project_root"] = root
                core.config["task_result_dir"] = root
                core.config["Runtime"]["redis"] = redis_config
                core.config["Runtime"]["Watchdog"] = {"enabled": False}
                core.runtime = core.config["Runtime"]
                core._populate_info_from_config()

                outcome = core.run(write_run_summary=False)
                self.assertEqual(outcome, 8)
                records = SimpleHDF5Writer(
                    os.path.join(root, "DATABASE", "samples.hdf5")
                ).read_records()
                self.assertEqual(len(records), 8)
                for row in records:
                    self.assertIsInstance(row["binomial"], int)
                    self.assertIsInstance(row["poisson"], int)
                    self.assertGreaterEqual(row["beta"], 0.0)
                    self.assertLessEqual(row["beta"], 1.0)
                    self.assertGreaterEqual(row["exponential"], 0.0)
                    self.assertGreaterEqual(row["gamma"], 0.0)
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()
