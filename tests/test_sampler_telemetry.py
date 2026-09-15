from __future__ import annotations

import json
import unittest
from types import SimpleNamespace

from jarvishep2.sampling.telemetry import (
    SAMPLER_STATUS_MAX_BYTES,
    sampler_status_json,
    sampler_status_snapshot,
)


class _Registry:
    def __init__(self, iterations: list[int]) -> None:
        self._chains = [
            SimpleNamespace(engine=SimpleNamespace(iterations=value))
            for value in iterations
        ]

    def all(self):
        return list(self._chains)


class SamplerTelemetryTests(unittest.TestCase):
    def test_random_reports_candidate_cursor_not_completed_samples(self) -> None:
        sampler = SimpleNamespace(
            method="Random",
            _index=23,
            _maxp=100,
            _accepted_index=19,
            _seed=7,
        )
        status = sampler_status_snapshot(sampler)
        self.assertEqual(status["progress"]["current"], 23)
        self.assertEqual(status["progress"]["target"], 100)
        self.assertEqual(status["metrics"]["accepted"], 19)

    def test_mcmc_progress_is_the_slowest_chain_floor(self) -> None:
        sampler = SimpleNamespace(
            method="MCMC",
            _registry=_Registry([24, 19, 22, 21]),
            _nchains=4,
            _niters=500,
            _total_accepted=20,
            _total_proposed=80,
            _finished=False,
        )
        status = sampler_status_snapshot(sampler)
        self.assertEqual(status["progress"]["current"], 19)
        self.assertEqual(status["progress"]["target"], 500)
        self.assertEqual(status["metrics"]["accept_rate"], 0.25)

    def test_nested_status_does_not_require_results_history(self) -> None:
        class Native:
            it = 42
            ncall = 314

            @property
            def results(self):
                raise AssertionError("heartbeat must not materialize results")

        sampler = SimpleNamespace(
            method="Dynesty",
            _sampler=Native(),
            _nlive=100,
            _dlogz=0.5,
            _finished=False,
        )
        status = sampler_status_snapshot(sampler)
        self.assertEqual(status["metrics"]["niter"], 42)
        self.assertEqual(status["metrics"]["ncall"], 314)

    def test_json_is_bounded_and_scalar_only(self) -> None:
        encoded = sampler_status_json(
            SimpleNamespace(method="CustomSampler", _finished=False)
        )
        self.assertLessEqual(len(encoded.encode("utf-8")), SAMPLER_STATUS_MAX_BYTES)
        self.assertEqual(json.loads(encoded)["schema"], 1)


if __name__ == "__main__":
    unittest.main()
