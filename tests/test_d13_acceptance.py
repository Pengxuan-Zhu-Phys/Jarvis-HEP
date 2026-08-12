#!/usr/bin/env python3
"""D13.6 acceptance gates: registry, diagnostics contract, nested Dynamic smoke."""

from __future__ import annotations

import csv
import json
import os
import tempfile
import threading
import unittest

import numpy as np

from jarvishep2.Sampling.Source.Dynesty.py.dynesty import DynamicNestedSampler
from jarvishep2.Sampling.Source.Dynesty.py.dynesty.jarvis_uuid import looks_uuid_augmented
from jarvishep2.Sampling.diagnostics_export import (
    DATABASE_CHAIN_DIAGNOSTIC_COLUMNS,
    DATABASE_NESTED_RESULT_COLUMNS,
    export_mcmc_diagnostics,
    export_nested_diagnostics,
)
from jarvishep2.Sampling.dynesty_sampler import (
    DynestySampler,
    _jarvis_prior_transform,
    export_dynesty_results_csv,
)
from jarvishep2.Sampling.dram import DRAMSampler
from jarvishep2.Sampling.mcmc_sampler import MCMCBaseSampler, _effective_sample_size
from jarvishep2.Sampling.multinest_sampler import MultiNestSampler
from jarvishep2.distributor import STATELESS_METHODS, Distributor


# Full D13 method surface (aliases included).
D13_METHODS = (
    "MCMC",
    "ToyMCMC",
    "AMMCMC",
    "AM",
    "DRAM",
    "EnsembleMCMC",
    "Ensemble",
    "DEMCMC",
    "PTMCMC",
    "PT",
    "PTEnsemble",
    "Dynesty",
    "MultiNest",
)


class RegistryAcceptanceTests(unittest.TestCase):
    def test_all_d13_methods_registered(self) -> None:
        for name in D13_METHODS:
            with self.subTest(method=name):
                self.assertNotIn(name, STATELESS_METHODS)
                sampler = Distributor.set_method(name)
                self.assertIsNotNone(sampler)
                self.assertEqual(getattr(sampler, "method", name), name)

    def test_diagnostic_column_contract_constants(self) -> None:
        for col in ("chain_id", "step", "accepted", "weight"):
            self.assertIn(col, DATABASE_CHAIN_DIAGNOSTIC_COLUMNS)
        for col in ("log_weight", "log_Like", "log_Evidence"):
            self.assertIn(col, DATABASE_NESTED_RESULT_COLUMNS)


class EssDiagnosticTests(unittest.TestCase):
    def test_ess_independent_noise_near_n(self) -> None:
        """D13.7d: uncorrelated series → ESS close to usable length (N/2)."""
        rng = np.random.default_rng(0)
        series = rng.normal(size=200).tolist()
        ess = _effective_sample_size(series)
        self.assertIsNotNone(ess)
        assert ess is not None
        # Second-half burn-in → ~100 samples; uncorrelated → ESS ~ n.
        self.assertGreater(ess, 40.0)
        self.assertLessEqual(ess, 100.0 + 1e-6)

    def test_ess_highly_correlated_is_small(self) -> None:
        # Positively autocorrelated AR(1)-like series → ESS << N.
        x = [0.0]
        for _ in range(199):
            x.append(0.95 * x[-1] + 0.05)
        ess = _effective_sample_size(x)
        self.assertIsNotNone(ess)
        assert ess is not None
        self.assertLess(ess, 30.0)


class McmcDiagnosticsExportTests(unittest.TestCase):
    def test_export_writes_summary_and_history(self) -> None:
        class _Ev:
            def __init__(self, it, acc, logl):
                self.iter = it
                self.accepted = acc
                self.logl = logl
                self.temperature = 1.0
                self.state = "UPDATE"

        class _Hist:
            def all(self):
                return [_Ev(0, True, -1.0), _Ev(1, False, -2.0)]

        class _Chain:
            chain_id = 0
            history = _Hist()

        class _Reg:
            def all(self):
                return [_Chain()]

        with tempfile.TemporaryDirectory() as tmp:
            written = export_mcmc_diagnostics(
                tmp,
                summary={"method": "MCMC", "accept_rate": 0.5, "rhat_logl": 1.01},
                registry=_Reg(),
            )
            self.assertIn("sampler_summary", written)
            self.assertIn("chain_history", written)
            with open(written["sampler_summary"], encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(payload["method"], "MCMC")
            with open(written["chain_history"], encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["accepted"], "1")
            self.assertEqual(rows[0]["weight"], "1")
            self.assertEqual(rows[1]["accepted"], "0")
            self.assertEqual(rows[1]["weight"], "0")
            for col in ("chain_id", "step", "accepted", "weight", "logl"):
                self.assertIn(col, rows[0])


class NestedDiagnosticsExportTests(unittest.TestCase):
    def test_nested_summary_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            written = export_nested_diagnostics(
                tmp, summary={"method": "Dynesty", "logz": -1.2, "logzerr": 0.1}
            )
            self.assertIn("sampler_summary", written)
            with open(written["sampler_summary"], encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(payload["logz"], -1.2)

    def test_nested_summary_strips_samples_uid_list(self) -> None:
        """On-disk summary keeps count only — full uids live in dynesty_result.csv."""
        with tempfile.TemporaryDirectory() as tmp:
            written = export_nested_diagnostics(
                tmp,
                summary={
                    "method": "Dynesty",
                    "logz": -2.0,
                    "samples_uid": [f"u{i}" for i in range(50)],
                },
            )
            with open(written["sampler_summary"], encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertNotIn("samples_uid", payload)
            self.assertEqual(payload["samples_uid_count"], 50)


class NestedCsvSchemaAcceptanceTests(unittest.TestCase):
    def test_export_schema_has_nested_contract_columns(self) -> None:
        n = 4
        results = {
            "logl": np.linspace(-5, -1, n),
            "logwt": np.linspace(-2, -0.1, n),
            "logvol": np.linspace(0, -1, n),
            "logz": np.linspace(-4, -1.5, n),
            "logzerr": np.full(n, 0.05),
            "samples_n": np.full(n, 20),
            "ncall": np.arange(1, n + 1),
            "samples_it": np.arange(n),
            "samples_id": np.arange(n),
            "information": np.zeros(n),
            "samples": np.zeros((n, 2)),
            "samples_u": np.zeros((n, 2)),
            "samples_uid": [f"u{i}" for i in range(n)],
            "nlive": 20,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = export_dynesty_results_csv(
                results, os.path.join(tmp, "dynesty_result.csv")
            )
            with open(path, encoding="utf-8", newline="") as handle:
                fieldnames = list(csv.DictReader(handle).fieldnames or [])
        for col in DATABASE_NESTED_RESULT_COLUMNS:
            self.assertIn(col, fieldnames)


class DynamicNestedUuidAcceptanceTests(unittest.TestCase):
    def test_dynamic_nested_runs_with_uuid_prior(self) -> None:
        """D13.6: DynamicNestedSampler + UUID channel must complete (not TypeError)."""

        def loglike(x):
            arr = np.asarray(x, dtype=object).reshape(-1)
            if looks_uuid_augmented(arr):
                u = np.asarray(arr[:-1], dtype=float)
            else:
                u = np.asarray(arr, dtype=float)
            return float(-0.5 * np.sum((u - 0.5) ** 2) / 0.1**2)

        sampler = DynamicNestedSampler(
            loglikelihood=loglike,
            prior_transform=_jarvis_prior_transform,
            ndim=2,
            nlive=15,
            bound="single",
            sample="unif",
            rstate=np.random.default_rng(0),
        )
        sampler.run_nested(maxiter=25, dlogz_init=1.0, print_progress=False)
        res = sampler.results
        self.assertIn("logz", res.keys())
        self.assertGreater(len(res["logz"]), 0)
        if "samples_uid" in res.keys():
            uids = [str(u) for u in res["samples_uid"]]
            self.assertTrue(any(uids))

    def test_dynamic_nested_batch_strips_uuid(self) -> None:
        """Dynamic add_batch must strip UUID from live_v (shape ndim, not ndim+1)."""

        def loglike(x):
            arr = np.asarray(x, dtype=object).reshape(-1)
            if looks_uuid_augmented(arr):
                u = np.asarray(arr[:-1], dtype=float)
            else:
                u = np.asarray(arr, dtype=float)
            return float(-0.5 * np.sum((u - 0.5) ** 2) / 0.1**2)

        sampler = DynamicNestedSampler(
            loglikelihood=loglike,
            prior_transform=_jarvis_prior_transform,
            ndim=2,
            nlive=20,
            bound="single",
            sample="unif",
            rstate=np.random.default_rng(1),
        )
        # Force base convergence then at least one dynamic batch (hits the
        # live_v[i] = v assignment that previously broadcast shape (3,)→(2,)).
        sampler.run_nested(
            dlogz_init=0.5,
            maxcall_init=800,
            nlive_batch=15,
            maxbatch=2,
            print_progress=False,
        )
        res = sampler.results
        self.assertIn("logz", res.keys())
        samples = np.asarray(res["samples"])
        self.assertEqual(samples.shape[1], 2)
        # UUID channel must land in results for CSV export.
        self.assertIn("samples_uid", res.keys())
        uids = [str(u) for u in np.asarray(res["samples_uid"]).reshape(-1)]
        nonempty = [u for u in uids if u.strip()]
        self.assertGreater(len(nonempty), 0, msg="expected non-empty samples_uid")


class SamplerConfigAcceptanceTests(unittest.TestCase):
    def test_mcmc_bounds_surface(self) -> None:
        s = DRAMSampler()
        self.assertIsInstance(s, MCMCBaseSampler)
        s.set_config(
            {
                "Sampling": {
                    "Method": "DRAM",
                    "Bounds": {"num_chains": 2,
                        "num_iters": 5,
                        "dr.steps": 2,
                        "adapt.enabled": True, "seed": 1},
                    "Variables": [
                        {
                            "name": "x",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1},
                            },
                        }
                    ],
                }
            }
        )
        self.assertEqual(s._nchains, 2)
        self.assertEqual(s._niters, 5)
        self.assertEqual(s._dr_steps, 2)

    def test_dynesty_and_multinest_defaults(self) -> None:
        d = DynestySampler()
        d.set_config(
            {
                "Sampling": {
                    "Method": "Dynesty",
                    "Variables": [
                        {
                            "name": "x",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1},
                            },
                        }
                    ],
                    "Bounds": {"nlive": 50, "dlogz": 0.2},
                }
            }
        )
        self.assertTrue(d._use_dynamic)
        self.assertIn("dlogz_init", d._run_nested_kwargs)

        m = MultiNestSampler()
        m.set_config(
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
                        }
                    ],
                    "Bounds": {"nlive": 50, "dlogz": 0.2},
                }
            }
        )
        self.assertFalse(m._use_dynamic)
        self.assertIn("dlogz", m._run_nested_kwargs)


if __name__ == "__main__":
    unittest.main()
