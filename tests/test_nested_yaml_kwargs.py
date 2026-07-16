#!/usr/bin/env python3
"""Full official dynesty Bounds → constructor / run_nested wiring tests."""

from __future__ import annotations

import unittest

from jarvishep2.Sampling.dynesty_sampler import (
    DYNAMIC_RUN_NESTED_KEYS,
    NESTED_CONSTRUCTOR_USER_KEYS,
    STATIC_RUN_NESTED_KEYS,
    DynestySampler,
    extract_nested_constructor_kwargs,
    extract_run_nested_kwargs,
)
from jarvishep2.Sampling.multinest_sampler import MultiNestSampler


class ExtractConstructorKwargsTests(unittest.TestCase):
    def test_flat_and_nested_sampler_block(self) -> None:
        bounds = {
            "nlive": 200,
            "rseed": 1,
            "bound": "single",  # flat official key
            "walks": 10,
            "sampler": {
                "sample": "rwalk",
                "bootstrap": 5,
                "facc": 0.4,
                "loglikelihood": "must-strip",  # HEP-injected, banned
            },
            "run_nested": {"maxcall": 100},
            "unknown_flat": 99,
        }
        ctor = extract_nested_constructor_kwargs(bounds)
        self.assertEqual(ctor.get("bound"), "single")
        self.assertEqual(ctor.get("sample"), "rwalk")
        self.assertEqual(ctor.get("bootstrap"), 5)
        self.assertEqual(ctor.get("walks"), 10)
        self.assertEqual(ctor.get("facc"), 0.4)
        self.assertNotIn("loglikelihood", ctor)
        self.assertNotIn("unknown_flat", ctor)
        self.assertNotIn("nlive", ctor)  # nlive is meta; set by sampler

    def test_constructor_block_alias(self) -> None:
        bounds = {"constructor": {"bound": "balls", "slices": 3}}
        ctor = extract_nested_constructor_kwargs(bounds)
        self.assertEqual(ctor["bound"], "balls")
        self.assertEqual(ctor["slices"], 3)

    def test_all_official_constructor_keys_accepted(self) -> None:
        raw = {k: i for i, k in enumerate(sorted(NESTED_CONSTRUCTOR_USER_KEYS))}
        # first_update must be a dict in real use; int is fine for filter test
        ctor = extract_nested_constructor_kwargs({"sampler": raw})
        self.assertEqual(set(ctor.keys()), NESTED_CONSTRUCTOR_USER_KEYS)


class ExtractRunNestedKwargsTests(unittest.TestCase):
    def test_static_uses_dlogz(self) -> None:
        run = extract_run_nested_kwargs(
            {"dlogz": 0.2, "run_nested": {"maxcall": 10, "maxiter": 5}},
            dynamic=False,
            dlogz_default=0.5,
        )
        self.assertEqual(run["dlogz"], 0.2)
        self.assertEqual(run["maxcall"], 10)
        self.assertNotIn("dlogz_init", run)
        self.assertTrue(set(run.keys()) <= STATIC_RUN_NESTED_KEYS)

    def test_dynamic_maps_dlogz_to_dlogz_init(self) -> None:
        run = extract_run_nested_kwargs(
            {"dlogz": 0.3, "run_nested": {"maxcall": 20, "dlogz": 0.15}},
            dynamic=True,
            dlogz_default=0.5,
        )
        self.assertEqual(run["dlogz_init"], 0.15)  # run_nested.dlogz wins after map
        self.assertNotIn("dlogz", run)
        self.assertEqual(run["maxcall"], 20)
        self.assertTrue(set(run.keys()) <= DYNAMIC_RUN_NESTED_KEYS)

    def test_dynamic_accepts_dlogz_init_and_batch_keys(self) -> None:
        run = extract_run_nested_kwargs(
            {
                "run_nested": {
                    "dlogz_init": 0.01,
                    "nlive_batch": 50,
                    "maxbatch": 3,
                    "n_effective": 1000,
                    "use_stop": False,
                }
            },
            dynamic=True,
            dlogz_default=0.5,
        )
        self.assertEqual(run["dlogz_init"], 0.01)
        self.assertEqual(run["nlive_batch"], 50)
        self.assertEqual(run["maxbatch"], 3)
        self.assertEqual(run["n_effective"], 1000)
        self.assertIs(run["use_stop"], False)

    def test_unknown_run_keys_dropped(self) -> None:
        run = extract_run_nested_kwargs(
            {"run_nested": {"maxcall": 1, "not_a_dynesty_key": True}},
            dynamic=False,
            dlogz_default=0.5,
        )
        self.assertIn("maxcall", run)
        self.assertNotIn("not_a_dynesty_key", run)

    def test_static_maps_dlogz_init_alias(self) -> None:
        run = extract_run_nested_kwargs(
            {"run_nested": {"dlogz_init": 0.07}},
            dynamic=False,
            dlogz_default=0.5,
        )
        self.assertEqual(run["dlogz"], 0.07)
        self.assertNotIn("dlogz_init", run)


class SamplerSetConfigWiringTests(unittest.TestCase):
    def _vars(self):
        return [
            {
                "name": "x",
                "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}},
            }
        ]

    def test_dynesty_dynamic_default_and_full_sampler_block(self) -> None:
        s = DynestySampler()
        s.set_config(
            {
                "Sampling": {
                    "Method": "Dynesty",
                    "Variables": self._vars(),
                    "Bounds": {
                        "nlive": 80,
                        "rseed": 9,
                        "dlogz": 0.25,
                        "sampler": {
                            "bound": "multi",
                            "sample": "rwalk",
                            "walks": 40,
                            "bootstrap": 0,
                            "update_interval": 0.6,
                            "first_update": {"min_ncall": 100, "min_eff": 10.0},
                        },
                        "run_nested": {
                            "maxcall": 5000,
                            "maxiter": 200,
                            "nlive_batch": 40,
                            "print_progress": False,
                        },
                    },
                }
            }
        )
        self.assertTrue(s._use_dynamic)
        self.assertEqual(s._nlive, 80)
        self.assertEqual(s._constructor_kwargs["bound"], "multi")
        self.assertEqual(s._constructor_kwargs["sample"], "rwalk")
        self.assertEqual(s._constructor_kwargs["walks"], 40)
        self.assertEqual(s._constructor_kwargs["nlive"], 80)
        self.assertEqual(s._run_nested_kwargs["dlogz_init"], 0.25)
        self.assertNotIn("dlogz", s._run_nested_kwargs)
        self.assertEqual(s._run_nested_kwargs["nlive_batch"], 40)
        self.assertIs(s._run_nested_kwargs["print_progress"], False)

    def test_dynesty_static_mode(self) -> None:
        s = DynestySampler()
        s.set_config(
            {
                "Sampling": {
                    "Method": "Dynesty",
                    "Variables": self._vars(),
                    "Bounds": {
                        "nlive": 30,
                        "dynamic": False,
                        "dlogz": 0.4,
                        "bound": "single",
                        "run_nested": {"maxiter": 10},
                    },
                }
            }
        )
        self.assertFalse(s._use_dynamic)
        self.assertEqual(s._constructor_kwargs["bound"], "single")
        self.assertEqual(s._run_nested_kwargs["dlogz"], 0.4)
        self.assertNotIn("dlogz_init", s._run_nested_kwargs)

    def test_multinest_always_static_even_if_dynamic_true(self) -> None:
        s = MultiNestSampler()
        s.set_config(
            {
                "Sampling": {
                    "Method": "MultiNest",
                    "Variables": self._vars(),
                    "Bounds": {
                        "nlive": 40,
                        "dynamic": True,  # ignored
                        "dlogz": 0.12,
                        "sample": "unif",
                        "run_nested": {"maxcall": 1000, "add_live": False},
                    },
                }
            }
        )
        self.assertFalse(s._use_dynamic)
        self.assertEqual(s._constructor_kwargs["sample"], "unif")
        self.assertEqual(s._run_nested_kwargs["dlogz"], 0.12)
        self.assertNotIn("dlogz_init", s._run_nested_kwargs)
        self.assertIs(s._run_nested_kwargs["add_live"], False)

    def test_v1_eggbox_style_dynesty_card(self) -> None:
        """Historical EggBox card: Bounds.dlogz + run_nested.maxcall only."""
        s = DynestySampler()
        s.set_config(
            {
                "Sampling": {
                    "Method": "Dynesty",
                    "Variables": self._vars(),
                    "Bounds": {
                        "nlive": 1600,
                        "rseed": 21,
                        "dlogz": 0.1,
                        "run_nested": {
                            "print_progress": True,
                            "maxcall": 40000,
                        },
                    },
                }
            }
        )
        # Dynamic default — dlogz must become dlogz_init (would crash if left as dlogz)
        self.assertTrue(s._use_dynamic)
        self.assertEqual(s._run_nested_kwargs["dlogz_init"], 0.1)
        self.assertNotIn("dlogz", s._run_nested_kwargs)
        self.assertEqual(s._run_nested_kwargs["maxcall"], 40000)


if __name__ == "__main__":
    unittest.main()
