#!/usr/bin/env python3
"""D25.6: one meaning per sampler public name."""

from __future__ import annotations

import unittest

from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.dynesty_sampler import DynestySampler
from jarvishep2.Sampling.feedback_sampler import FeedbackSampler
from jarvishep2.Sampling.mcmc import MCMCSampler
from jarvishep2.Sampling.mcmc_sampler import MCMCBaseSampler
from jarvishep2.Sampling.multinest_sampler import MultiNestSampler
from jarvishep2.distributor import Distributor
from jarvishep2.sampler_catalog import names
from jarvishep2.Sampling.sampler import SamplingVirtual, SamplingVirtial


class SamplingPublicSurfaceTests(unittest.TestCase):
    def test_yaml_mcmc_is_the_profile_not_the_base(self) -> None:
        sampler = Distributor.set_method("MCMC")
        self.assertIsInstance(sampler, MCMCSampler)
        self.assertIsInstance(sampler, MCMCBaseSampler)
        self.assertIs(type(sampler), MCMCSampler)
        self.assertEqual(type(sampler).__module__, "jarvishep2.sampling.mcmc")

    def test_mcmc_sampler_module_has_no_alias_or_factory_barrel(self) -> None:
        import jarvishep2.Sampling.mcmc_sampler as module

        self.assertFalse(hasattr(module, "create_mcmc"))
        self.assertFalse(hasattr(module, "create_toymcmc"))
        self.assertFalse(hasattr(module, "MCMCSampler"))
        self.assertIsNot(MCMCSampler, MCMCBaseSampler)

    def test_sampling_package_exports_bases_and_distributor_catalog(self) -> None:
        import jarvishep2.Sampling as sampling

        self.assertIs(sampling.SamplingVirtial, SamplingVirtial)
        self.assertIs(sampling.SamplingVirtual, SamplingVirtual)
        self.assertIs(SamplingVirtial, SamplingVirtual)
        self.assertIs(sampling.CheckpointedSampler, CheckpointedSampler)
        self.assertIs(sampling.FeedbackSampler, FeedbackSampler)
        self.assertIn("MCMC", names())
        self.assertEqual(
            Distributor.set_method("Random").method,
            "Random",
        )
        from jarvishep2.Sampling import Bridson

        self.assertEqual(Bridson.__name__, "Bridson")

    def test_nested_samplers_are_checkpointed_not_generation_barrier(self) -> None:
        self.assertTrue(issubclass(DynestySampler, CheckpointedSampler))
        self.assertFalse(issubclass(DynestySampler, FeedbackSampler))
        self.assertTrue(issubclass(MultiNestSampler, CheckpointedSampler))
        self.assertFalse(issubclass(MultiNestSampler, FeedbackSampler))
        self.assertFalse(hasattr(DynestySampler, "propose_generation"))
        self.assertFalse(hasattr(DynestySampler, "absorb_generation"))


if __name__ == "__main__":
    unittest.main()
