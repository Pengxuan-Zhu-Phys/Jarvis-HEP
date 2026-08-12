"""Baseline Metropolis MCMC sampler."""

from jarvishep2.Sampling.mcmc_sampler import MCMCBaseSampler


class MCMCSampler(MCMCBaseSampler):
    """Standard random-walk Metropolis sampler."""

    method = "MCMC"


def create_mcmc() -> MCMCSampler:
    return MCMCSampler()


__all__ = ["MCMCSampler", "create_mcmc"]
