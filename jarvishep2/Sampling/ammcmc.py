"""Adaptive Metropolis MCMC sampler."""

from jarvishep2.Sampling.adaptive_mcmc import AdaptiveMCMCBase
from jarvishep2.Sampling.Source.MCMC.ammcmc_chain import AMMCMCChain


class AMMCMCSampler(AdaptiveMCMCBase):
    method = "AMMCMC"

    def _make_engine(self, chain_id, scale, rng):  # noqa: ANN001
        initial = rng.random(self._dim)
        return AMMCMCChain(
            initial,
            scale,
            self._niters,
            adapt_enabled=self._adapt_enabled,
            adapt_start_iter=self._adapt_start_iter,
            adapt_window=self._adapt_window,
            adapt_eps=self._adapt_eps,
            adapt_scale=self._adapt_scale,
            rng=rng,
        )


def create_ammcmc() -> AMMCMCSampler:
    return AMMCMCSampler()


__all__ = ["AMMCMCSampler", "create_ammcmc"]
