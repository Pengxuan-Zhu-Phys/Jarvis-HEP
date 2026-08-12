"""AM alias sampler."""

from jarvishep2.Sampling.ammcmc import AMMCMCSampler


class AMSampler(AMMCMCSampler):
    method = "AM"


def create_am() -> AMSampler:
    return AMSampler()


__all__ = ["AMSampler", "create_am"]
