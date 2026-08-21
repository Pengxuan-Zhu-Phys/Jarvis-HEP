"""Jarvis-HEP V2 sampler package.

YAML ``Sampling.Method`` names are constructed by
``jarvishep2.distributor.Distributor.set_method`` (see
``jarvishep2.sampler_catalog``). This module exports sampler *bases* and keeps
the historic coverage-sampler names importable.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.feedback_sampler import FeedbackSampler
from jarvishep2.Sampling.sampler import SamplingVirtial

_LAZY_ATTRS: dict[str, tuple[str, str]] = {
    "Bridson": ("jarvishep2.Sampling.bridson", "Bridson"),
    "CSVSampler": ("jarvishep2.Sampling.csv_sampler", "CSVSampler"),
    "Grid": ("jarvishep2.Sampling.grid", "Grid"),
    "MCMCBaseSampler": ("jarvishep2.Sampling.mcmc_sampler", "MCMCBaseSampler"),
    "RandomS": ("jarvishep2.Sampling.randoms", "RandomS"),
}

__all__ = [
    "Bridson",
    "CSVSampler",
    "CheckpointedSampler",
    "FeedbackSampler",
    "Grid",
    "MCMCBaseSampler",
    "RandomS",
    "SamplingVirtial",
]


def __getattr__(name: str) -> Any:
    target = _LAZY_ATTRS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr = target
    value = getattr(import_module(module_name), attr)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(__all__) | set(globals()))
