#!/usr/bin/env python3
"""Sampler factory dispatch via a register table (V1-compatible method names)."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from typing import Any

from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.sampler_catalog import builtin_specs

SamplerFactory = Callable[[], SamplingVirtial]


@dataclass(frozen=True)
class SamplerRegistration:
    """One entry in the Distributor sampler registry."""

    name: str
    factory: SamplerFactory
    stateless: bool = True
    resume: str = "implemented"


# Process-local registry: name → registration.
_REGISTRY: dict[str, SamplerRegistration] = {}


class _StatelessMethodsView:
    """Live view of stateless method names (safe under ``from … import``)."""

    def __contains__(self, item: object) -> bool:
        registration = _REGISTRY.get(str(item).strip())
        return registration is not None and registration.stateless

    def __iter__(self) -> Iterator[str]:
        return (reg.name for reg in _REGISTRY.values() if reg.stateless)

    def __len__(self) -> int:
        return sum(1 for reg in _REGISTRY.values() if reg.stateless)

    def __repr__(self) -> str:
        return f"frozenset({sorted(self)!r})"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, (set, frozenset, _StatelessMethodsView)):
            return set(self) == set(other)
        return NotImplemented


# Backward-compatible name used by core.run(); always reflects the live registry.
STATELESS_METHODS = _StatelessMethodsView()


class Distributor:
    """Maps ``Sampling.Method`` strings to concrete sampler instances."""

    # Dynamic mirror for any code that still reads the class attribute.
    RESUME_SUPPORT_STATUS: dict[str, str] = {}

    @classmethod
    def register(
        cls,
        name: str,
        factory: SamplerFactory,
        *,
        stateless: bool = True,
        resume: str = "implemented",
        override: bool = False,
    ) -> None:
        """Register a sampler factory under a YAML method name."""
        method = str(name).strip()
        if not method:
            raise ValueError("sampler method name must not be empty")
        if not callable(factory):
            raise TypeError(f"sampler factory for '{method}' must be callable")
        if method in _REGISTRY and not override:
            raise ValueError(
                f"sampler method '{method}' is already registered; pass override=True to replace"
            )
        _REGISTRY[method] = SamplerRegistration(
            name=method,
            factory=factory,
            stateless=bool(stateless),
            resume=str(resume),
        )
        cls.RESUME_SUPPORT_STATUS = cls.resume_support_status()

    @classmethod
    def unregister(cls, name: str) -> None:
        """Remove a registration (tests / dynamic cleanup)."""
        method = str(name).strip()
        if method in _REGISTRY:
            del _REGISTRY[method]
            cls.RESUME_SUPPORT_STATUS = cls.resume_support_status()

    @classmethod
    def clear_registry_for_tests(cls) -> None:
        """Drop all registrations (tests only)."""
        _REGISTRY.clear()
        cls.RESUME_SUPPORT_STATUS = {}

    @classmethod
    def available_methods(cls) -> list[str]:
        return sorted(_REGISTRY.keys())

    @classmethod
    def get_resume_status(cls, method: str) -> str:
        registration = _REGISTRY.get(str(method).strip())
        if registration is None:
            return "intentionally unsupported"
        return registration.resume

    @classmethod
    def resume_support_status(cls) -> dict[str, str]:
        return {name: reg.resume for name, reg in _REGISTRY.items()}

    @staticmethod
    def set_method(method: str) -> SamplingVirtial:
        key = str(method).strip()
        registration = _REGISTRY.get(key)
        if registration is None:
            available = ", ".join(Distributor.available_methods()) or "none"
            raise NotImplementedError(
                f"Sampling.Method '{method}' is not implemented in Jarvis-HEP V2 yet. "
                f"Available: {available}"
            )
        sampler = registration.factory()
        if not isinstance(sampler, SamplingVirtial):
            raise TypeError(
                f"sampler factory for '{key}' must return SamplingVirtial, got {type(sampler)}"
            )
        return sampler


def _factory_bridson() -> SamplingVirtial:
    from jarvishep2.Sampling.bridson import Bridson

    return Bridson()


def _factory_random() -> SamplingVirtial:
    from jarvishep2.Sampling.randoms import RandomS

    return RandomS()


def _factory_grid() -> SamplingVirtial:
    from jarvishep2.Sampling.grid import Grid

    return Grid()


def _factory_csv() -> SamplingVirtial:
    from jarvishep2.Sampling.csv_sampler import CSVSampler

    return CSVSampler()


def _factory_adaptive_bridson() -> SamplingVirtial:
    from jarvishep2.Sampling.adaptive_bridson import AdaptiveBridsonSampler

    return AdaptiveBridsonSampler()


def _factory_mcmc() -> SamplingVirtial:
    from jarvishep2.Sampling.mcmc import create_mcmc

    return create_mcmc()


def _factory_toymcmc() -> SamplingVirtial:
    from jarvishep2.Sampling.toymcmc import create_toymcmc

    return create_toymcmc()


def _factory_ammcmc() -> SamplingVirtial:
    from jarvishep2.Sampling.ammcmc import create_ammcmc

    return create_ammcmc()


def _factory_am() -> SamplingVirtial:
    from jarvishep2.Sampling.am import create_am

    return create_am()


def _factory_dram() -> SamplingVirtial:
    from jarvishep2.Sampling.dram import create_dram

    return create_dram()


def _factory_ensemble_mcmc() -> SamplingVirtial:
    from jarvishep2.Sampling.ensemble_mcmc import create_ensemble

    return create_ensemble()


def _factory_ensemble() -> SamplingVirtial:
    from jarvishep2.Sampling.ensemble_mcmc import create_ensemble_alias

    return create_ensemble_alias()


def _factory_demcmc() -> SamplingVirtial:
    from jarvishep2.Sampling.demcmc import create_demcmc

    return create_demcmc()


def _factory_ptmcmc() -> SamplingVirtial:
    from jarvishep2.Sampling.ptmcmc import create_ptmcmc

    return create_ptmcmc()


def _factory_ptensemble() -> SamplingVirtial:
    from jarvishep2.Sampling.ptensemble import create_pt_ensemble

    return create_pt_ensemble()


def _factory_dynesty() -> SamplingVirtial:
    from jarvishep2.Sampling.dynesty_sampler import create_dynesty

    return create_dynesty()


def _factory_multinest() -> SamplingVirtial:
    from jarvishep2.Sampling.multinest_sampler import create_multinest

    return create_multinest()


def _builtin_factories() -> dict[str, SamplerFactory]:
    """YAML-name → lazy factory. Keys must match ``sampler_catalog``."""
    return {
        "Bridson": _factory_bridson,
        "Random": _factory_random,
        "CSV": _factory_csv,
        "Grid": _factory_grid,
        "AdaptiveBridson": _factory_adaptive_bridson,
        "MCMC": _factory_mcmc,
        "ToyMCMC": _factory_toymcmc,
        "AMMCMC": _factory_ammcmc,
        "AM": _factory_am,
        "DRAM": _factory_dram,
        "EnsembleMCMC": _factory_ensemble_mcmc,
        "Ensemble": _factory_ensemble,
        "DEMCMC": _factory_demcmc,
        "PTMCMC": _factory_ptmcmc,
        "PTEnsemble": _factory_ptensemble,
        "Dynesty": _factory_dynesty,
        "MultiNest": _factory_multinest,
    }


def register_builtin_samplers(*, override: bool = True) -> None:
    """Install the built-in V2 sampler set from ``sampler_catalog``."""
    factories = _builtin_factories()
    missing = [spec.name for spec in builtin_specs() if spec.name not in factories]
    extra = sorted(set(factories) - {spec.name for spec in builtin_specs()})
    if missing or extra:
        raise RuntimeError(
            "sampler factory map drifted from sampler_catalog: "
            f"missing={missing or 'none'} extra={extra or 'none'}"
        )
    for spec in builtin_specs():
        Distributor.register(
            spec.name,
            factories[spec.name],
            stateless=spec.stateless,
            resume=spec.resume,
            override=override,
        )


register_builtin_samplers()


__all__ = [
    "Distributor",
    "STATELESS_METHODS",
    "SamplerRegistration",
    "register_builtin_samplers",
]
