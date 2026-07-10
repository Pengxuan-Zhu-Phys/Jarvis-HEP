#!/usr/bin/env python3
"""Sampler factory dispatch via a register table (V1-compatible method names)."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from typing import Any

from jarvishep2.Sampling.sampler import SamplingVirtial

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


def _factory_adaptive_level_set() -> SamplingVirtial:
    from jarvishep2.Sampling.adaptive_level_set import AdaptiveLevelSetSampler

    return AdaptiveLevelSetSampler()


def register_builtin_samplers(*, override: bool = True) -> None:
    """Install the built-in V2 sampler set."""
    Distributor.register(
        "Bridson", _factory_bridson, stateless=True, resume="implemented", override=override
    )
    Distributor.register(
        "Random", _factory_random, stateless=True, resume="implemented", override=override
    )
    Distributor.register(
        "Grid", _factory_grid, stateless=True, resume="implemented", override=override
    )
    Distributor.register(
        "CSV", _factory_csv, stateless=True, resume="implemented", override=override
    )
    Distributor.register(
        "AdaptiveLevelSet",
        _factory_adaptive_level_set,
        stateless=False,
        resume="implemented",
        override=override,
    )


register_builtin_samplers()


__all__ = [
    "Distributor",
    "STATELESS_METHODS",
    "SamplerRegistration",
    "register_builtin_samplers",
]
