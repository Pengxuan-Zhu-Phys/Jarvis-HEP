#!/usr/bin/env python3
"""Builtin Sampling.Method catalog (names and flags only).

This module is import-safe for ``Jarvis validate`` / ``Jarvis man``: it must not
import sampler implementations, Redis, or Workers. Factories live in
``distributor.py`` and stay lazy.

D25.1: Distributor, JSON Schema manifest keys, contracts method families,
``STATISTICAL_METHODS``, and man sampler families all derive from this table.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator


_SCHEMA_PREFIX = "https://jarvis-hep.org/schema/v2/sampling/methods"


@dataclass(frozen=True)
class SamplerSpec:
    """One builtin YAML ``Sampling.Method``."""

    name: str
    stateless: bool = True
    resume: str = "implemented"
    needs_variables: bool = True
    statistical: bool = False
    nested: bool = False
    mcmc: bool = False
    mcmc_pipeline: str | None = None
    family: str = "simple"
    schema_id: str = ""


def _method(
    name: str,
    *,
    schema: str,
    stateless: bool = True,
    needs_variables: bool = True,
    statistical: bool = False,
    nested: bool = False,
    mcmc: bool = False,
    mcmc_pipeline: str | None = None,
    family: str = "simple",
    resume: str = "implemented",
) -> SamplerSpec:
    return SamplerSpec(
        name=name,
        stateless=stateless,
        resume=resume,
        needs_variables=needs_variables,
        statistical=statistical,
        nested=nested,
        mcmc=mcmc,
        mcmc_pipeline=mcmc_pipeline,
        family=family,
        schema_id=f"{_SCHEMA_PREFIX}/{schema}",
    )


_SPECS: tuple[SamplerSpec, ...] = (
    _method("Bridson", schema="bridson.json"),
    _method("Random", schema="random.json"),
    _method("CSV", schema="csv.json", needs_variables=False),
    _method("Grid", schema="grid.json"),
    _method(
        "AdaptiveBridson",
        schema="adaptive_bridson.json",
        stateless=False,
        family="adaptive",
    ),
    _method(
        "MCMC",
        schema="mcmc.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="async_independent",
        family="MCMC",
    ),
    _method(
        "ToyMCMC",
        schema="toymcmc.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="async_independent",
        family="MCMC",
    ),
    _method(
        "AMMCMC",
        schema="ammcmc.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="async_independent",
        family="MCMC",
    ),
    _method(
        "AM",
        schema="am.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="async_independent",
        family="MCMC",
    ),
    _method(
        "DRAM",
        schema="dram.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="async_independent",
        family="MCMC",
    ),
    _method(
        "EnsembleMCMC",
        schema="ensemble_mcmc.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="barrier_coupled",
        family="MCMC",
    ),
    _method(
        "Ensemble",
        schema="ensemble_base.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="barrier_coupled",
        family="MCMC",
    ),
    _method(
        "DEMCMC",
        schema="demcmc.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="barrier_coupled",
        family="MCMC",
    ),
    _method(
        "PTMCMC",
        schema="ptmcmc.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="barrier_coupled",
        family="MCMC",
    ),
    _method(
        "PTEnsemble",
        schema="ptensemble.json",
        stateless=False,
        statistical=True,
        mcmc=True,
        mcmc_pipeline="barrier_coupled",
        family="MCMC",
    ),
    _method(
        "Dynesty",
        schema="dynesty.json",
        stateless=False,
        statistical=True,
        nested=True,
        family="nested",
    ),
    _method(
        "MultiNest",
        schema="multinest.json",
        stateless=False,
        statistical=True,
        nested=True,
        family="nested",
    ),
)

_BY_NAME: dict[str, SamplerSpec] = {spec.name: spec for spec in _SPECS}


def builtin_specs() -> tuple[SamplerSpec, ...]:
    return _SPECS


def spec_for(name: str) -> SamplerSpec | None:
    return _BY_NAME.get(str(name).strip())


def is_known(name: str) -> bool:
    return spec_for(name) is not None


def names() -> frozenset[str]:
    return frozenset(_BY_NAME)


def available_methods() -> list[str]:
    """Sorted YAML method names (validation / error messages)."""
    return sorted(_BY_NAME)


def schema_ids() -> dict[str, str]:
    return {spec.name: spec.schema_id for spec in _SPECS}


def is_stateless(name: str) -> bool:
    spec = spec_for(name)
    return bool(spec is not None and spec.stateless)


def methods_need_variables() -> frozenset[str]:
    return frozenset(spec.name for spec in _SPECS if spec.needs_variables)


def nested_names() -> frozenset[str]:
    return frozenset(spec.name for spec in _SPECS if spec.nested)


def mcmc_names() -> frozenset[str]:
    return frozenset(spec.name for spec in _SPECS if spec.mcmc)


def statistical_names() -> frozenset[str]:
    return frozenset(spec.name for spec in _SPECS if spec.statistical)


def mcmc_independent_names() -> frozenset[str]:
    return frozenset(
        spec.name for spec in _SPECS if spec.mcmc_pipeline == "async_independent"
    )


def mcmc_coupled_names() -> frozenset[str]:
    return frozenset(
        spec.name for spec in _SPECS if spec.mcmc_pipeline == "barrier_coupled"
    )


def family_of(name: str) -> str:
    spec = spec_for(name)
    return spec.family if spec is not None else "simple"


def iter_specs() -> Iterator[SamplerSpec]:
    return iter(_SPECS)


__all__ = [
    "SamplerSpec",
    "available_methods",
    "builtin_specs",
    "family_of",
    "is_known",
    "is_stateless",
    "iter_specs",
    "mcmc_coupled_names",
    "mcmc_independent_names",
    "mcmc_names",
    "methods_need_variables",
    "names",
    "nested_names",
    "schema_ids",
    "spec_for",
    "statistical_names",
]
