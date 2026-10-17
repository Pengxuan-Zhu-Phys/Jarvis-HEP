#!/usr/bin/env python3
"""Sampling.Method catalog (names and flags only).

This module is import-safe for ``Jarvis validate`` / ``Jarvis man``: it must not
import sampler implementations, Redis, or Workers. Factories live in
``distributor.py`` and stay lazy.

D25.1: Distributor, JSON Schema manifest keys, contracts method families,
``STATISTICAL_METHODS``, and man sampler families all derive from this table.

Three layers share one ``Sampling.Method`` namespace (external sampler S1,
``Jarvis-Books Docs/Jarvis-HEP/specs/sampler-integration.md`` §6.1):

* ``builtin`` — the samplers implemented in this package (``_SPECS``);
* ``builtin-plugin`` — samplers shipped with Jarvis-HEP but written in the
  public plug-in format (``_BUILTIN_PLUGINS``), run through ``PluginSampler``;
* ``user`` — plug-ins loaded with ``Jarvis sampler load``, read from the
  Jarvis-HEP-Sampler registry (``~/.jarvis/samplers/sources.json``).

Built-in names cannot be overridden. Lookups merge the layers at call time, so
callers must not freeze name sets at import.
"""

from __future__ import annotations

from collections.abc import Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterator

from jarvis_hep_sampler.errors import PluginConfigError, RegistryError

SOURCE_BUILTIN = "builtin"
SOURCE_BUILTIN_PLUGIN = "builtin-plugin"
SOURCE_USER = "user"
SOURCE_SCAN = "scan"
PLUGIN_FAMILY = "plugin"

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
    source: str = SOURCE_BUILTIN

    @property
    def is_plugin(self) -> bool:
        return self.source != SOURCE_BUILTIN


@dataclass(frozen=True)
class PluginRecord:
    """What validation, man pages and the run need to know about a plug-in method.

    For user plug-ins this comes from the registry (no plug-in code runs). For
    built-in plug-ins it comes from loading the shipped file once.
    """

    method: str
    source: str
    version: str
    profile: str
    description: str
    statistical: bool
    options_schema: Mapping[str, Any] | None
    hooks: tuple[str, ...]
    source_path: str
    snapshot_path: str
    source_sha256: str
    loaded_at: str = ""


class SamplerRegistryPreflightError(PluginConfigError):
    """A registry problem that must stop a run before runtime work begins."""


@dataclass(frozen=True)
class RegistryReportEntry:
    """One method in the merged, read-only sampler registry report."""

    method: str
    source: str
    status: str
    problems: tuple[str, ...] = ()
    repair_command: str | None = None


@dataclass(frozen=True)
class RegistryReport:
    """Built-in and external names, with an optional registry-wide read error."""

    entries: tuple[RegistryReportEntry, ...]
    error: str | None = None

    @property
    def conflicts(self) -> tuple[RegistryReportEntry, ...]:
        return tuple(entry for entry in self.entries if entry.status == "name-conflict")

    def matching(self, method: str) -> tuple[RegistryReportEntry, ...]:
        folded = str(method).strip().casefold()
        return tuple(entry for entry in self.entries if entry.method.casefold() == folded)

    def for_method(self, method: str, *, source: str | None = None) -> RegistryReportEntry | None:
        for entry in self.matching(method):
            if source is None or entry.source == source:
                return entry
        return None


_SCAN_RECORD: ContextVar[PluginRecord | None] = ContextVar("sampler_scan_record", default=None)


@contextmanager
def scan_plugin_context(record: PluginRecord | None):
    """Resolve one scan's saved metadata without changing the user's registry."""
    token = _SCAN_RECORD.set(record)
    try:
        yield
    finally:
        _SCAN_RECORD.reset(token)


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

# Samplers shipped with Jarvis-HEP in the public plug-in format (dogfooding:
# they exercise the same path as user plug-ins). Flags here must match the
# file's PLUGIN dict; a test checks it.
_PLUGIN_DIR = Path(__file__).resolve().parent / "sampling" / "plugins"
_BUILTIN_PLUGINS: tuple[tuple[SamplerSpec, str], ...] = (
    (
        SamplerSpec(
            name="LatinHypercube",
            family=PLUGIN_FAMILY,
            source=SOURCE_BUILTIN_PLUGIN,
        ),
        "latin_hypercube.py",
    ),
)
_BUILTIN_PLUGIN_BY_NAME: dict[str, tuple[SamplerSpec, str]] = {
    spec.name: (spec, filename) for spec, filename in _BUILTIN_PLUGINS
}


def builtin_specs() -> tuple[SamplerSpec, ...]:
    """Samplers implemented in this package (each has a Distributor factory)."""
    return _SPECS


def builtin_plugin_specs() -> tuple[SamplerSpec, ...]:
    return tuple(spec for spec, _ in _BUILTIN_PLUGINS)


def reserved_method_names() -> frozenset[str]:
    """Names a user plug-in may not take (compared case-insensitively)."""
    from jarvis_hep_sampler import metadata

    # The SDK is the sole maintained reserved list. Include the HEP catalog
    # defensively while old SDK installations are still in use.
    sdk_names = getattr(metadata, "RESERVED_METHOD_NAMES", ())
    return frozenset(sdk_names) | frozenset(_BY_NAME) | frozenset(_BUILTIN_PLUGIN_BY_NAME)


# --------------------------------------------------------------- user layer
_USER_CACHE: dict[str, Any] = {"key": None, "specs": {}, "records": {}, "entries": {},
                               "error": None}


def _user_layer() -> tuple[dict[str, SamplerSpec], dict[str, PluginRecord], str | None]:
    """Read the Jarvis-HEP-Sampler registry, cached on its path and file signature."""
    try:
        from jarvis_hep_sampler.registry import SamplerRegistry
    except ImportError as exc:  # pragma: no cover - hard dependency
        from jarvishep2.sampler_sdk import SDK_DISTRIBUTION, SDK_SPECIFIER

        return {}, {}, (
            f"{SDK_DISTRIBUTION}{SDK_SPECIFIER} is not installed "
            f"({exc}). Install it with: python3 -m pip install "
            f'"{SDK_DISTRIBUTION}{SDK_SPECIFIER}"'
        )
    from jarvishep2.sampler_sdk import SamplerSDKError, require_sampler_sdk

    try:
        require_sampler_sdk()
    except SamplerSDKError as exc:
        return {}, {}, str(exc)
    registry = SamplerRegistry()
    try:
        key = (str(registry.home), registry.signature())
    except RegistryError as exc:
        return {}, {}, str(exc)
    if _USER_CACHE["key"] == key:
        return _USER_CACHE["specs"], _USER_CACHE["records"], _USER_CACHE["error"]
    specs: dict[str, SamplerSpec] = {}
    records: dict[str, PluginRecord] = {}
    error: str | None = None
    try:
        entries = registry.entries()
    except RegistryError as exc:
        entries, error = {}, str(exc)
    reserved = {name.casefold() for name in reserved_method_names()}
    for method, entry in entries.items():
        if method.casefold() in reserved:
            continue
        specs[method] = SamplerSpec(
            name=method,
            stateless=entry.profile == "independent",
            statistical=entry.statistical,
            family=PLUGIN_FAMILY,
            source=SOURCE_USER,
        )
        try:
            snapshot_path = str(registry.snapshot_path(entry))
        except RegistryError:
            # Keep the registration visible. The report will identify the bad
            # snapshot and preflight will reject it if this method is selected.
            snapshot_path = ""
        records[method] = PluginRecord(
            method=method,
            source=SOURCE_USER,
            version=entry.version,
            profile=entry.profile,
            description=entry.description,
            statistical=entry.statistical,
            options_schema=entry.options_schema,
            hooks=entry.hooks,
            source_path=entry.source_path,
            snapshot_path=snapshot_path,
            source_sha256=entry.source_sha256,
            loaded_at=entry.loaded_at,
        )
    _USER_CACHE.update(key=key, specs=specs, records=records, entries=entries, error=error)
    return specs, records, error


def registry_report() -> RegistryReport:
    """Read all names and recheck every external snapshot without importing it."""
    from jarvis_hep_sampler.registry import SamplerRegistry

    _, _, error = _user_layer()
    builtins = (*builtin_specs(), *builtin_plugin_specs())
    rows = [
        RegistryReportEntry(spec.name, spec.source, "built-in")
        for spec in builtins
    ]
    if error is not None:
        return RegistryReport(tuple(rows), error)
    registry = SamplerRegistry()
    reserved = {spec.name.casefold(): spec.name for spec in builtins}
    for method, entry in _USER_CACHE["entries"].items():
        problems = tuple(registry.problems(entry))
        builtin = reserved.get(method.casefold())
        if builtin is not None:
            rows.append(RegistryReportEntry(
                method, SOURCE_USER, "name-conflict",
                (f"conflicts with Jarvis-HEP built-in {builtin!r}", *problems),
                f"Jarvis sampler remove {method} --yes",
            ))
        elif problems:
            rows.append(RegistryReportEntry(method, SOURCE_USER, "snapshot-broken", problems))
        else:
            rows.append(RegistryReportEntry(method, SOURCE_USER, "ok"))
    return RegistryReport(tuple(rows))


def user_specs() -> tuple[SamplerSpec, ...]:
    specs, _, _ = _user_layer()
    return tuple(specs[name] for name in sorted(specs))


def user_registry_error() -> str | None:
    """Why the user registry could not be read, or None."""
    _, _, error = _user_layer()
    return error


@lru_cache(maxsize=None)
def builtin_plugin_path(name: str) -> Path | None:
    item = _BUILTIN_PLUGIN_BY_NAME.get(str(name).strip())
    return _PLUGIN_DIR / item[1] if item is not None else None


@lru_cache(maxsize=None)
def _builtin_plugin_record(name: str) -> PluginRecord:
    from jarvis_hep_sampler.loader import load_plugin_file

    path = builtin_plugin_path(name)
    from jarvishep2.sampler_sdk import require_sampler_sdk

    require_sampler_sdk()
    plugin = load_plugin_file(path)
    meta = plugin.metadata
    return PluginRecord(
        method=meta.method,
        source=SOURCE_BUILTIN_PLUGIN,
        version=meta.version,
        profile=meta.profile,
        description=meta.description,
        statistical=meta.statistical,
        options_schema=plugin.options_schema,
        hooks=plugin.hooks,
        source_path=str(path),
        snapshot_path=str(path),
        source_sha256=plugin.source_sha256,
    )


def plugin_record(name: str) -> PluginRecord | None:
    """Plug-in details for *name*, or None for a built-in or unknown method."""
    key = str(name).strip()
    saved = _SCAN_RECORD.get()
    if saved is not None and saved.method == key:
        return saved
    return registered_plugin_record(key)


def registered_plugin_record(name: str) -> PluginRecord | None:
    """Current registration, deliberately ignoring a scan's contextual binding."""
    key = str(name).strip()
    if key in _BUILTIN_PLUGIN_BY_NAME:
        return _builtin_plugin_record(key)
    _, records, _ = _user_layer()
    return records.get(key)


# ------------------------------------------------------------------- lookup
def iter_specs() -> Iterator[SamplerSpec]:
    """Every method: built-ins, built-in plug-ins, then user plug-ins."""
    yield from _SPECS
    saved = _SCAN_RECORD.get()
    for spec in (*builtin_plugin_specs(), *user_specs()):
        if saved is None or spec.name != saved.method:
            yield spec
    if saved is not None:
        yield SamplerSpec(name=saved.method, stateless=saved.profile == "independent",
                          statistical=saved.statistical, family=PLUGIN_FAMILY, source=SOURCE_SCAN)


def spec_for(name: str) -> SamplerSpec | None:
    key = str(name).strip()
    spec = _BY_NAME.get(key)
    if spec is not None:
        return spec
    saved = _SCAN_RECORD.get()
    if saved is not None and saved.method == key:
        return SamplerSpec(name=key, stateless=saved.profile == "independent",
                           statistical=saved.statistical, family=PLUGIN_FAMILY, source=SOURCE_SCAN)
    item = _BUILTIN_PLUGIN_BY_NAME.get(key)
    if item is not None:
        return item[0]
    specs, _, _ = _user_layer()
    return specs.get(key)


def canonical_name(name: str) -> str | None:
    """Exact spelling of a method name matched case-insensitively, or None."""
    folded = str(name).strip().casefold()
    for spec in iter_specs():
        if spec.name.casefold() == folded:
            return spec.name
    return None


def is_known(name: str) -> bool:
    return spec_for(name) is not None


def is_plugin(name: str) -> bool:
    spec = spec_for(name)
    return bool(spec is not None and spec.is_plugin)


def names() -> frozenset[str]:
    return frozenset(spec.name for spec in iter_specs())


def available_methods() -> list[str]:
    """Sorted YAML method names (validation / error messages)."""
    return sorted(names())


def schema_ids() -> dict[str, str]:
    """Bundled method schemas. Only built-ins have one; plug-ins bring their own."""
    return {spec.name: spec.schema_id for spec in _SPECS}


def is_stateless(name: str) -> bool:
    spec = spec_for(name)
    return bool(spec is not None and spec.stateless)


def methods_need_variables() -> frozenset[str]:
    return frozenset(spec.name for spec in iter_specs() if spec.needs_variables)


def nested_names() -> frozenset[str]:
    return frozenset(spec.name for spec in iter_specs() if spec.nested)


def mcmc_names() -> frozenset[str]:
    return frozenset(spec.name for spec in iter_specs() if spec.mcmc)


def statistical_names() -> frozenset[str]:
    return frozenset(spec.name for spec in iter_specs() if spec.statistical)


def is_statistical(name: str) -> bool:
    spec = spec_for(name)
    return bool(spec is not None and spec.statistical)


def mcmc_independent_names() -> frozenset[str]:
    return frozenset(
        spec.name for spec in iter_specs() if spec.mcmc_pipeline == "async_independent"
    )


def mcmc_coupled_names() -> frozenset[str]:
    return frozenset(
        spec.name for spec in iter_specs() if spec.mcmc_pipeline == "barrier_coupled"
    )


def family_of(name: str) -> str:
    spec = spec_for(name)
    return spec.family if spec is not None else "simple"


__all__ = [
    "PLUGIN_FAMILY",
    "PluginRecord",
    "RegistryReport",
    "RegistryReportEntry",
    "SamplerRegistryPreflightError",
    "SOURCE_BUILTIN",
    "SOURCE_BUILTIN_PLUGIN",
    "SOURCE_SCAN",
    "SOURCE_USER",
    "SamplerSpec",
    "available_methods",
    "builtin_plugin_path",
    "builtin_plugin_specs",
    "builtin_specs",
    "canonical_name",
    "family_of",
    "is_known",
    "is_plugin",
    "is_stateless",
    "is_statistical",
    "iter_specs",
    "mcmc_coupled_names",
    "mcmc_independent_names",
    "mcmc_names",
    "methods_need_variables",
    "names",
    "nested_names",
    "plugin_record",
    "registered_plugin_record",
    "registry_report",
    "reserved_method_names",
    "schema_ids",
    "scan_plugin_context",
    "spec_for",
    "statistical_names",
    "user_registry_error",
    "user_specs",
]
