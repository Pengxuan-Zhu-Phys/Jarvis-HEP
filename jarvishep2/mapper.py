#!/usr/bin/env python3
"""u → physical-parameter mappers for control process and Workers (D22).

Single implementation path: :class:`MapperPipeline` applies
``Variables[].distribution`` first, then optional ``Sampling.Mapper``
list entries (``{name, expression}``). Control-side selection and
Worker-side ``bind_params`` share the same pipeline so uuid ↔ coordinate
binding stays pure under D21 resume replay.
"""

from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import sympy as sp

from jarvishep2.Sampling.variables import Variable
from jarvishep2.expression import ExpressionContext
from jarvishep2.inner_func import EXPRESSION_CONSTANTS
from jarvishep2.sample import UMapperProtocol

# DATABASE / Sample reserved column names (must not appear as derive targets).
_DATABASE_RESERVED: frozenset[str] = frozenset(
    {"uuid", "sample_index", "status", "LogL"}
)
_CONSTANT_RESERVED: frozenset[str] = frozenset(EXPRESSION_CONSTANTS)

# Statistical samplers: dimension expansion / nonlinear reparameterization
# warnings (JV2-MAP-050 / 051). Coverage-style methods stay silent.
STATISTICAL_METHODS: frozenset[str] = frozenset(
    {
        "Dynesty",
        "MultiNest",
        "MCMC",
        "ToyMCMC",
        "AMMCMC",
        "AM",
        "DRAM",
        "EnsembleMCMC",
        "Ensemble",
        "DEMCMC",
        "PTMCMC",
        "PTEnsemble",
    }
)


@dataclass(frozen=True)
class VariableSpec:
    """Picklable sampling-variable descriptor (no callables)."""

    name: str
    distribution: str
    parameters: Mapping[str, Any]
    description: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "distribution": {
                "type": self.distribution,
                "parameters": dict(self.parameters),
            },
        }

    @classmethod
    def from_mapping(cls, item: Mapping[str, Any]) -> "VariableSpec":
        name = str(item.get("name", "")).strip()
        dist_block = item.get("distribution")
        if not isinstance(dist_block, Mapping):
            dist_block = {}
        return cls(
            name=name,
            description=str(item.get("description", "") or ""),
            distribution=str(dist_block.get("type", "Flat")).strip() or "Flat",
            parameters=dict(dist_block.get("parameters") or {}),
        )


@dataclass(frozen=True)
class MapperSpec:
    """Pure data for u→params; picklable; compiled callables rebuilt per process."""

    variables: tuple[VariableSpec, ...]
    derive_order: tuple[str, ...] = ()
    derive_exprs: Mapping[str, str] = field(default_factory=dict)
    # name → free symbols used by that expression (subset of namespace at eval time)
    derive_symbols: Mapping[str, tuple[str, ...]] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "variables": [var.to_dict() for var in self.variables],
            "derive_order": list(self.derive_order),
            "derive_exprs": dict(self.derive_exprs),
            "derive_symbols": {
                name: list(syms) for name, syms in self.derive_symbols.items()
            },
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "MapperSpec":
        variables = tuple(
            VariableSpec.from_mapping(item)
            for item in (payload.get("variables") or [])
            if isinstance(item, Mapping)
        )
        derive_order = tuple(str(n) for n in (payload.get("derive_order") or ()))
        derive_exprs = {
            str(k): str(v)
            for k, v in dict(payload.get("derive_exprs") or {}).items()
        }
        raw_syms = payload.get("derive_symbols") or {}
        derive_symbols = {
            str(k): tuple(str(s) for s in (v or ()))
            for k, v in dict(raw_syms).items()
        }
        return cls(
            variables=variables,
            derive_order=derive_order,
            derive_exprs=derive_exprs,
            derive_symbols=derive_symbols,
        )

    def fingerprint_payload(self) -> dict[str, Any]:
        """Stable subset used for checkpoint mapper_hash (Mapper + var names)."""
        return {
            "variables": [
                {
                    "name": v.name,
                    "distribution": v.distribution,
                    "parameters": dict(v.parameters),
                }
                for v in self.variables
            ],
            # Sorted by name so hash is order-insensitive for equivalent sets;
            # evaluation order still follows derive_order / DAG.
            "mapper": [
                {"name": name, "expression": self.derive_exprs[name]}
                for name in sorted(self.derive_order)
            ],
        }


class MapperError(ValueError):
    """Load-time mapper configuration error with optional JV2-MAP code."""

    def __init__(
        self,
        message: str,
        *,
        code: str = "JV2-MAP-001",
        path: str = "Sampling.Mapper",
    ) -> None:
        super().__init__(message)
        self.code = code
        self.path = path
        self.message = message


def _variable_specs_from_config(config: Mapping[str, Any]) -> tuple[VariableSpec, ...]:
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    raw = sampling.get("Variables") or []
    if not isinstance(raw, list):
        return ()
    specs: list[VariableSpec] = []
    for item in raw:
        if not isinstance(item, Mapping):
            continue
        name = str(item.get("name", "")).strip()
        if not name:
            continue
        specs.append(VariableSpec.from_mapping(item))
    return tuple(specs)


def _mapper_needs_operas_context(derive_exprs: Mapping[str, str]) -> bool:
    """True when any Mapper expression contains a qualified Operas reference."""
    if not derive_exprs:
        return False
    try:
        from jarvishep2.operas_functions import _string_uses_operas_reference
    except Exception:  # pragma: no cover
        return False
    return any(_string_uses_operas_reference(str(text)) for text in derive_exprs.values())


def _expression_context_for_mapper_exprs(
    derive_exprs: Mapping[str, str],
    *,
    expression_context: ExpressionContext | None = None,
) -> ExpressionContext:
    """Return a context that can resolve Operas constants/functions if needed."""
    if expression_context is not None:
        return expression_context
    if not _mapper_needs_operas_context(derive_exprs):
        return ExpressionContext()
    try:
        from jarvishep2.operas_functions import build_operas_expression_context

        return build_operas_expression_context(required=True)
    except ImportError as exc:
        raise MapperError(
            "Mapper expressions reference Jarvis-Operas names but Jarvis-Operas "
            f"is not importable: {exc}",
            code="JV2-MAP-007",
            path="Sampling.Mapper",
        ) from exc


def _compile_mapper_error_message(name: str, expr_text: str, exc: BaseException) -> str:
    """Prefer actionable Operas-namespace guidance over a bare sympy traceback."""
    detail = str(exc)
    # D23.9: bare ExpressionContext turns ``pdg.mZ`` into Symbol.pdg AttributeError.
    if "has no attribute" in detail or "UndefinedFunction" in detail:
        return (
            f"cannot compile Mapper expression for {name!r}: {exc}; "
            f"expression={expr_text!r}. If this uses a Jarvis-Operas name "
            f"(e.g. pdg.mZ or math.add(...)), ensure Jarvis-Operas is installed "
            f"and the name is registered; constants use bare form without ()."
        )
    return (
        f"cannot compile Mapper expression for {name!r}: {exc}; "
        f"expression={expr_text!r}"
    )


def _topo_sort_derive(
    derive_exprs: Mapping[str, str],
    *,
    sample_var_names: frozenset[str],
    expression_context: ExpressionContext | None = None,
) -> tuple[tuple[str, ...], Mapping[str, tuple[str, ...]]]:
    """Return (evaluate order, name→free-symbol tuple) or raise MapperError."""
    if not derive_exprs:
        return (), {}

    ctx = _expression_context_for_mapper_exprs(
        derive_exprs, expression_context=expression_context
    )
    all_derive = frozenset(derive_exprs)
    # Symbols that may appear: sampling vars + all derive names (DAG may reorder).
    declared = tuple(sorted(sample_var_names | all_derive))
    free_by_name: dict[str, tuple[str, ...]] = {}
    deps: dict[str, set[str]] = {}

    for name, expr_text in derive_exprs.items():
        try:
            compiled = ctx.compile(str(expr_text), symbols=declared)
        except Exception as exc:
            raise MapperError(
                _compile_mapper_error_message(name, str(expr_text), exc),
                code="JV2-MAP-007",
                path=f"Sampling.Mapper.{name}",
            ) from exc
        free = tuple(compiled.variable_names)
        free_by_name[name] = free
        unknown = [s for s in free if s not in sample_var_names and s not in all_derive]
        if unknown:
            available = sorted(sample_var_names | all_derive)
            raise MapperError(
                f"Mapper entry {name!r} references unknown symbol(s) {unknown}; "
                f"available: {available}",
                code="JV2-MAP-002",
                path=f"Sampling.Mapper.{name}",
            )
        # Depend only on other derive names (sample vars are always present).
        deps[name] = {s for s in free if s in all_derive and s != name}

    # Kahn topological sort; stable by original key order among ready nodes.
    original_order = list(derive_exprs.keys())
    indegree = {n: 0 for n in derive_exprs}
    dependents: dict[str, list[str]] = defaultdict(list)
    for name, dep_set in deps.items():
        for dep in dep_set:
            dependents[dep].append(name)
            indegree[name] += 1

    ready = deque([n for n in original_order if indegree[n] == 0])
    order: list[str] = []
    while ready:
        node = ready.popleft()
        order.append(node)
        for child in dependents[node]:
            indegree[child] -= 1
            if indegree[child] == 0:
                # preserve relative original order among newly ready
                ready.append(child)
        # re-sort ready by original order for determinism
        if len(ready) > 1:
            ready = deque(sorted(ready, key=lambda n: original_order.index(n)))

    if len(order) != len(derive_exprs):
        remaining = [n for n, d in indegree.items() if d > 0]
        raise MapperError(
            f"Mapper dependency cycle involving: {remaining}",
            code="JV2-MAP-004",
            path="Sampling.Mapper",
        )
    return tuple(order), free_by_name


def _is_affine_in_sample_vars(expr_text: str, sample_vars: Sequence[str]) -> bool:
    """True when expression is affine-linear in the sampling variables."""
    if not sample_vars:
        return True
    try:
        symbols = {name: sp.Symbol(name) for name in sample_vars}
        parse_locals = dict(symbols)
        parse_locals.update(EXPRESSION_CONSTANTS)
        parsed = sp.sympify(str(expr_text), locals=parse_locals)
        free = [symbols[n] for n in sample_vars if symbols[n] in parsed.free_symbols]
        if not free:
            return True
        poly = sp.Poly(sp.expand(parsed), *free)
        return poly.total_degree() <= 1
    except Exception:
        return False


def _parse_mapper_list(mapper_block: Any) -> dict[str, str]:
    """Parse ``Sampling.Mapper`` list of ``{name, expression}`` into ordered dict.

    List form keeps JSON Schema keys fixed (``name`` / ``expression``) so user
    symbols never become object property names.
    """
    if isinstance(mapper_block, Mapping):
        raise MapperError(
            "Sampling.Mapper must be a list of {name, expression} mappings, "
            "not a free-form name → expression object "
            '(example: Mapper: [{name: x, expression: "cos(t)"}])',
            code="JV2-MAP-001",
            path="Sampling.Mapper",
        )
    if not isinstance(mapper_block, list):
        raise MapperError(
            f"Sampling.Mapper must be a list of {{name, expression}} mappings, "
            f"got {type(mapper_block).__name__}",
            code="JV2-MAP-001",
            path="Sampling.Mapper",
        )

    derive_exprs: dict[str, str] = {}
    for index, item in enumerate(mapper_block):
        path = f"Sampling.Mapper[{index}]"
        if not isinstance(item, Mapping):
            raise MapperError(
                f"Mapper entry must be a mapping with name and expression, "
                f"got {type(item).__name__}",
                code="JV2-MAP-001",
                path=path,
            )
        name = str(item.get("name") or "").strip()
        if not name:
            raise MapperError(
                "Mapper entry name must be a non-empty string",
                code="JV2-MAP-001",
                path=f"{path}.name",
            )
        raw_expr = item.get("expression")
        if not isinstance(raw_expr, str) or not str(raw_expr).strip():
            raise MapperError(
                f"Mapper entry {name!r} expression must be a non-empty string",
                code="JV2-MAP-001",
                path=f"{path}.expression",
            )
        if name in derive_exprs:
            raise MapperError(
                f"duplicate Mapper name {name!r}",
                code="JV2-MAP-003",
                path=f"{path}.name",
            )
        derive_exprs[name] = str(raw_expr).strip()
    return derive_exprs


def build_mapper_spec_from_config(config: Mapping[str, Any]) -> MapperSpec:
    """Build a :class:`MapperSpec` from a task card (no Mapper ⇒ empty expressions).

    Raises :class:`MapperError` for structural / closed-namespace failures.
    Callers that only need distribution mapping with no YAML Mapper still get
    a valid empty-expression spec.

    YAML shape (list of fixed-key mappings)::

        Sampling:
          Mapper:
            - name: x
              expression: "cos(t)"
            - name: y
              expression: "sin(t)"
    """
    variables = _variable_specs_from_config(config)
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    method = str(sampling.get("Method") or "").strip()
    mapper_block = sampling.get("Mapper") if isinstance(sampling, Mapping) else None

    if mapper_block is None:
        return MapperSpec(variables=variables)

    if method == "CSV":
        raise MapperError(
            "Sampling.Mapper is not supported for Method CSV in v1 "
            "(CSV parameters come from columns; Worker uses type=none)",
            code="JV2-MAP-010",
            path="Sampling.Mapper",
        )

    derive_exprs = _parse_mapper_list(mapper_block)
    if not derive_exprs:
        return MapperSpec(variables=variables)

    sample_names = [v.name for v in variables]
    sample_set = frozenset(sample_names)

    for name in derive_exprs:
        if name in sample_set:
            raise MapperError(
                f"Mapper name {name!r} conflicts with Sampling.Variables name",
                code="JV2-MAP-003",
                path=f"Sampling.Mapper",
            )
        if name in _CONSTANT_RESERVED:
            raise MapperError(
                f"Mapper name {name!r} is a reserved expression constant "
                f"({sorted(_CONSTANT_RESERVED)})",
                code="JV2-MAP-005",
                path="Sampling.Mapper",
            )
        if name in _DATABASE_RESERVED:
            raise MapperError(
                f"Mapper name {name!r} is a reserved DATABASE column "
                f"({sorted(_DATABASE_RESERVED)})",
                code="JV2-MAP-006",
                path="Sampling.Mapper",
            )

    order, free_by_name = _topo_sort_derive(
        derive_exprs, sample_var_names=sample_set
    )
    return MapperSpec(
        variables=variables,
        derive_order=order,
        derive_exprs=dict(derive_exprs),
        derive_symbols=free_by_name,
    )


def mapper_block_fingerprint(config: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalized Mapper (+ variable names) payload for checkpoint hashing."""
    if not isinstance(config, Mapping):
        return {}
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    mapper = sampling.get("Mapper") if isinstance(sampling, Mapping) else None
    variables = []
    for item in sampling.get("Variables") or [] if isinstance(sampling, Mapping) else []:
        if isinstance(item, Mapping) and str(item.get("name") or "").strip():
            dist = item.get("distribution") if isinstance(item.get("distribution"), Mapping) else {}
            variables.append(
                {
                    "name": str(item.get("name")).strip(),
                    "distribution": str(dist.get("type") or "Flat"),
                    "parameters": dict(dist.get("parameters") or {}),
                }
            )
    entries: list[dict[str, str]] = []
    if isinstance(mapper, list):
        for item in mapper:
            if not isinstance(item, Mapping):
                continue
            name = str(item.get("name") or "").strip()
            expr = item.get("expression")
            if name and isinstance(expr, str) and expr.strip():
                entries.append({"name": name, "expression": expr.strip()})
        entries.sort(key=lambda e: e["name"])
    return {"variables": variables, "mapper": entries}


class MapperPipeline:
    """Single u→params implementation for control process and Workers."""

    def __init__(
        self,
        spec: MapperSpec,
        *,
        expression_context: ExpressionContext | None = None,
    ) -> None:
        self._spec = spec
        self._variables: list[Variable] = [
            Variable(
                name=v.name,
                description=v.description,
                distribution=v.distribution,
                parameters=dict(v.parameters),
            )
            for v in spec.variables
        ]
        # D23.9: inject Operas-aware context when Mapper expressions need it
        # (pdg.mZ, math.add, …). A bare ExpressionContext treats namespaces as
        # plain Symbols and fails with a misleading AttributeError.
        self._expression_context = _expression_context_for_mapper_exprs(
            spec.derive_exprs, expression_context=expression_context
        )
        self._compiled: dict[str, Any] = {}
        # Pre-declare symbols available at each derive step for cache keys.
        known: list[str] = [v.name for v in self._variables]
        for name in spec.derive_order:
            expr = spec.derive_exprs[name]
            # Namespace at evaluation: sample vars + prior derives.
            symbols = tuple(known)
            try:
                self._compiled[name] = self._expression_context.compile(
                    expr, symbols=symbols
                )
            except Exception as exc:
                raise MapperError(
                    _compile_mapper_error_message(name, expr, exc),
                    code="JV2-MAP-007",
                    path=f"Sampling.Mapper.{name}",
                ) from exc
            known.append(name)
        self._output_names = tuple(known)

    @classmethod
    def from_spec(
        cls,
        spec: MapperSpec,
        *,
        expression_context: ExpressionContext | None = None,
    ) -> "MapperPipeline":
        return cls(spec, expression_context=expression_context)

    @classmethod
    def from_config(
        cls,
        config: Mapping[str, Any],
        *,
        expression_context: ExpressionContext | None = None,
    ) -> "MapperPipeline":
        return cls.from_spec(
            build_mapper_spec_from_config(config),
            expression_context=expression_context,
        )

    @property
    def spec(self) -> MapperSpec:
        return self._spec

    @property
    def output_names(self) -> tuple[str, ...]:
        return self._output_names

    @property
    def variable_names(self) -> tuple[str, ...]:
        return tuple(v.name for v in self._variables)

    @property
    def derive_names(self) -> tuple[str, ...]:
        return self._spec.derive_order

    def map(self, u_coords: np.ndarray | Sequence[float]) -> dict[str, float]:
        coords = np.asarray(u_coords, dtype=np.float64).reshape(-1)
        n_vars = len(self._variables)
        if n_vars == 0:
            return {}
        if len(coords) != n_vars:
            raise ValueError(
                f"u_coords length {len(coords)} must equal variable count {n_vars}"
            )
        mapped: dict[str, float] = {}
        for index, var in enumerate(self._variables):
            value = var.map_standard_random_to_distribution(float(coords[index]))
            mapped[var.name] = float(value) if not isinstance(value, int) else value
        for name in self._spec.derive_order:
            compiled = self._compiled[name]
            raw = compiled.evaluate(mapped)
            mapped[name] = float(np.asarray(raw, dtype=np.float64).reshape(-1)[0])
        return mapped


class DistributionUMapper:
    """Map u_coords in [0, 1] using V1-compatible variable distributions."""

    def __init__(self, variables: Sequence[Mapping[str, Any]]) -> None:
        self._variables: list[Variable] = []
        for var in variables:
            if not isinstance(var, Mapping):
                continue
            name = str(var.get("name", "")).strip()
            if not name:
                continue
            dist_block = dict(var.get("distribution") or {})
            self._variables.append(
                Variable(
                    name=name,
                    description=str(var.get("description", "")),
                    distribution=str(dist_block.get("type", "Flat")).strip(),
                    parameters=dict(dist_block.get("parameters") or {}),
                )
            )

    def map(self, u_coords: np.ndarray) -> dict[str, float]:
        coords = np.asarray(u_coords, dtype=np.float64).reshape(-1)
        if len(coords) != len(self._variables):
            raise ValueError(
                f"u_coords length {len(coords)} must equal variable count "
                f"{len(self._variables)}"
            )
        mapped: dict[str, float] = {}
        for index, var in enumerate(self._variables):
            mapped[var.name] = var.map_standard_random_to_distribution(float(coords[index]))
        return mapped


class FlatUMapper:
    """Map normalized u_coords in [0, 1] to flat-distributed physical parameters.

    Internal / test helper — not reachable from task-card YAML (D22).
    """

    def __init__(
        self,
        variables: Sequence[Mapping[str, Any]],
    ) -> None:
        self._names: list[str] = []
        self._bounds: list[tuple[float, float]] = []
        for var in variables:
            name = str(var.get("name", "")).strip()
            if not name:
                continue
            params = var.get("distribution", {}).get("parameters", {}) or {}
            lo = float(params.get("min", 0.0))
            hi = float(params.get("max", 1.0))
            self._names.append(name)
            self._bounds.append((lo, hi))

    def map(self, u_coords: np.ndarray) -> dict[str, float]:
        coords = np.asarray(u_coords, dtype=np.float64).reshape(-1)
        if len(coords) != len(self._names):
            raise ValueError(
                f"u_coords length {len(coords)} must equal variable count {len(self._names)}"
            )
        mapped: dict[str, float] = {}
        for index, name in enumerate(self._names):
            lo, hi = self._bounds[index]
            mapped[name] = float(lo + coords[index] * (hi - lo))
        return mapped


class IdentityParamMapper:
    """Pass through params already embedded in opera_params (test helper)."""

    def __init__(self, keys: Sequence[str] | None = None) -> None:
        self._keys = tuple(keys or ())

    def map(self, u_coords: np.ndarray) -> dict[str, float]:
        coords = np.asarray(u_coords, dtype=np.float64).reshape(-1)
        if not self._keys:
            return {"u": float(coords[0]) if coords.size else 0.0}
        if len(coords) != len(self._keys):
            raise ValueError(
                f"u_coords length {len(coords)} must equal identity key count "
                f"{len(self._keys)}"
            )
        return {key: float(coords[index]) for index, key in enumerate(self._keys)}


def build_mapper(config: Mapping[str, Any] | MapperSpec | None) -> UMapperProtocol | None:
    """Construct a mapper from picklable Worker config or :class:`MapperSpec`."""
    if config is None:
        return None
    if isinstance(config, MapperSpec):
        return MapperPipeline.from_spec(config)
    if not isinstance(config, Mapping):
        return None
    # Pipeline payload produced by worker_config.
    if config.get("type") == "pipeline" or "derive_order" in config or (
        "variables" in config and "derive_exprs" in config
    ):
        if config.get("type") == "pipeline":
            payload = config.get("spec") or {}
            if not isinstance(payload, Mapping):
                return None
            return MapperPipeline.from_spec(MapperSpec.from_dict(payload))
        return MapperPipeline.from_spec(MapperSpec.from_dict(config))
    mapper_type = str(config.get("type", "flat")).strip().lower()
    if mapper_type == "none":
        return None
    if mapper_type == "identity":
        return IdentityParamMapper(config.get("keys") or ())
    if mapper_type == "distribution":
        variables = config.get("variables") or []
        if not variables:
            return None
        # Prefer pipeline so Worker path matches control when derive is empty.
        specs = tuple(
            VariableSpec.from_mapping(item)
            for item in variables
            if isinstance(item, Mapping)
        )
        return MapperPipeline.from_spec(MapperSpec(variables=specs))
    variables = config.get("variables") or []
    if not variables:
        return None
    return FlatUMapper(variables)


__all__ = [
    "STATISTICAL_METHODS",
    "DistributionUMapper",
    "FlatUMapper",
    "IdentityParamMapper",
    "MapperError",
    "MapperPipeline",
    "MapperSpec",
    "VariableSpec",
    "build_mapper",
    "build_mapper_spec_from_config",
    "mapper_block_fingerprint",
    "_is_affine_in_sample_vars",
]
