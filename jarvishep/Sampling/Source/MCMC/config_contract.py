#!/usr/bin/env python3
from __future__ import annotations

from typing import Any, Iterable, Mapping


_MISSING = object()


def _dotted_get(mapping: Mapping[str, Any], dotted_key: str) -> Any:
    cur: Any = mapping
    for part in dotted_key.split("."):
        if not isinstance(cur, Mapping) or part not in cur:
            return _MISSING
        cur = cur[part]
    return cur


def bounds_get(
    bounds: Mapping[str, Any],
    key: str,
    *,
    aliases: Iterable[str] = (),
    default: Any = _MISSING,
) -> Any:
    for candidate in (key, *aliases):
        if "." in candidate:
            value = _dotted_get(bounds, candidate)
            if value is not _MISSING:
                return value
            continue
        if candidate in bounds:
            return bounds[candidate]
    if default is not _MISSING:
        return default
    raise KeyError(f"Missing required bounds key: {key}")


def bounds_get_int(
    bounds: Mapping[str, Any],
    key: str,
    *,
    aliases: Iterable[str] = (),
    default: Any = _MISSING,
    minimum: int | None = None,
) -> int:
    value = bounds_get(bounds, key, aliases=aliases, default=default)
    ivalue = int(value)
    if minimum is not None and ivalue < int(minimum):
        raise ValueError(f"{key} must be >= {int(minimum)}, got {ivalue}")
    return ivalue


def bounds_get_float(
    bounds: Mapping[str, Any],
    key: str,
    *,
    aliases: Iterable[str] = (),
    default: Any = _MISSING,
    minimum: float | None = None,
) -> float:
    value = bounds_get(bounds, key, aliases=aliases, default=default)
    fvalue = float(value)
    if minimum is not None and fvalue < float(minimum):
        raise ValueError(f"{key} must be >= {float(minimum)}, got {fvalue}")
    return fvalue


def bounds_get_bool(
    bounds: Mapping[str, Any],
    key: str,
    *,
    aliases: Iterable[str] = (),
    default: Any = _MISSING,
) -> bool:
    value = bounds_get(bounds, key, aliases=aliases, default=default)
    return bool(value)


def bounds_get_list(
    bounds: Mapping[str, Any],
    key: str,
    *,
    aliases: Iterable[str] = (),
    default: Any = _MISSING,
) -> list[Any]:
    value = bounds_get(bounds, key, aliases=aliases, default=default)
    if isinstance(value, list):
        return list(value)
    if isinstance(value, tuple):
        return list(value)
    return [value]


def parse_common_chain_counts(bounds: Mapping[str, Any]) -> tuple[int, int]:
    nchains = bounds_get_int(bounds, "num_chains", aliases=("chains",), minimum=1)
    niters = bounds_get_int(bounds, "num_iters", aliases=("iterations",), minimum=1)
    return nchains, niters


def parse_proposal_scale_value(bounds: Mapping[str, Any], *, default: float = 0.1) -> float | list[float]:
    scalar = _MISSING
    if "proposal_scale" in bounds:
        scalar = bounds["proposal_scale"]
    else:
        scalar = _dotted_get(bounds, "proposal.scale")

    scales = _MISSING
    if "proposal_scales" in bounds:
        scales = bounds["proposal_scales"]
    else:
        scales = _dotted_get(bounds, "proposal.scales")

    if scalar is _MISSING and scales is _MISSING:
        return float(default)
    if scales is _MISSING:
        return float(scalar)
    if scalar is not _MISSING and not isinstance(scales, (list, tuple)):
        # Explicit scalar + non-list fallback keeps scalar behavior.
        return float(scalar)

    values = [float(x) for x in (scales if isinstance(scales, (list, tuple)) else [scales])]
    if len(values) == 1 and scalar is _MISSING:
        return float(values[0])
    if scalar is not _MISSING and len(values) <= 1:
        return float(scalar)
    return values


def normalize_proposal_scales(
    value: float | list[float],
    *,
    nchains: int,
    sampler_method: str,
) -> list[float]:
    if isinstance(value, (int, float)):
        return [float(value) for _ in range(int(nchains))]
    values = [float(x) for x in value]
    if len(values) == 1:
        return [values[0] for _ in range(int(nchains))]
    if len(values) != int(nchains):
        raise ValueError(
            f"proposal_scale size mismatch for {sampler_method}: "
            f"expect 1 or {int(nchains)}, got {len(values)}"
        )
    return values
