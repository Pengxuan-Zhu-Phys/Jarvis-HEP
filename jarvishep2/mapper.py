#!/usr/bin/env python3
"""Worker-held u → x mapper implementations."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np

from jarvishep2.Sampling.variables import Variable
from jarvishep2.sample import UMapperProtocol


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
        if len(coords) < len(self._variables):
            raise ValueError(
                f"u_coords length {len(coords)} is smaller than variable count {len(self._variables)}"
            )
        mapped: dict[str, float] = {}
        for index, var in enumerate(self._variables):
            mapped[var.name] = var.map_standard_random_to_distribution(float(coords[index]))
        return mapped


class FlatUMapper:
    """Map normalized u_coords in [0, 1] to flat-distributed physical parameters."""

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
        if len(coords) < len(self._names):
            raise ValueError(
                f"u_coords length {len(coords)} is smaller than variable count {len(self._names)}"
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
        if len(coords) < len(self._keys):
            raise ValueError("u_coords shorter than configured identity keys")
        return {key: float(coords[index]) for index, key in enumerate(self._keys)}


def build_mapper(config: Mapping[str, Any] | None) -> UMapperProtocol | None:
    """Construct a mapper from picklable Worker config."""
    if not isinstance(config, Mapping):
        return None
    mapper_type = str(config.get("type", "flat")).strip().lower()
    if mapper_type == "none":
        return None
    if mapper_type == "identity":
        return IdentityParamMapper(config.get("keys") or ())
    if mapper_type == "distribution":
        variables = config.get("variables") or []
        if variables:
            return DistributionUMapper(variables)
        return None
    variables = config.get("variables") or []
    if not variables:
        return None
    return FlatUMapper(variables)


__all__ = ["DistributionUMapper", "FlatUMapper", "IdentityParamMapper", "build_mapper"]