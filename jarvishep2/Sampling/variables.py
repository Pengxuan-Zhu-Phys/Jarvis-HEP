#!/usr/bin/env python3
"""Sampling variables and distribution mapping (V1-compatible)."""

from __future__ import annotations

from statistics import NormalDist

import numpy as np

_STANDARD_NORMAL = NormalDist()


def _unit_interval(value: float) -> float:
    return float(np.clip(float(value), np.finfo(float).tiny, 1.0 - np.finfo(float).eps))


class Variable:
    def __init__(self, name: str, description: str, distribution: str, parameters: dict) -> None:
        self._name = name
        self._description = description
        self._distribution = distribution
        self._parameters = parameters

    @property
    def name(self) -> str:
        return self._name

    @property
    def description(self) -> str:
        return self._description

    @property
    def distribution(self) -> str:
        return self._distribution

    @property
    def parameters(self) -> dict:
        return self._parameters

    def map_standard_random_to_distribution(self, std_rand: float) -> float:
        if self.distribution == "Flat":
            lo = self.parameters["min"]
            hi = self.parameters["max"]
            return float(lo + (hi - lo) * std_rand)
        if self.distribution == "Log":
            log_min = np.log(self.parameters["min"])
            log_max = np.log(self.parameters["max"])
            return float(np.exp(log_min + (log_max - log_min) * std_rand))
        if self.distribution == "Normal":
            mean = self.parameters["mean"]
            stddev = self.parameters["stddev"]
            return float(mean + stddev * _STANDARD_NORMAL.inv_cdf(_unit_interval(std_rand)))
        if self.distribution == "Log-Normal":
            mean = self.parameters["mean"]
            stddev = self.parameters["stddev"]
            return float(np.exp(mean + stddev * _STANDARD_NORMAL.inv_cdf(_unit_interval(std_rand))))
        if self.distribution == "Logit":
            location = self.parameters.get("location", 0)
            scale = self.parameters.get("scale", 1)
            prob = _unit_interval(std_rand)
            return float((np.log(prob) - np.log1p(-prob)) * scale + location)
        raise ValueError(f"Unsupported distribution type: {self.distribution}")


def load_variables(config: dict) -> list[Variable]:
    sampling = dict(config.get("Sampling") or {})
    raw_vars = sampling.get("Variables") or []
    variables: list[Variable] = []
    for item in raw_vars:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", "")).strip()
        if not name:
            continue
        dist_block = dict(item.get("distribution") or {})
        dist_type = str(dist_block.get("type", "Flat")).strip()
        parameters = dict(dist_block.get("parameters") or {})
        variables.append(
            Variable(
                name=name,
                description=str(item.get("description", "")),
                distribution=dist_type,
                parameters=parameters,
            )
        )
    if not variables:
        raise ValueError("Sampling.Variables must define at least one variable")
    return variables


__all__ = ["Variable", "load_variables"]