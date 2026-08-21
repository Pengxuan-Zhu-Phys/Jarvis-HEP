#!/usr/bin/env python3
"""Sampling variables and distribution mapping (V1-compatible)."""

from __future__ import annotations

from statistics import NormalDist

import numpy as np
from scipy import special, stats

_STANDARD_NORMAL = NormalDist()


def _unit_interval(value: float) -> float:
    return float(np.clip(float(value), np.finfo(float).tiny, 1.0 - np.finfo(float).eps))


class Variable:
    def __init__(self, name: str, description: str, distribution: str, parameters: dict) -> None:
        self._name = name
        self._description = description
        self._distribution = distribution
        # JSON Schema deliberately accepts numeric YAML strings for V1 cards.
        # Normalize only those values here, so a validated ``"0.05"`` reaches
        # the exact same mapper semantics as YAML's native ``0.05``.
        self._parameters = {
            key: _numeric_string_to_float(value) for key, value in parameters.items()
        }

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

    def generate(self) -> float | int:
        """Draw one V1-compatible value from this distribution.

        V2's distributed samplers normally use
        :meth:`map_standard_random_to_distribution` so that samples can be
        replayed from their deterministic unit-cube coordinates.  Keeping this
        direct-draw helper makes the complete V1 ``Variable`` surface
        available for local callers as well.
        """
        if self.distribution == "Flat":
            return float(np.random.uniform(self.parameters["min"], self.parameters["max"]))
        if self.distribution == "Log":
            return float(
                np.exp(
                    np.random.uniform(
                        np.log(self.parameters["min"]), np.log(self.parameters["max"])
                    )
                )
            )
        if self.distribution == "Normal":
            return float(np.random.normal(self.parameters["mean"], self.parameters["stddev"]))
        if self.distribution == "Log-Normal":
            return float(np.random.lognormal(self.parameters["mean"], self.parameters["stddev"]))
        if self.distribution == "Logit":
            prob = float(np.random.uniform(0.0, 1.0))
            return float(
                self.parameters.get("location", 0.0)
                + self.parameters.get("scale", 1.0) * np.log(prob / (1.0 - prob))
            )
        if self.distribution == "Binomial":
            return int(np.random.binomial(int(self.parameters["n"]), self.parameters["p"]))
        if self.distribution == "Poisson":
            return int(np.random.poisson(self.parameters["lambda"]))
        if self.distribution == "Beta":
            return float(np.random.beta(self.parameters["alpha"], self.parameters["beta"]))
        if self.distribution == "Exponential":
            return float(np.random.exponential(scale=1.0 / self.parameters["rate"]))
        if self.distribution == "Gamma":
            return float(np.random.gamma(self.parameters["shape"], scale=self.parameters["scale"]))
        raise ValueError(f"Unsupported distribution type: {self.distribution}")

    def map_standard_random_to_distribution(self, std_rand: float) -> float | int:
        """Map a deterministic unit-cube coordinate through an inverse CDF.

        The five discrete/non-Gaussian V1 distributions need this V2-specific
        mapping because Workers receive deterministic ``u`` coordinates rather
        than an independent numpy RNG stream.
        """
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
        if self.distribution == "Binomial":
            return int(
                stats.binom.ppf(
                    _unit_interval(std_rand), int(self.parameters["n"]), self.parameters["p"]
                )
            )
        if self.distribution == "Poisson":
            return int(stats.poisson.ppf(_unit_interval(std_rand), self.parameters["lambda"]))
        if self.distribution == "Beta":
            return float(
                special.betaincinv(
                    self.parameters["alpha"], self.parameters["beta"], _unit_interval(std_rand)
                )
            )
        if self.distribution == "Exponential":
            return float(-np.log1p(-_unit_interval(std_rand)) / self.parameters["rate"])
        if self.distribution == "Gamma":
            return float(
                special.gammaincinv(self.parameters["shape"], _unit_interval(std_rand))
                * self.parameters["scale"]
            )
        raise ValueError(f"Unsupported distribution type: {self.distribution}")


def load_variables(config: dict) -> list[Variable]:
    """Build variables from ``Sampling.Variables``.

    Expects cards that already passed :func:`jarvishep2.task_validation.validate_task_config`.
    Still hard-fails on malformed entries (no silent drop / Flat default).
    """
    sampling = dict(config.get("Sampling") or {})
    raw_vars = sampling.get("Variables") or []
    if not isinstance(raw_vars, list):
        raise ValueError("Sampling.Variables must be a list")
    variables: list[Variable] = []
    allowed = frozenset(
        {
            "Flat", "Log", "Logit", "Normal", "Log-Normal", "Binomial", "Poisson",
            "Beta", "Exponential", "Gamma",
        }
    )
    for index, item in enumerate(raw_vars):
        if not isinstance(item, dict):
            raise ValueError(
                f"Sampling.Variables[{index}] must be a mapping, got {type(item).__name__}"
            )
        name = str(item.get("name", "")).strip()
        if not name:
            raise ValueError(
                f"Sampling.Variables[{index}].name is required "
                "(nameless entries are not silently dropped)"
            )
        dist_block = item.get("distribution")
        if not isinstance(dist_block, dict):
            raise ValueError(
                f"Sampling.Variables[{index}].distribution is required and must be a mapping"
            )
        dist_type = str(dist_block.get("type", "")).strip()
        if dist_type not in allowed:
            raise ValueError(
                f"Sampling.Variables[{index}].distribution.type {dist_type!r} "
                f"is not supported; expected one of: {', '.join(sorted(allowed))}"
            )
        parameters = dist_block.get("parameters")
        if not isinstance(parameters, dict):
            raise ValueError(
                f"Sampling.Variables[{index}].distribution.parameters is required"
            )
        variables.append(
            Variable(
                name=name,
                description=str(item.get("description", "")),
                distribution=dist_type,
                parameters=dict(parameters),
            )
        )
    if not variables:
        raise ValueError("Sampling.Variables must define at least one variable")
    return variables


def _numeric_string_to_float(value: object) -> object:
    if not isinstance(value, str):
        return value
    try:
        return float(value)
    except ValueError:
        return value


__all__ = ["Variable", "load_variables"]
