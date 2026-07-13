#!/usr/bin/env python3
"""V1-compatible, dependency-light functions for the V2 expression runtime."""

from __future__ import annotations

from functools import reduce
from typing import Any

import numpy as np
import sympy as sp


def _scalar_or_array(value: Any) -> Any:
    array = np.asarray(value)
    if array.ndim == 0:
        return array.item()
    return array


def Gauss(xx: Any, mean: Any, err: Any) -> Any:
    """V1 unnormalised Gaussian: ``exp(-0.5 * ((x-mean)/err)**2)``."""
    value = np.exp(-0.5 * ((np.asarray(xx) - mean) / err) ** 2)
    return _scalar_or_array(value)


def Normal(xx: Any, mean: Any, err: Any) -> Any:
    """V1 normalised Gaussian probability density."""
    value = Gauss(xx, mean, err) / (np.asarray(err) * np.sqrt(2.0 * np.pi))
    return _scalar_or_array(value)


def LogGauss(xx: Any, mean: Any, err: Any) -> Any:
    """V1 Gaussian log-kernel without the normalisation constant."""
    value = -0.5 * ((np.asarray(xx) - mean) / err) ** 2
    return _scalar_or_array(value)


def Heaviside(value: Any, zero_value: Any = 0.5) -> Any:
    """SymPy-compatible Heaviside with the V1 default ``H(0) = 1/2``."""
    result = np.heaviside(np.asarray(value), np.asarray(zero_value))
    return _scalar_or_array(result)


def _sec(value: Any) -> Any:
    return _scalar_or_array(1.0 / np.cos(value))


def _csc(value: Any) -> Any:
    return _scalar_or_array(1.0 / np.sin(value))


def _cot(value: Any) -> Any:
    return _scalar_or_array(1.0 / np.tan(value))


def _asec(value: Any) -> Any:
    return _scalar_or_array(np.arccos(1.0 / np.asarray(value)))


def _acsc(value: Any) -> Any:
    return _scalar_or_array(np.arcsin(1.0 / np.asarray(value)))


def _acot(value: Any) -> Any:
    return _scalar_or_array(np.arctan(1.0 / np.asarray(value)))


def _sech(value: Any) -> Any:
    return _scalar_or_array(1.0 / np.cosh(value))


def _csch(value: Any) -> Any:
    return _scalar_or_array(1.0 / np.sinh(value))


def _coth(value: Any) -> Any:
    return _scalar_or_array(1.0 / np.tanh(value))


def _asech(value: Any) -> Any:
    return _scalar_or_array(np.arccosh(1.0 / np.asarray(value)))


def _acsch(value: Any) -> Any:
    return _scalar_or_array(np.arcsinh(1.0 / np.asarray(value)))


def _acoth(value: Any) -> Any:
    return _scalar_or_array(np.arctanh(1.0 / np.asarray(value)))


def _minimum(*values: Any) -> Any:
    if not values:
        raise TypeError("Min requires at least one argument")
    return _scalar_or_array(reduce(np.minimum, values))


def _maximum(*values: Any) -> Any:
    if not values:
        raise TypeError("Max requires at least one argument")
    return _scalar_or_array(reduce(np.maximum, values))


def _root(value: Any, degree: Any, branch: Any = 0) -> Any:
    principal = np.emath.power(value, 1.0 / np.asarray(degree))
    result = principal * np.exp(2j * np.pi * np.asarray(branch) / np.asarray(degree))
    result = np.real_if_close(result)
    return _scalar_or_array(result)


# This is the frozen V1 ``inner_func._Inner_FCs`` surface plus the four V1 helper
# functions. Values are efficient numerical implementations used by lambdify.
NUMERIC_MODULES: dict[str, Any] = {
    "log": np.log,
    "exp": np.exp,
    "ln": np.log,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "sec": _sec,
    "csc": _csc,
    "cot": _cot,
    # SymPy's NumPy printer emits ``sinc(x / pi)``; numpy.sinc is the
    # normalized sin(pi*x)/(pi*x), yielding the intended SymPy sin(x)/x.
    "sinc": np.sinc,
    "asin": np.arcsin,
    "acos": np.arccos,
    "atan": np.arctan,
    "asec": _asec,
    "acsc": _acsc,
    "acot": _acot,
    "atan2": np.arctan2,
    "sinh": np.sinh,
    "cosh": np.cosh,
    "tanh": np.tanh,
    "sech": _sech,
    "csch": _csch,
    "coth": _coth,
    "asinh": np.arcsinh,
    "acosh": np.arccosh,
    "atanh": np.arctanh,
    "acoth": _acoth,
    "asech": _asech,
    "acsch": _acsch,
    "sqrt": np.sqrt,
    "Min": _minimum,
    "Max": _maximum,
    "root": _root,
    "Abs": np.abs,
    "Gauss": Gauss,
    "Normal": Normal,
    "LogGauss": LogGauss,
    "Heaviside": Heaviside,
}

EXPRESSION_FUNCTION_NAMES = frozenset(NUMERIC_MODULES)


EXPRESSION_CONSTANTS: dict[str, Any] = {
    "Pi": sp.pi,
    "pi": sp.pi,
    "PI": sp.pi,
    "E": sp.E,
    "Inf": sp.oo,
}


__all__ = [
    "EXPRESSION_CONSTANTS",
    "EXPRESSION_FUNCTION_NAMES",
    "Gauss",
    "Heaviside",
    "LogGauss",
    "NUMERIC_MODULES",
    "Normal",
]
