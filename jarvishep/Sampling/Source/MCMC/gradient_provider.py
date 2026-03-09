#!/usr/bin/env python3
from __future__ import annotations

from typing import Any, Mapping, Protocol, runtime_checkable

import numpy as np


@runtime_checkable
class GradientProviderProtocol(Protocol):
    """Protocol for gradient-enabled MCMC engines (future strict path)."""

    def gradient(
        self,
        point: np.ndarray,
        *,
        context: Mapping[str, Any] | None = None,
    ) -> np.ndarray: ...


def validate_gradient_provider(provider: Any, *, method_name: str) -> None:
    if provider is None:
        return
    gradient_fn = getattr(provider, "gradient", None)
    if not callable(gradient_fn):
        raise TypeError(
            f"{method_name} gradient provider must expose callable 'gradient(point, *, context=None)'."
        )
