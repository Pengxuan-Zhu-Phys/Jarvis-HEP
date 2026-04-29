#!/usr/bin/env python3
from __future__ import annotations

import inspect
from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class MCMCEngineProtocol(Protocol):
    """Runtime contract for ChainRuntime.engine objects."""

    param: Any
    proposed_param: Any
    proposal_scale: float
    n_iterations: int
    iterations: int
    last_loglikelihood: float | None

    def __next__(self): ...

    def update(self, new_loglikelihood: float, beta: float = 1.0) -> bool: ...


class MCMCEngineContractError(TypeError):
    """Raised when a chain engine does not satisfy the required contract."""


def _validate_update_signature(update_method: Any) -> bool:
    try:
        sig = inspect.signature(update_method)
    except (TypeError, ValueError):
        return True

    params = list(sig.parameters.values())
    positional = [
        p
        for p in params
        if p.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)
    ]
    var_positional = any(p.kind == inspect.Parameter.VAR_POSITIONAL for p in params)

    # Bound methods should accept at least one positional argument: new_loglikelihood.
    return bool(len(positional) >= 1 or var_positional)


def validate_engine_contract(
    engine: Any,
    *,
    sampler_method: str,
    chain_id: int,
) -> None:
    if engine is None:
        raise MCMCEngineContractError(
            f"{sampler_method} chain[{chain_id}] engine is None; expected MCMC engine object."
        )

    required_attrs = (
        "param",
        "proposed_param",
        "proposal_scale",
        "n_iterations",
        "iterations",
        "last_loglikelihood",
    )
    missing = [name for name in required_attrs if not hasattr(engine, name)]
    if missing:
        raise MCMCEngineContractError(
            f"{sampler_method} chain[{chain_id}] engine missing required fields: {missing}"
        )

    next_fn = getattr(engine, "__next__", None)
    if not callable(next_fn):
        raise MCMCEngineContractError(
            f"{sampler_method} chain[{chain_id}] engine missing callable '__next__'."
        )

    update_fn = getattr(engine, "update", None)
    if not callable(update_fn):
        raise MCMCEngineContractError(
            f"{sampler_method} chain[{chain_id}] engine missing callable 'update(...)'."
        )
    if not _validate_update_signature(update_fn):
        raise MCMCEngineContractError(
            f"{sampler_method} chain[{chain_id}] engine.update signature is invalid; "
            "expected at least one positional argument for new_loglikelihood."
        )
