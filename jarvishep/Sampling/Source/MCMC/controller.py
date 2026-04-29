#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Protocol, Sequence, Tuple


@dataclass
class MCMCControlPatch:
    temperature_ladder: Optional[Sequence[float]] = None
    exchange_interval: Optional[int] = None
    exchange_rule: Optional[str] = None
    swap_pairing_mode: Optional[str] = None
    proposal_scales: Optional[Sequence[float] | float] = None
    chain_priority_weight: Optional[Sequence[float]] = None


class MCMCControllerProtocol(Protocol):
    def on_run_start(self, context: Dict[str, Any]) -> Optional[MCMCControlPatch]:
        ...

    def on_pre_step(self, snapshot: Dict[str, Any]) -> Optional[MCMCControlPatch]:
        ...

    def on_post_step(
        self, snapshot: Dict[str, Any], outcome: Dict[str, Any]
    ) -> Optional[MCMCControlPatch]:
        ...

    def on_pre_exchange(self, snapshot: Dict[str, Any]) -> Optional[MCMCControlPatch]:
        ...

    def on_post_exchange(
        self, snapshot: Dict[str, Any], exchange_metrics: Dict[str, Any]
    ) -> Optional[MCMCControlPatch]:
        ...


class NoopMCMCController:
    def on_run_start(self, context: Dict[str, Any]) -> Optional[MCMCControlPatch]:
        return None

    def on_pre_step(self, snapshot: Dict[str, Any]) -> Optional[MCMCControlPatch]:
        return None

    def on_post_step(
        self, snapshot: Dict[str, Any], outcome: Dict[str, Any]
    ) -> Optional[MCMCControlPatch]:
        return None

    def on_pre_exchange(self, snapshot: Dict[str, Any]) -> Optional[MCMCControlPatch]:
        return None

    def on_post_exchange(
        self, snapshot: Dict[str, Any], exchange_metrics: Dict[str, Any]
    ) -> Optional[MCMCControlPatch]:
        return None


class MCMCControlGuard:
    @staticmethod
    def validate_patch(patch: MCMCControlPatch, nchains: int) -> Tuple[bool, str]:
        if patch is None:
            return True, "empty patch"

        nchains = int(max(1, nchains))

        if patch.temperature_ladder is not None:
            ladder = [float(x) for x in patch.temperature_ladder]
            if len(ladder) != nchains:
                return False, "temperature_ladder size mismatch"
            if any(x <= 0.0 for x in ladder):
                return False, "temperature_ladder must be positive"
            for i in range(1, len(ladder)):
                if ladder[i] < ladder[i - 1]:
                    return False, "temperature_ladder must be non-decreasing"

        if patch.exchange_interval is not None:
            if int(patch.exchange_interval) <= 0:
                return False, "exchange_interval must be > 0"

        if patch.proposal_scales is not None:
            scales = patch.proposal_scales
            if isinstance(scales, (int, float)):
                if float(scales) <= 0.0:
                    return False, "proposal_scales scalar must be > 0"
            else:
                vals = [float(x) for x in scales]
                if len(vals) != nchains:
                    return False, "proposal_scales size mismatch"
                if any(v <= 0.0 for v in vals):
                    return False, "proposal_scales values must be > 0"

        if patch.chain_priority_weight is not None:
            weights = [float(x) for x in patch.chain_priority_weight]
            if len(weights) != nchains:
                return False, "chain_priority_weight size mismatch"
            if any(w <= 0.0 for w in weights):
                return False, "chain_priority_weight values must be > 0"

        return True, "ok"

