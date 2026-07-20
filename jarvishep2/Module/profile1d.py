#!/usr/bin/env python3
"""Profile1D golden-section nuisance profiler (D13.4 / V1 Profile1D science)."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from math import sqrt
from typing import Any

from jarvishep2.Module.nuisance import (
    NuisanceExpressionRegistry,
    NuisancePassConditionRegistry,
    extract_nuisance_config,
    parse_nuisance_variable,
)
from jarvishep2.expression import ExpressionContext
from jarvishep2.logging import get_jarvis_logger

# Golden-section constants (V1 Profile1D).
_PHI = 0.5 + sqrt(5.0) / 2.0
_INVPHI = 1.0 / _PHI
_INVPHI2 = _INVPHI * _INVPHI


@dataclass
class ProfileProbeResult:
    """One probe evaluation at a nuisance value."""

    z: float
    logl: float
    terms: dict[str, float] = field(default_factory=dict)
    pass_ok: bool = True
    pass_terms: dict[str, bool] = field(default_factory=dict)
    observables: dict[str, Any] = field(default_factory=dict)


@dataclass
class Profile1DResult:
    """Final profile outcome applied back onto the Sample."""

    var_name: str
    best_z: float
    best_logl: float
    best_terms: dict[str, float]
    pass_ok: bool
    pass_terms: dict[str, bool]
    n_attempts: int
    history_z: list[float]
    history_logl: list[float]
    mode: str
    status: str  # Accept | Failed


EvaluateFn = Callable[[float], ProfileProbeResult]


class Profile1DProfiler:
    """Golden-section search over one nuisance variable (V1 Profile1D).

    ``evaluate(z)`` is injected by the Worker: it may re-run calc/opera with the
    free nuisance parameter set, then evaluate expression LogL + pass conditions
    (compile-once registries — no recompile in the inner loop).
    """

    def __init__(
        self,
        *,
        var_name: str,
        zmin: float,
        zmax: float,
        logl_registry: NuisanceExpressionRegistry,
        pass_registry: NuisancePassConditionRegistry,
        mode: str = "min",
        max_attempt: int = 50,
        # D13.7c: default True matches V1 Profile1D, which re-ran the full
        # sample pipeline per NAttempt probe (not expression-only).
        re_run_physics: bool = True,
    ) -> None:
        self.var_name = str(var_name)
        self.zmin = float(zmin)
        self.zmax = float(zmax)
        self.logl_registry = logl_registry
        self.pass_registry = pass_registry
        self.mode = str(mode or "min").lower()
        self.max_attempt = max(1, int(max_attempt))
        self.re_run_physics = bool(re_run_physics)
        self._logger = get_jarvis_logger("nuisance.profile1d")

    @classmethod
    def from_config(
        cls,
        config: Mapping[str, Any],
        *,
        context: ExpressionContext | None = None,
    ) -> Profile1DProfiler | None:
        block = extract_nuisance_config(config)
        if block is None:
            return None
        method = str(block.get("Method") or block.get("method") or "Profile1D").strip()
        if method and method not in ("Profile1D", "profile1d", "profile_1d"):
            # Only Profile1D ships in D13.4; unknown methods still try Profile1D shape.
            pass
        ctx = context or ExpressionContext()
        logl_reg = NuisanceExpressionRegistry(ctx)
        pass_reg = NuisancePassConditionRegistry(ctx)
        logl_reg.load_from_config(block.get("LogLikelihood") or block.get("loglikelihood"))
        pass_reg.load_from_config(block.get("PassCondition") or block.get("pass_condition"))
        name, zmin, zmax = parse_nuisance_variable(block)
        mode = str(block.get("TargetMode") or block.get("target_mode") or "min")
        max_attempt = int(block.get("MaxAttempt") or block.get("max_attempt") or 50)
        # Default True: V1 re-executed each probe as a full sample (NAttempt
        # card). Set re_run_physics/rerun_physics: false for pure expression
        # nuisances that do not feed calculators/Operas.
        re_run = block.get("re_run_physics", block.get("rerun_physics", True))
        return cls(
            var_name=name,
            zmin=zmin,
            zmax=zmax,
            logl_registry=logl_reg,
            pass_registry=pass_reg,
            mode=mode,
            max_attempt=max_attempt,
            re_run_physics=bool(re_run),
        )

    def _prefer(self, a: float, b: float) -> bool:
        """Return True if score ``a`` is better than ``b`` under TargetMode."""
        if self.mode in ("max", "maximize"):
            return a > b
        return a < b

    def optimize(self, evaluate: EvaluateFn) -> Profile1DResult:
        """Run golden-section search; ``evaluate(z)`` probes one nuisance value."""
        a, b = self.zmin, self.zmax
        c = b - (b - a) * _INVPHI
        d = a + (b - a) * _INVPHI
        history_z: list[float] = []
        history_logl: list[float] = []
        best: ProfileProbeResult | None = None

        def _consider(result: ProfileProbeResult) -> ProfileProbeResult:
            nonlocal best
            if best is None:
                best = result
            elif self._prefer(result.logl, best.logl):
                best = result
            elif result.logl == best.logl and result.pass_ok and not best.pass_ok:
                best = result
            return result

        def probe(z: float) -> ProfileProbeResult:
            result = evaluate(float(z))
            history_z.append(float(z))
            history_logl.append(float(result.logl))
            return _consider(result)

        fc = probe(c)
        fd = probe(d)
        attempts = 2
        tol = 1e-12 * max(1.0, abs(a), abs(b))

        while attempts < self.max_attempt and abs(b - a) > tol:
            # V1: gn = -fn for TargetMode=max, else gn = fn; if gn[c] < gn[d] keep [a,d].
            if self.mode in ("max", "maximize"):
                gn_c, gn_d = -fc.logl, -fd.logl
            else:
                gn_c, gn_d = fc.logl, fd.logl
            if gn_c < gn_d:
                b, d, fd = d, c, fc
                c = b - (b - a) * _INVPHI
                fc = probe(c)
            else:
                a, c, fc = c, d, fd
                d = a + (b - a) * _INVPHI
                fd = probe(d)
            attempts += 1

        assert best is not None
        status = "Accept" if best.pass_ok else "Failed"
        self._logger.info(
            "Profile1D done var=%s best_z=%.6g best_logl=%.6g pass=%s attempts=%d",
            self.var_name,
            best.z,
            best.logl,
            best.pass_ok,
            attempts,
        )
        return Profile1DResult(
            var_name=self.var_name,
            best_z=float(best.z),
            best_logl=float(best.logl),
            best_terms=dict(best.terms),
            pass_ok=bool(best.pass_ok),
            pass_terms=dict(best.pass_terms),
            n_attempts=int(attempts),
            history_z=list(history_z),
            history_logl=list(history_logl),
            mode=self.mode,
            status=status,
        )


__all__ = [
    "Profile1DProfiler",
    "Profile1DResult",
    "ProfileProbeResult",
]
