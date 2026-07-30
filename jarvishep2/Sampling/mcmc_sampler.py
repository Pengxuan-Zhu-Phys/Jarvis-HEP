#!/usr/bin/env python3
"""Feedback-driven MCMC family samplers (D13.2 + D13.3).

Ports V1 chain science (``Sampling/Source/MCMC/*`` engines) onto
:class:`FeedbackSampler`: one generation = one Metropolis step across all
chains (DRAM delayed-rejection stages are intra-generation mini-batches;
ensemble methods use half-ensemble barriers; PT swaps at barriers).
Never imports ``jarvishep``.
"""

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from jarvishep2.Sampling.Source.MCMC.ammcmc_chain import AMMCMCChain
from jarvishep2.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep2.Sampling.Source.MCMC.config_contract import (
    bounds_get,
    bounds_get_bool,
    bounds_get_float,
    bounds_get_int,
    bounds_get_list,
    normalize_proposal_scales,
    parse_common_chain_counts,
    parse_proposal_scale_value,
)
from jarvishep2.Sampling.Source.MCMC.dram_chain import DRAMChain
from jarvishep2.Sampling.Source.MCMC.engine_demcmc import DEMCMCChain
from jarvishep2.Sampling.Source.MCMC.engine_ensemble import EnsembleChain
from jarvishep2.Sampling.Source.MCMC.mcmc_chain import MCMCChain
from jarvishep2.Sampling.feedback_sampler import FeedbackSampler
from jarvishep2.Sampling.sampling_utils import evaluate_selection, map_u_to_physical
from jarvishep2.Sampling.variables import load_variables
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample

_FAILED_LOGL = -1.0e300

# Methods that use complementary-half ensemble proposals (emcee-style).
_ENSEMBLE_METHODS = frozenset(
    {"EnsembleMCMC", "Ensemble", "DEMCMC", "PTEnsemble"}
)
# Methods that run temperature-exchange after a full step.
_PT_METHODS = frozenset({"PTMCMC", "PT", "PTEnsemble"})


def mcmc_sample_uuid(
    *,
    method: str,
    seed: int,
    chain_id: int,
    step: int,
    stage: int,
) -> str:
    """Deterministic uuid from (method, seed, chain, step, stage) — D13 §3.3."""
    digest = hashlib.sha256(
        f"{method}:{int(seed)}:{int(chain_id)}:{int(step)}:{int(stage)}".encode("utf-8")
    ).hexdigest()
    return (
        f"{digest[0:8]}-{digest[8:12]}-{digest[12:16]}-{digest[16:20]}-{digest[20:32]}"
    )


class MCMCSampler(FeedbackSampler):
    """Metropolis / ensemble / PT family on the V2 redis/feedback runtime.

    YAML ``Sampling.Method``: ``MCMC`` | ``AMMCMC`` | ``AM`` | ``DRAM`` |
    ``EnsembleMCMC`` | ``Ensemble`` | ``DEMCMC`` | ``PTMCMC`` | ``PT`` |
    ``PTEnsemble``.
    """

    method = "MCMC"

    def __init__(self, method_name: str = "MCMC") -> None:
        super().__init__()
        self.method = str(method_name)
        self._logger = get_jarvis_logger(f"sampler.{self.method.lower()}")
        self.vars: list = []
        self._dim = 0
        self._nchains = 1
        self._niters = 1
        self._proposal_scales: float | list[float] = 0.1
        self._selectionexp: str | None = None
        self._registry: ChainRegistry | None = None
        self._uuid_to_meta: dict[str, dict[str, int]] = {}
        self._finished = False
        self._total_accepted = 0
        self._total_proposed = 0
        self._failed_uuids: list[str] = []
        # AM / DRAM extras
        self._adapt_enabled = True
        self._adapt_start_iter = 100
        self._adapt_window = 25
        self._adapt_eps = 1e-6
        self._adapt_scale = 2.38
        self._dr_steps = 2
        self._dr_scale_factors = [1.0, 0.5]
        # Ensemble / DE (D13.3)
        self._stretch_a = 2.0
        self._de_gamma = 0.0
        self._de_noise = 1.0e-3
        self._de_crossover = 1.0
        # PT (D13.3)
        self._temperature_ladder: list[float] = [1.0]
        self._exchange_interval = 1
        self._exchange_offset = 0
        self._swap_attempts = 0
        self._swap_accepts = 0
        self._exchange_rng: np.random.Generator | None = None
        self._summary: dict[str, Any] | None = None

    # ------------------------------------------------------------------ config
    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._dim = len(self.vars)
        if self._dim < 1:
            raise ValueError(f"{self.method} requires at least one Sampling.Variable")
        bounds = sampling.get("Bounds") or {}
        if not isinstance(bounds, Mapping):
            raise ValueError("Sampling.Bounds must be a mapping for MCMC methods")
        self._nchains, self._niters = parse_common_chain_counts(bounds)
        self._proposal_scales = parse_proposal_scale_value(bounds, default=0.1)
        self._selectionexp = sampling.get("selection")
        seed = int(sampling.get("Seed", sampling.get("seed", 0)) or 0)
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        self._init_seed_sequence(seed)
        on_fail = bounds.get("on_failure", sampling.get("on_failure", "reject"))
        self._on_failure = str(on_fail or "reject").lower()

        if self.method in ("AMMCMC", "AM", "DRAM"):
            self._adapt_enabled = bounds_get_bool(
                bounds, "adapt_enabled", aliases=("adapt.enabled",), default=True
            )
            self._adapt_start_iter = bounds_get_int(
                bounds,
                "adapt_start_iter",
                aliases=("adapt.start_iter",),
                default=100,
                minimum=1,
            )
            self._adapt_window = bounds_get_int(
                bounds,
                "adapt_window",
                aliases=("adapt.window",),
                default=25,
                minimum=1,
            )
            self._adapt_eps = bounds_get_float(
                bounds,
                "adapt_eps",
                aliases=("adapt.eps",),
                default=1e-6,
                minimum=0.0,
            )
            self._adapt_scale = bounds_get_float(
                bounds,
                "adapt_scale",
                aliases=("adapt.scale",),
                default=2.38,
                minimum=1e-9,
            )
        if self.method == "DRAM":
            self._dr_steps = bounds_get_int(
                bounds,
                "dr_steps",
                aliases=("dr.steps",),
                default=2,
                minimum=1,
            )
            factors = bounds_get_list(
                bounds,
                "dr_scale_factors",
                aliases=("dr.scale_factors",),
                default=[1.0, 0.5],
            )
            if isinstance(factors, (int, float)):
                factors = [1.0, float(factors)]
            self._dr_scale_factors = [max(1e-6, float(x)) for x in factors]

        if self.method in ("EnsembleMCMC", "Ensemble", "PTEnsemble"):
            self._stretch_a = bounds_get_float(
                bounds,
                "stretch_a",
                aliases=("ensemble.stretch_a",),
                default=2.0,
                minimum=1.01,
            )
        if self.method == "DEMCMC":
            self._de_gamma = bounds_get_float(
                bounds,
                "de_gamma",
                aliases=("de.gamma",),
                default=0.0,
                minimum=0.0,
            )
            self._de_noise = bounds_get_float(
                bounds,
                "de_noise",
                aliases=("de.noise",),
                default=1.0e-3,
                minimum=0.0,
            )
            self._de_crossover = bounds_get_float(
                bounds,
                "de_crossover",
                aliases=("de.crossover",),
                default=1.0,
                minimum=0.0,
            )
        if self.method in _PT_METHODS:
            self._exchange_interval = bounds_get_int(
                bounds,
                "exchange_interval",
                aliases=("exchange.interval",),
                default=1,
                minimum=1,
            )
            ladder = bounds_get(
                bounds,
                "temperature_ladder",
                aliases=("temperature.ladder",),
                default=None,
            )
            if ladder is None:
                ladder = [1.0 + float(ii) for ii in range(self._nchains)]
            self._temperature_ladder = [float(x) for x in ladder]
            if len(self._temperature_ladder) != self._nchains:
                raise ValueError(
                    f"temperature_ladder size mismatch for {self.method}: "
                    f"expect {self._nchains}, got {len(self._temperature_ladder)}"
                )
        if self.method in _ENSEMBLE_METHODS and self._nchains < 2:
            raise ValueError(f"{self.method} requires num_chains >= 2")

        if self._nchains < workers:
            self._logger.warning(
                "%s: chains (%d) < workers (%d) — workers may idle at barriers; "
                "prefer chains >= workers",
                self.method,
                self._nchains,
                workers,
            )

    def _normalize_scales(self) -> list[float]:
        return normalize_proposal_scales(
            self._proposal_scales,
            nchains=self._nchains,
            sampler_method=self.method,
        )

    def _uses_half_ensemble(self) -> bool:
        return self.method in _ENSEMBLE_METHODS

    def _uses_pt(self) -> bool:
        return self.method in _PT_METHODS

    def _population_for(self, chain_id: int) -> list[np.ndarray]:
        """Complementary population for ensemble / DE moves.

        Half-ensemble: partners are the other half of walkers (emcee semantics).
        """
        assert self._registry is not None
        ids = self._registry.ids()
        mid = max(1, len(ids) // 2)
        if int(chain_id) < mid:
            partner_ids = [j for j in ids if j >= mid]
        else:
            partner_ids = [j for j in ids if j < mid]
        out: list[np.ndarray] = []
        for j in partner_ids:
            eng = self._registry.get(j).engine
            if eng is None:
                continue
            out.append(np.asarray(eng.param, dtype=float))
        return out

    def _make_engine(self, chain_id: int, scale: float, rng: np.random.Generator):
        initial = rng.random(self._dim)
        if self.method == "DRAM":
            return DRAMChain(
                initial,
                scale,
                self._niters,
                adapt_enabled=self._adapt_enabled,
                adapt_start_iter=self._adapt_start_iter,
                adapt_window=self._adapt_window,
                adapt_eps=self._adapt_eps,
                adapt_scale=self._adapt_scale,
                dr_steps=self._dr_steps,
                dr_scale_factors=self._dr_scale_factors,
                rng=rng,
            )
        if self.method in ("AMMCMC", "AM"):
            return AMMCMCChain(
                initial,
                scale,
                self._niters,
                adapt_enabled=self._adapt_enabled,
                adapt_start_iter=self._adapt_start_iter,
                adapt_window=self._adapt_window,
                adapt_eps=self._adapt_eps,
                adapt_scale=self._adapt_scale,
                rng=rng,
            )
        if self.method in ("EnsembleMCMC", "Ensemble", "PTEnsemble"):
            return EnsembleChain(
                initial,
                scale,
                self._niters,
                chain_id=chain_id,
                population_getter=self._population_for,
                stretch_a=self._stretch_a,
                rng=rng,
            )
        if self.method == "DEMCMC":
            return DEMCMCChain(
                initial,
                scale,
                self._niters,
                chain_id=chain_id,
                population_getter=self._population_for,
                de_gamma=self._de_gamma,
                de_noise=self._de_noise,
                de_crossover=self._de_crossover,
                rng=rng,
            )
        return MCMCChain(initial, scale, self._niters, rng=rng)

    def _ensure_registry(self) -> ChainRegistry:
        if self._registry is not None:
            return self._registry
        root = self._ensure_seed_sequence()
        # One child SeedSequence per chain → worker-count-independent trajectories.
        n_spawn = self._nchains + (1 if self._uses_pt() else 0)
        children = root.spawn(max(1, n_spawn))
        scales = self._normalize_scales()
        ladder = (
            list(self._temperature_ladder)
            if self._uses_pt()
            else [1.0] * self._nchains
        )
        cold_temp = min(ladder) if ladder else 1.0
        cold_idx = 0
        for ii, temp in enumerate(ladder):
            if float(temp) == float(cold_temp):
                cold_idx = ii
                break
        chains: list[ChainRuntime] = []
        for cid in range(self._nchains):
            rng = np.random.default_rng(children[cid])
            engine = self._make_engine(cid, scales[cid], rng)
            chains.append(
                ChainRuntime(
                    chain_id=cid,
                    engine=engine,
                    temperature=float(ladder[cid] if cid < len(ladder) else 1.0),
                    is_cold=(cid == cold_idx),
                    open_stage=None,
                )
            )
        self._registry = ChainRegistry(chains)
        if self._uses_pt() and self._exchange_rng is None:
            self._exchange_rng = np.random.default_rng(children[self._nchains])
        return self._registry

    # --------------------------------------------------------------- proposals
    def _propose_u_for_chain(self, chain: ChainRuntime, stage: int) -> np.ndarray:
        engine = chain.engine
        if hasattr(engine, "propose_stage"):
            u = np.asarray(engine.propose_stage(stage), dtype=np.float64)
        else:
            if stage != 0:
                raise RuntimeError(f"{self.method} does not support stage={stage}")
            u = np.asarray(next(engine), dtype=np.float64)
        # Selection filter: redraw (selection-only; does not advance iterations).
        if self._selectionexp:
            for _ in range(4096):
                physical = map_u_to_physical(u, self.vars)
                if evaluate_selection(
                    self._selectionexp,
                    physical,
                    context=self._expression_context,
                ):
                    break
                if hasattr(engine, "propose_stage"):
                    u = np.asarray(engine.propose_stage(stage), dtype=np.float64)
                else:
                    # re-propose without advancing: undo is hard for plain MCMC;
                    # re-call next would skip — for selection we redraw from same state
                    # by manually sampling like engine.__next__ without bumping iters.
                    # Safer: only re-run propose_stage path for DRAM/AM; for MCMC use
                    # a local redraw of the same engine next after temporarily
                    # decrementing — not available. Redraw gaussian from current param.
                    if engine.iterations == 0:
                        u = engine._rng.random(self._dim)
                    else:
                        step = engine._rng.normal(
                            0.0, engine.proposal_scale, size=self._dim
                        )
                        u = np.clip(engine.param + step, 0.0, 1.0)
                    engine.proposed_param = u
            else:
                self._logger.warning(
                    "chain %s stage %s: selection never passed; using last proposal",
                    chain.chain_id,
                    stage,
                )
        return u

    def _build_chain_sample(
        self, chain: ChainRuntime, stage: int, u: np.ndarray
    ) -> Sample:
        step = int(chain.engine.iterations)
        uuid = mcmc_sample_uuid(
            method=self.method,
            seed=self._seed,
            chain_id=chain.chain_id,
            step=step,
            stage=stage,
        )
        sample = self._build_sample(u)
        sample.uuid = uuid
        # Additive diagnostic columns for DATABASE (design goal 6).
        sample.observables = dict(sample.observables or {})
        sample.observables.setdefault("chain_id", int(chain.chain_id))
        sample.observables.setdefault("step", int(step))
        sample.observables.setdefault("stage", int(stage))
        self._uuid_to_meta[uuid] = {
            "chain_id": int(chain.chain_id),
            "step": int(step),
            "stage": int(stage),
        }
        return sample

    def _half_chain_ids(self, half: int) -> list[int]:
        ids = self._ensure_registry().ids()
        mid = max(1, len(ids) // 2)
        if int(half) == 0:
            return [i for i in ids if i < mid]
        return [i for i in ids if i >= mid]

    def _propose_for_stages(
        self,
        *,
        stages: str,
        chain_ids: Sequence[int] | None = None,
    ) -> list[Sample]:
        """Build proposals.

        ``stages``:
          * ``"step"`` — stage-0 for every unfinished chain whose open_stage is
            idle (``None``); starts a new Metropolis step for that chain.
          * ``"followup"`` — only delayed-rejection stages (``open_stage > 0``).

        ``chain_ids``: optional subset (half-ensemble); default = all chains.
        """
        registry = self._ensure_registry()
        allowed = None if chain_ids is None else {int(c) for c in chain_ids}
        samples: list[Sample] = []
        for chain in registry.all():
            if allowed is not None and int(chain.chain_id) not in allowed:
                continue
            engine = chain.engine
            if engine.iterations >= engine.n_iterations:
                chain.open_stage = None
                continue
            if stages == "followup":
                if chain.open_stage is None or int(chain.open_stage) <= 0:
                    continue
                stage = int(chain.open_stage)
            else:
                # New step: skip chains still mid delayed-rejection.
                if chain.open_stage is not None and int(chain.open_stage) > 0:
                    continue
                stage = 0
                chain.open_stage = 0
            u = self._propose_u_for_chain(chain, stage)
            samples.append(self._build_chain_sample(chain, stage, u))
        return samples

    def propose_generation(self) -> Sequence[Sample] | None:
        """Propose stage-0 for the next Metropolis step across active chains."""
        samples = self._propose_for_stages(stages="step")
        if not samples:
            return None
        return samples

    # ----------------------------------------------------------------- absorb
    def _extract_logl(self, record: Mapping[str, Any]) -> float | None:
        """Read flat ``logL`` from feedback; −∞ / non-finite / missing → failed."""
        from jarvishep2.feedback_return import (
            WIRE_LOGL_KEY,
            feedback_logl,
            is_unusable_logl,
            normalize_feedback_record,
        )

        flat = normalize_feedback_record(record)
        if WIRE_LOGL_KEY not in flat:
            return None
        logl = feedback_logl(flat)
        if is_unusable_logl(logl):
            return None
        return float(logl)

    def absorb_generation(self, results: Sequence[Mapping[str, Any]]) -> None:
        registry = self._ensure_registry()
        # Deterministic absorb order (not worker arrival order).
        ordered = sorted(
            results,
            key=lambda r: (
                self._uuid_to_meta.get(str(r.get("uuid", "")), {}).get("chain_id", 10**9),
                self._uuid_to_meta.get(str(r.get("uuid", "")), {}).get("stage", 0),
            ),
        )
        for record in ordered:
            uuid = str(record.get("uuid", ""))
            meta = self._uuid_to_meta.get(uuid)
            if meta is None:
                continue
            chain = registry.get(int(meta["chain_id"]))
            stage = int(meta["stage"])
            logl = self._extract_logl(record)
            if logl is None:
                self._failed_uuids.append(uuid)
                if self._failure_policy_halt(record):
                    raise RuntimeError(
                        f"{self.method}: Failed sample {uuid} with on_failure=halt"
                    )
                logl = _FAILED_LOGL
            beta = 1.0
            if chain.temperature > 0:
                beta = 1.0 / float(chain.temperature)
            engine = chain.engine
            self._total_proposed += 1
            if hasattr(engine, "consume_stage_result"):
                outcome = engine.consume_stage_result(stage, logl, beta=beta)
                accepted = bool(outcome.get("accepted", False))
                if accepted:
                    self._total_accepted += 1
                    chain.accepted += 1
                else:
                    chain.rejected += 1
                if outcome.get("iteration_done"):
                    chain.iter = int(engine.iterations)
                    chain.last_logl = getattr(engine, "last_loglikelihood", logl)
                    # Idle until the outer loop starts the next step.
                    chain.open_stage = None
                    chain.window_iter = int(chain.window_iter) + 1
                    chain.history.append_from_values(
                        iter=chain.iter,
                        state="UPDATE",
                        proposal=getattr(engine, "proposed_param", None),
                        logl=chain.last_logl,
                        accepted=accepted,
                        temperature=chain.temperature,
                    )
                else:
                    next_stage = outcome.get("next_stage")
                    chain.open_stage = (
                        int(next_stage) if next_stage is not None else None
                    )
            else:
                accepted = bool(engine.update(logl, beta=beta))
                if accepted:
                    self._total_accepted += 1
                    chain.accepted += 1
                else:
                    chain.rejected += 1
                chain.iter = int(engine.iterations)
                chain.last_logl = getattr(engine, "last_loglikelihood", logl)
                chain.open_stage = None
                chain.window_iter = int(chain.window_iter) + 1
                chain.history.append_from_values(
                    iter=chain.iter,
                    state="UPDATE",
                    proposal=getattr(engine, "proposed_param", None),
                    logl=chain.last_logl,
                    accepted=accepted,
                    temperature=chain.temperature,
                )
            self._uuid_to_meta.pop(uuid, None)

    def _needs_stage_followup(self, chain_ids: Sequence[int] | None = None) -> bool:
        registry = self._ensure_registry()
        allowed = None if chain_ids is None else {int(c) for c in chain_ids}
        for chain in registry.all():
            if allowed is not None and int(chain.chain_id) not in allowed:
                continue
            if chain.open_stage is not None and int(chain.open_stage) > 0:
                if chain.engine.iterations < chain.engine.n_iterations:
                    return True
        return False

    def _all_finished(self) -> bool:
        registry = self._ensure_registry()
        return all(
            chain.engine.iterations >= chain.engine.n_iterations
            for chain in registry.all()
        )

    # ----------------------------------------------------------------- PT swap
    def _pair_chain_ids(self) -> list[tuple[int, int]]:
        ids = self._ensure_registry().ids()
        if len(ids) < 2:
            return []
        offset = int(self._exchange_offset) % len(ids)
        rotated = ids[offset:] + ids[:offset]
        pairs: list[tuple[int, int]] = []
        for ii in range(0, len(rotated) - 1, 2):
            pairs.append((rotated[ii], rotated[ii + 1]))
        self._exchange_offset = (self._exchange_offset + 1) % max(len(ids), 1)
        return pairs

    def _swap_acceptance(
        self, logl1: float, temp1: float, logl2: float, temp2: float
    ) -> bool:
        beta1 = 1.0 / float(temp1) if float(temp1) > 0 else 1.0
        beta2 = 1.0 / float(temp2) if float(temp2) > 0 else 1.0
        delta = (beta1 - beta2) * (float(logl2) - float(logl1))
        if delta >= 0.0:
            return True
        rng = self._exchange_rng or np.random.default_rng()
        return bool(rng.random() < np.exp(np.clip(delta, -700.0, 0.0)))

    def _attempt_swap(self, cid1: int, cid2: int) -> bool:
        registry = self._ensure_registry()
        chain1 = registry.get(cid1)
        chain2 = registry.get(cid2)
        if chain1.last_logl is None or chain2.last_logl is None:
            return False
        accepted = self._swap_acceptance(
            float(chain1.last_logl),
            chain1.temperature,
            float(chain2.last_logl),
            chain2.temperature,
        )
        if not accepted:
            return False
        # Swap state (params + logL); temperatures stay on the ladder slots.
        e1, e2 = chain1.engine, chain2.engine
        e1.param, e2.param = (
            np.asarray(e2.param, dtype=float).copy(),
            np.asarray(e1.param, dtype=float).copy(),
        )
        e1.last_loglikelihood, e2.last_loglikelihood = (
            e2.last_loglikelihood,
            e1.last_loglikelihood,
        )
        chain1.last_logl, chain2.last_logl = chain2.last_logl, chain1.last_logl
        return True

    def _should_exchange(self) -> bool:
        if not self._uses_pt():
            return False
        if int(self._exchange_interval) <= 0:
            return False
        registry = self._ensure_registry()
        for chain in registry.all():
            if chain.engine.iterations >= int(self._niters):
                continue
            if int(chain.window_iter) < int(self._exchange_interval):
                return False
            if chain.last_logl is None:
                return False
        return True

    def _do_exchange(self) -> dict[str, int]:
        attempted = 0
        accepted = 0
        for cid1, cid2 in self._pair_chain_ids():
            attempted += 1
            self._swap_attempts += 1
            if self._attempt_swap(cid1, cid2):
                accepted += 1
                self._swap_accepts += 1
        for chain in self._ensure_registry().all():
            chain.window_iter = 0
        self._logger.info(
            "%s exchange: attempted=%d accepted=%d",
            self.method,
            attempted,
            accepted,
        )
        return {"attempted": attempted, "accepted": accepted}

    def _run_half_or_all(
        self,
        *,
        timeout: float,
        chain_ids: Sequence[int] | None,
    ) -> int:
        """Propose → barrier → absorb for one subset of chains; return submitted count."""
        batch = self._propose_for_stages(stages="step", chain_ids=chain_ids)
        if not batch:
            return 0
        submitted = self._submit_sample_batch(batch)
        results = self.wait_for_generation(timeout=timeout)
        self.absorb_generation(results)
        while self._needs_stage_followup(chain_ids):
            stage_batch = self._propose_for_stages(
                stages="followup", chain_ids=chain_ids
            )
            if not stage_batch:
                break
            submitted += self._submit_sample_batch(stage_batch)
            results = self.wait_for_generation(timeout=timeout)
            self.absorb_generation(results)
        return submitted

    # ----------------------------------------------------------------- driver
    def run_adaptive(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Generation loop with a per-generation timeout."""
        if timeout is not None:
            generation_timeout = timeout
        timeout = generation_timeout
        self._require_redis(f"{self.method}.run_adaptive")
        self._ensure_seed_sequence()
        self._ensure_registry()
        if self._finished and self._all_finished() and not self._pending_uuids:
            return 0

        total_submitted = 0
        # Resume mid-barrier: drain pending first.
        if self._pending_uuids:
            results = self.wait_for_generation(timeout=timeout)
            self.absorb_generation(results)

        while not self._all_finished():
            if self._uses_half_ensemble():
                # Half 0 then half 1 — complementary partners stay fixed within a half.
                for half in (0, 1):
                    total_submitted += self._run_half_or_all(
                        timeout=timeout,
                        chain_ids=self._half_chain_ids(half),
                    )
            else:
                total_submitted += self._run_half_or_all(
                    timeout=timeout, chain_ids=None
                )

            if self._should_exchange():
                self._do_exchange()

            self._generation = max(
                (int(c.engine.iterations) for c in self._ensure_registry().all()),
                default=0,
            )
            self.checkpoint_at_barrier(reason=f"{self.method}_step_{self._generation}")

        self._finished = True
        self._summary = self._build_summary()
        # D13.6: additive DATABASE diagnostics (summary JSON + chain_history.csv).
        try:
            from jarvishep2.Sampling.diagnostics_export import export_mcmc_diagnostics

            task_dir = str(
                self.config.get("task_result_dir")
                or (self.config.get("Runtime") or {}).get("task_result_dir")
                or ""
            ).strip()
            if task_dir:
                written = export_mcmc_diagnostics(
                    task_dir,
                    summary=self._summary,
                    registry=self._registry,
                )
                if written and self._summary is not None:
                    self._summary["diagnostics"] = written
                for key, path in written.items():
                    self._logger.info("%s diagnostics %s → %s", self.method, key, path)
        except Exception as exc:
            self._logger.warning("failed to export MCMC diagnostics: %s", exc)
        self._logger.info(
            "%s finished: proposed=%d accepted=%d accept_rate=%.4f",
            self.method,
            self._total_proposed,
            self._total_accepted,
            (
                self._total_accepted / self._total_proposed
                if self._total_proposed
                else 0.0
            ),
        )
        return total_submitted

    def _build_summary(self) -> dict[str, Any]:
        registry = self._ensure_registry()
        chains = []
        logl_series: list[list[float]] = []
        ess_values: list[float] = []
        for chain in registry.all():
            accept = int(chain.accepted)
            reject = int(chain.rejected)
            total = accept + reject
            series = [
                float(ev.logl)
                for ev in chain.history.all()
                if ev.logl is not None and ev.accepted
            ]
            logl_series.append(series)
            ess = _effective_sample_size(series)
            if ess is not None:
                ess_values.append(float(ess))
            chains.append(
                {
                    "chain_id": chain.chain_id,
                    "iterations": int(chain.engine.iterations),
                    "accepted": accept,
                    "rejected": reject,
                    "accept_rate": (accept / total) if total else 0.0,
                    "last_logl": chain.last_logl,
                    "ess_logl": ess,
                }
            )
        rhat = _gelman_rubin_rhat(logl_series)
        payload: dict[str, Any] = {
            "method": self.method,
            "n_chains": self._nchains,
            "n_iters": self._niters,
            "total_proposed": self._total_proposed,
            "total_accepted": self._total_accepted,
            "accept_rate": (
                self._total_accepted / self._total_proposed
                if self._total_proposed
                else 0.0
            ),
            "rhat_logl": rhat,
            "ess_logl_mean": (
                float(sum(ess_values) / len(ess_values)) if ess_values else None
            ),
            "failed_samples": len(self._failed_uuids),
            "chains": chains,
        }
        if self._uses_half_ensemble():
            payload["half_ensemble"] = True
        if self._uses_pt():
            payload["swap_attempts"] = int(self._swap_attempts)
            payload["swap_accepts"] = int(self._swap_accepts)
            payload["temperature_ladder"] = list(self._temperature_ladder)
        return payload

    def summary(self) -> dict[str, Any]:
        if self._summary is None:
            self._summary = self._build_summary()
        return dict(self._summary)

    # ----------------------------------------------------------- checkpoint
    def repropose_unfinished(self) -> list[str]:
        if not self._pending_uuids:
            return []
        registry = self._ensure_registry()
        requeued: list[str] = []
        for uuid in sorted(self._pending_uuids):
            meta = self._uuid_to_meta.get(uuid)
            if meta is None:
                continue
            chain = registry.get(int(meta["chain_id"]))
            stage = int(meta["stage"])
            engine = chain.engine
            u = getattr(engine, "proposed_param", None)
            if u is None and hasattr(engine, "_stage_proposals"):
                u = engine._stage_proposals.get(stage)
            if u is None:
                continue
            sample = self._build_sample(np.asarray(u, dtype=np.float64))
            sample.uuid = uuid
            self._submit(sample)
            requeued.append(uuid)
        return requeued

    def export_runtime_state(self) -> dict[str, Any]:
        state = self._feedback_export_state()
        registry = self._ensure_registry()
        state.update(
            {
                "method": self.method,
                "nchains": self._nchains,
                "niters": self._niters,
                "finished": self._finished,
                "total_accepted": self._total_accepted,
                "total_proposed": self._total_proposed,
                "failed_uuids": list(self._failed_uuids),
                "uuid_to_meta": dict(self._uuid_to_meta),
                "summary": self._summary,
                "swap_attempts": int(self._swap_attempts),
                "swap_accepts": int(self._swap_accepts),
                "exchange_offset": int(self._exchange_offset),
                "temperature_ladder": list(self._temperature_ladder),
                "chains": [
                    {
                        "chain_id": c.chain_id,
                        "temperature": c.temperature,
                        "is_cold": c.is_cold,
                        "iter": c.iter,
                        "accepted": c.accepted,
                        "rejected": c.rejected,
                        "window_iter": int(c.window_iter),
                        "last_logl": c.last_logl,
                        "open_stage": c.open_stage,
                        "engine": c.engine.export_state(),
                    }
                    for c in registry.all()
                ],
            }
        )
        return state

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        self._feedback_import_state(state)
        self._finished = bool(state.get("finished", False))
        self._total_accepted = int(state.get("total_accepted", 0) or 0)
        self._total_proposed = int(state.get("total_proposed", 0) or 0)
        self._failed_uuids = list(state.get("failed_uuids") or [])
        self._uuid_to_meta = {
            str(k): dict(v) for k, v in dict(state.get("uuid_to_meta") or {}).items()
        }
        self._summary = state.get("summary")
        self._swap_attempts = int(state.get("swap_attempts", 0) or 0)
        self._swap_accepts = int(state.get("swap_accepts", 0) or 0)
        self._exchange_offset = int(state.get("exchange_offset", 0) or 0)
        ladder = state.get("temperature_ladder")
        if ladder is not None:
            self._temperature_ladder = [float(x) for x in ladder]
        raw_chains = list(state.get("chains") or [])
        if not raw_chains:
            self._registry = None
            return
        scales = self._normalize_scales() if self._nchains else []
        n_spawn = max(1, len(raw_chains) + (1 if self._uses_pt() else 0))
        children = self._ensure_seed_sequence().spawn(n_spawn)
        chains: list[ChainRuntime] = []
        for item in raw_chains:
            if not isinstance(item, Mapping):
                continue
            cid = int(item.get("chain_id", 0))
            rng = np.random.default_rng(children[min(cid, len(children) - 1)])
            scale = scales[cid] if cid < len(scales) else 0.1
            engine = self._make_engine(cid, scale, rng)
            eng_state = dict(item.get("engine") or {})
            if eng_state:
                engine.import_state(eng_state)
            chains.append(
                ChainRuntime(
                    chain_id=cid,
                    engine=engine,
                    temperature=float(item.get("temperature", 1.0)),
                    is_cold=bool(item.get("is_cold", cid == 0)),
                    iter=int(item.get("iter", 0) or 0),
                    accepted=int(item.get("accepted", 0) or 0),
                    rejected=int(item.get("rejected", 0) or 0),
                    window_iter=int(item.get("window_iter", 0) or 0),
                    last_logl=item.get("last_logl"),
                    open_stage=item.get("open_stage"),
                )
            )
        self._registry = ChainRegistry(chains)
        if self._uses_pt() and len(children) > len(raw_chains):
            self._exchange_rng = np.random.default_rng(children[len(raw_chains)])


def _gelman_rubin_rhat(series: list[list[float]]) -> float | None:
    """Simple Gelman–Rubin R-hat on scalar logL traces (diagnostic only)."""
    cleaned = [np.asarray(s, dtype=float) for s in series if len(s) >= 4]
    if len(cleaned) < 2:
        return None
    n = min(len(s) for s in cleaned)
    if n < 4:
        return None
    # Use second half (burn-in discard heuristic).
    start = n // 2
    chains = np.stack([s[start:n] for s in cleaned], axis=0)
    m, n_eff = chains.shape
    chain_means = chains.mean(axis=1)
    chain_vars = chains.var(axis=1, ddof=1)
    b = n_eff * float(np.var(chain_means, ddof=1))
    w = float(np.mean(chain_vars))
    if w <= 0:
        return None
    var_hat = ((n_eff - 1) / n_eff) * w + b / n_eff
    return float(np.sqrt(var_hat / w))


def _effective_sample_size(series: Sequence[float] | list[float]) -> float | None:
    """Autocorrelation-time ESS on a scalar logL chain (D13.7d).

    Uses the positive initial sequence of lag autocorrelations:
    ``ESS ≈ N / (1 + 2 Σ ρ_t)`` until the first non-positive ρ. Diagnostic
    only — not a substitute for multi-parameter ESS tools (ArviZ/emcee).
    """
    x = np.asarray(list(series), dtype=float)
    n = int(x.size)
    if n < 4:
        return None
    # Burn-in heuristic aligned with R-hat (second half).
    x = x[n // 2 :]
    n = int(x.size)
    if n < 4:
        return None
    x = x - float(x.mean())
    var0 = float(np.dot(x, x) / n)
    if var0 <= 1e-16:
        return float(n)
    rho_sum = 0.0
    max_lag = n // 2
    for lag in range(1, max_lag + 1):
        cov = float(np.dot(x[:-lag], x[lag:]) / n)
        rho = cov / var0
        if rho <= 0.0:
            break
        rho_sum += rho
    tau = 1.0 + 2.0 * rho_sum
    tau = max(1.0, min(tau, float(n)))
    return float(n / tau)


def create_mcmc() -> MCMCSampler:
    return MCMCSampler("MCMC")


def create_ammcmc() -> MCMCSampler:
    return MCMCSampler("AMMCMC")


def create_am() -> MCMCSampler:
    return MCMCSampler("AM")


def create_dram() -> MCMCSampler:
    return MCMCSampler("DRAM")


def create_ensemble() -> MCMCSampler:
    return MCMCSampler("EnsembleMCMC")


def create_ensemble_alias() -> MCMCSampler:
    return MCMCSampler("Ensemble")


def create_demcmc() -> MCMCSampler:
    return MCMCSampler("DEMCMC")


def create_ptmcmc() -> MCMCSampler:
    return MCMCSampler("PTMCMC")


def create_pt() -> MCMCSampler:
    return MCMCSampler("PT")


def create_pt_ensemble() -> MCMCSampler:
    return MCMCSampler("PTEnsemble")


__all__ = [
    "MCMCSampler",
    "create_am",
    "create_ammcmc",
    "create_demcmc",
    "create_dram",
    "create_ensemble",
    "create_ensemble_alias",
    "create_mcmc",
    "create_pt",
    "create_pt_ensemble",
    "create_ptmcmc",
    "mcmc_sample_uuid",
]
