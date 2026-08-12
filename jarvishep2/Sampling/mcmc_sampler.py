#!/usr/bin/env python3
"""Feedback-driven MCMC family samplers (D13.2 + D13.3).

Ports V1 chain science (``Sampling/Source/MCMC/*`` engines) onto
:class:`FeedbackSampler`.

**Base runtime design (independent chains)** — ToyMCMC / MCMC / AM / DRAM and
any method that leaves ``_uses_half_ensemble`` / ``_uses_pt`` false:

* single shared Redis **task** queue (workers steal freely);
* per-chain **feedback** shards (``hep:feedback:chain:{id}``);
* **async event loop**: at most one in-flight evaluation per chain; as soon as
  any chain returns, absorb it and immediately re-submit that chain;
* no generation barrier across chains (Markov order stays *inside* each chain).

Coupled methods keep the classical barrier path:

* ensemble / DE → half-ensemble generation barriers;
* PT → full barriers at temperature-exchange points.

Concrete methods only override thin hooks (``_make_engine``,
``_configure_method``, progress observers). Never imports ``jarvishep``.
"""

from __future__ import annotations

import hashlib
import pickle
import time
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from jarvishep2.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep2.Sampling.Source.MCMC.chain_history import ChainHistory
from jarvishep2.Sampling.Source.MCMC.config_contract import (
    normalize_proposal_scales,
    parse_common_chain_counts,
    parse_proposal_scale_value,
)
from jarvishep2.Sampling.Source.MCMC.mcmc_chain import MCMCChain
from jarvishep2.Sampling.feedback_sampler import FeedbackSampler
from jarvishep2.Sampling.sampling_utils import evaluate_selection, physical_from_u
from jarvishep2.Sampling.variables import load_variables
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.redis_queue import FEEDBACK_QUEUE, chain_feedback_queue
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample

_FAILED_LOGL = -1.0e300
_MCMC_NATIVE_STATE_FORMAT = "jarvis-hep.mcmc-native"
_MCMC_NATIVE_STATE_VERSION = 1


def _pickle_mcmc_engine(engine: Any) -> tuple[bytes, bool]:
    """Pickle one engine after detaching callbacks into the live sampler.

    Ensemble/DE engines hold a bound ``population_getter``. Pickling that
    callback would recursively capture the sampler, Redis, and logger. The
    callback is restored immediately after serialization and rebound to the
    new sampler after deserialization.
    """
    attrs = getattr(engine, "__dict__", None)
    missing = object()
    getter = missing
    had_population_getter = False
    if isinstance(attrs, dict) and "_population_getter" in attrs:
        getter = attrs.pop("_population_getter")
        had_population_getter = True
    try:
        return (
            pickle.dumps(engine, protocol=pickle.HIGHEST_PROTOCOL),
            had_population_getter,
        )
    finally:
        if isinstance(attrs, dict) and getter is not missing:
            attrs["_population_getter"] = getter


def _unpickle_mcmc_engine(
    blob: bytes,
    sampler: "MCMCBaseSampler",
    *,
    requires_population_getter: bool = False,
) -> Any:
    """Restore one engine and reconnect any sampler-owned population callback."""
    engine = pickle.loads(blob)
    if requires_population_getter or type(engine).__name__ in {
        "EnsembleChain",
        "DEMCMCChain",
    }:
        engine._population_getter = sampler._population_for
    return engine

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


class MCMCBaseSampler(FeedbackSampler):
    """Common feedback/checkpoint runtime for one concrete MCMC sampler.

    Concrete methods live in separate modules and override the small method
    hooks below (configuration, engine construction, ensemble/temperature
    capabilities). The base class owns transport, the independent-chain async
    pipeline (base design), barrier mode for coupled methods, checkpointing,
    diagnostics, and shared MH bookkeeping.

    **Contract for new independent MCMC methods**

    1. Subclass ``MCMCBaseSampler`` (or a thin science base like AM).
    2. Set ``method`` and implement ``_make_engine`` (+ optional Bounds hooks).
    3. Leave ``_uses_half_ensemble`` / ``_uses_pt`` False — the async per-chain
       pipeline is the default (and only) independent-chain design.
    4. Optional observers: ``_on_run_started``, ``_on_generation_completed``
       (fires after each absorbed transition in async mode).
    """

    method = "MCMC"
    _checkpoint_saved_attributes = frozenset(
        {
            "method",
            "_pending_uuids",
            "_seed",
            "_generation",
            "_on_failure",
            "_uuid_to_meta",
            "_finished",
            "_total_accepted",
            "_total_proposed",
            "_failed_uuids",
            "_summary",
        }
    )
    _checkpoint_excluded_attributes = frozenset(
        {
            "_logger",
            "vars",
            "_dim",
            "_nchains",
            "_niters",
            "_proposal_scales",
            "_selectionexp",
            "_registry",
            "_seed_seq",
            "_batch_size",
            "_sampling_checkpoint_interval_sec",
            "_last_sampling_checkpoint_at",
            "_last_sampling_checkpoint_generation",
        }
    )

    def __init__(self) -> None:
        super().__init__()
        self.method = str(type(self).method)
        self._logger = get_jarvis_logger(f"sampler.{self.method.lower()}")
        self.vars: list = []
        self._mapper_pipeline = None
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
        self._summary: dict[str, Any] | None = None
        # MCMC checkpointing is time-cadenced at a feedback barrier. The
        # durable archive barrier is reserved for explicit/final checkpoints.
        self._sampling_checkpoint_interval_sec = 30.0
        self._last_sampling_checkpoint_at = 0.0
        self._last_sampling_checkpoint_generation = -1

    # ------------------------------------------------------------------ config
    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        from jarvishep2.mapper import MapperPipeline

        self._mapper_pipeline = MapperPipeline.from_config(self.config)
        self._dim = len(self.vars)
        if self._dim < 1:
            raise ValueError(f"{self.method} requires at least one Sampling.Variable")
        bounds = sampling.get("Bounds") or {}
        if not isinstance(bounds, Mapping):
            raise ValueError("Sampling.Bounds must be a mapping for MCMC methods")
        self._nchains, self._niters = parse_common_chain_counts(bounds)
        self._proposal_scales = parse_proposal_scale_value(bounds, default=0.1)
        self._selectionexp = sampling.get("selection")
        # seed lives under Bounds only.
        seed = int(bounds.get("seed", 0) or 0)
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        self._sampling_checkpoint_interval_sec = max(
            1.0,
            float(runtime.get("checkpoint_heartbeat_sec", 30.0) or 30.0),
        )
        self._last_sampling_checkpoint_at = 0.0
        self._last_sampling_checkpoint_generation = -1
        self._init_seed_sequence(seed)
        on_fail = bounds.get("on_failure", sampling.get("on_failure", "reject"))
        self._on_failure = str(on_fail or "reject").lower()
        if self._nchains < workers:
            self._logger.warning(
                "%s: chains (%d) < workers (%d) — workers may idle; "
                "prefer chains >= workers for the async pipeline",
                self.method,
                self._nchains,
                workers,
            )

        self._configure_method(bounds)

    def _normalize_scales(self) -> list[float]:
        return normalize_proposal_scales(
            self._proposal_scales,
            nchains=self._nchains,
            sampler_method=self.method,
        )

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        """Apply method-specific Bounds configuration in a concrete subclass."""
        return None

    def _on_run_started(self) -> None:
        """Hook for method-specific runtime observers such as progress logs."""
        return None

    def _on_generation_completed(self) -> None:
        """Hook called after feedback is absorbed at a generation barrier."""
        return None

    def _on_runtime_state_imported(self) -> None:
        """Reset process-local observers after checkpoint state is restored."""
        return None

    def _export_method_state(self) -> dict[str, Any]:
        """Return checkpoint fields owned by the concrete sampler."""
        return {}

    def _import_method_state(self, state: Mapping[str, Any]) -> None:
        """Restore checkpoint fields owned by the concrete sampler."""
        return None

    def checkpoint_runtime_state(self, *, safe: bool) -> dict[str, Any]:
        """Export current MCMC state at a feedback barrier.

        ``safe=False`` means the Archiver has not acknowledged every submitted
        UUID yet. Independent-chain async mode may checkpoint with in-flight
        samples: ``pending_uuids`` + ``uuid_to_meta`` + engine state are all in
        the snapshot and resume re-drains outstanding feedback. Coupled barrier
        methods still require an empty pending set for a fresh export.
        """
        if safe:
            return super().checkpoint_runtime_state(safe=True)
        if self.at_safe_barrier() or self._uses_async_chain_pipeline():
            return self.export_runtime_state()
        return super().checkpoint_runtime_state(safe=False)

    def _export_native_mcmc_state(self) -> bytes | None:
        """Serialize chain engines/history without Redis or sampler callbacks."""
        registry = self._ensure_registry()
        chains: list[dict[str, Any]] = []
        try:
            for chain in registry.all():
                engine = chain.engine
                engine_blob, requires_population_getter = _pickle_mcmc_engine(engine)
                chains.append(
                    {
                        "chain_id": int(chain.chain_id),
                        "temperature": float(chain.temperature),
                        "is_cold": bool(chain.is_cold),
                        "iter": int(chain.iter),
                        "accepted": int(chain.accepted),
                        "rejected": int(chain.rejected),
                        "window_iter": int(chain.window_iter),
                        "last_logl": chain.last_logl,
                        "open_stage": chain.open_stage,
                        "meta": dict(chain.meta),
                        "history": chain.history,
                        "engine_type": (
                            f"{type(engine).__module__}.{type(engine).__qualname__}"
                        ),
                        "engine_blob": engine_blob,
                        "requires_population_getter": requires_population_getter,
                    }
                )
            exchange_rng = getattr(self, "_exchange_rng", None)
            payload = {
                "format": _MCMC_NATIVE_STATE_FORMAT,
                "version": _MCMC_NATIVE_STATE_VERSION,
                "method": self.method,
                "generation": int(self._generation),
                "chains": chains,
                "method_state": self._export_method_state(),
                "exchange_rng_state": (
                    None
                    if exchange_rng is None
                    else exchange_rng.bit_generator.state
                ),
            }
            return pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as exc:
            self._logger.warning(
                "%s: native MCMC checkpoint pickle failed; using explicit state: %s",
                self.method,
                exc,
            )
            return None

    def _import_native_mcmc_state(self, blob: bytes) -> bool:
        """Restore native chain state and rebind sampler-owned callbacks."""
        try:
            payload = pickle.loads(blob)
            if not isinstance(payload, Mapping):
                raise ValueError("native MCMC checkpoint is not a mapping")
            if payload.get("format") != _MCMC_NATIVE_STATE_FORMAT:
                raise ValueError("unknown native MCMC checkpoint format")
            if int(payload.get("version", -1)) != _MCMC_NATIVE_STATE_VERSION:
                raise ValueError("unsupported native MCMC checkpoint version")
            if str(payload.get("method", self.method)) != self.method:
                raise ValueError("native MCMC checkpoint method mismatch")

            raw_chains = list(payload.get("chains") or [])
            chains: list[ChainRuntime] = []
            for item in raw_chains:
                if not isinstance(item, Mapping):
                    raise ValueError("native MCMC chain entry is not a mapping")
                engine_blob = item.get("engine_blob")
                if not isinstance(engine_blob, (bytes, bytearray)):
                    raise ValueError("native MCMC chain is missing engine blob")
                engine = _unpickle_mcmc_engine(
                    bytes(engine_blob),
                    self,
                    requires_population_getter=bool(
                        item.get("requires_population_getter", False)
                    ),
                )
                engine_type = (
                    f"{type(engine).__module__}.{type(engine).__qualname__}"
                )
                if str(item.get("engine_type", engine_type)) != engine_type:
                    raise ValueError("native MCMC engine type mismatch")
                history = item.get("history")
                if not isinstance(history, ChainHistory):
                    raise ValueError("native MCMC chain history is invalid")
                chains.append(
                    ChainRuntime(
                        chain_id=int(item.get("chain_id", 0)),
                        engine=engine,
                        temperature=float(item.get("temperature", 1.0)),
                        is_cold=bool(item.get("is_cold", False)),
                        iter=int(item.get("iter", 0) or 0),
                        accepted=int(item.get("accepted", 0) or 0),
                        rejected=int(item.get("rejected", 0) or 0),
                        window_iter=int(item.get("window_iter", 0) or 0),
                        last_logl=item.get("last_logl"),
                        history=history,
                        meta=dict(item.get("meta") or {}),
                        open_stage=item.get("open_stage"),
                    )
                )
            if not chains:
                raise ValueError("native MCMC checkpoint has no chains")
            self._registry = ChainRegistry(chains)

            method_state = payload.get("method_state")
            if isinstance(method_state, Mapping):
                self._import_method_state(method_state)
            exchange_rng_state = payload.get("exchange_rng_state")
            if self._uses_pt() and exchange_rng_state is not None:
                self._exchange_rng = np.random.default_rng()
                self._exchange_rng.bit_generator.state = exchange_rng_state
            return True
        except Exception as exc:
            self._logger.warning(
                "%s: native MCMC checkpoint restore failed; using explicit state: %s",
                self.method,
                exc,
            )
            return False

    def _checkpoint_at_sampling_barrier(self, *, reason: str) -> bool:
        """Save MCMC state without waiting for Archiver acknowledgements."""
        # Async independent chains almost always have in-flight work; still
        # allow time-cadenced snapshots (pending uuids are part of the state).
        if not self.at_safe_barrier() and not self._uses_async_chain_pipeline():
            return False
        now = time.monotonic()
        if (
            self._last_sampling_checkpoint_generation >= 0
            and now - self._last_sampling_checkpoint_at
            < self._sampling_checkpoint_interval_sec
        ):
            return False
        written = self.persist_runtime_checkpoint(force=True, reason=reason)
        if written:
            self._last_sampling_checkpoint_at = now
            self._last_sampling_checkpoint_generation = int(self._generation)
        return bool(written)

    def _uses_half_ensemble(self) -> bool:
        return False

    def _uses_pt(self) -> bool:
        return False

    def _chains_are_independent(self) -> bool:
        """True when chains share no cross-chain science barrier."""
        return not self._uses_half_ensemble() and not self._uses_pt()

    def _uses_async_chain_pipeline(self) -> bool:
        """Per-chain 1-inflight event loop — base design for independent chains."""
        return self._chains_are_independent()

    def _uses_sharded_feedback(self) -> bool:
        """Per-chain feedback shards — base design with the async pipeline."""
        return self._chains_are_independent()

    def _feedback_queue_for_chain(self, chain_id: int) -> str:
        return chain_feedback_queue(int(chain_id))

    def _pending_feedback_queues(self) -> list[str]:
        """Feedback lists that may currently hold a pending result."""
        if not self._uses_sharded_feedback():
            return [FEEDBACK_QUEUE]
        queues: list[str] = []
        seen: set[str] = set()
        for uuid in self._pending_uuids:
            meta = self._uuid_to_meta.get(str(uuid))
            if meta is None:
                continue
            key = self._feedback_queue_for_chain(int(meta["chain_id"]))
            if key not in seen:
                seen.add(key)
                queues.append(key)
        return queues or [FEEDBACK_QUEUE]

    def _inflight_chain_ids(self) -> set[int]:
        out: set[int] = set()
        for uuid in self._pending_uuids:
            meta = self._uuid_to_meta.get(str(uuid))
            if meta is not None:
                out.add(int(meta["chain_id"]))
        return out

    def _ready_chain_ids(self) -> list[int]:
        """Chains that can accept a new proposal (no in-flight evaluation)."""
        inflight = self._inflight_chain_ids()
        ready: list[int] = []
        for chain in self._ensure_registry().all():
            cid = int(chain.chain_id)
            if cid in inflight:
                continue
            if int(chain.engine.iterations) >= int(chain.engine.n_iterations):
                continue
            ready.append(cid)
        return ready

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
                physical = physical_from_u(u, self.vars, self._mapper_pipeline)
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
        # MCMC wire contract (control → Redis task → Worker):
        #   uuid + u_coords + chain_id [+ step/stage] (+ optional feedback_queue)
        # Observables do *not* cross Redis; the Worker stamps identity into
        # sample info / DATABASE after bind_params (Sample.apply_mcmc_identity).
        sample.chain_id = int(chain.chain_id)
        sample.step = int(step)
        sample.stage = int(stage)
        if self._uses_sharded_feedback():
            sample.feedback_queue = self._feedback_queue_for_chain(chain.chain_id)
        else:
            sample.feedback_queue = None
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
        results = self.wait_for_generation(
            timeout=timeout, queues=self._pending_feedback_queues()
        )
        self.absorb_generation(results)
        while self._needs_stage_followup(chain_ids):
            stage_batch = self._propose_for_stages(
                stages="followup", chain_ids=chain_ids
            )
            if not stage_batch:
                break
            submitted += self._submit_sample_batch(stage_batch)
            results = self.wait_for_generation(
                timeout=timeout, queues=self._pending_feedback_queues()
            )
            self.absorb_generation(results)
        return submitted

    def _top_up_ready_chains(self) -> int:
        """Submit at most one evaluation for every chain without in-flight work.

        Order: delayed-rejection follow-ups first, then new Metropolis steps.
        Ready sets exclude chains that already have a pending Redis evaluation,
        so each chain stays strictly sequential while chains run in parallel.
        """
        ready = self._ready_chain_ids()
        if not ready:
            return 0
        submitted = 0
        followup = self._propose_for_stages(stages="followup", chain_ids=ready)
        if followup:
            submitted += self._submit_sample_batch(followup)
        # Recompute ready: follow-up submit marks those chains in-flight.
        ready = self._ready_chain_ids()
        if not ready:
            return submitted
        steps = self._propose_for_stages(stages="step", chain_ids=ready)
        if steps:
            submitted += self._submit_sample_batch(steps)
        return submitted

    def _sync_generation_from_chains(self) -> None:
        self._generation = max(
            (int(c.engine.iterations) for c in self._ensure_registry().all()),
            default=0,
        )

    def _run_async_independent_chains(self, *, timeout: float) -> int:
        """Base high-throughput driver for independent multi-chain MCMC.

        Pipeline shape::

            while work remains:
                top-up every idle chain  →  shared hep:task_queue
                wait for ANY feedback   →  per-chain shard
                absorb that one chain
                immediately allow that chain to re-enter top-up

        In-flight depth ≈ number of unfinished chains (≤ num_chains). Workers
        stay busy whenever evaluations dominate control-plane cost.
        """
        total_submitted = 0
        self._logger.info(
            "%s independent-chain pipeline: chains=%d iters=%d sharded_feedback=True",
            self.method,
            self._nchains,
            self._niters,
        )

        while True:
            total_submitted += self._top_up_ready_chains()

            if not self._pending_uuids:
                if self._all_finished():
                    break
                # No pending and not finished: invariant broken (no proposals
                # can be made). Never mark the run successful.
                if not self._ready_chain_ids():
                    iters = [
                        (
                            int(c.chain_id),
                            int(c.engine.iterations),
                            int(c.engine.n_iterations),
                            c.open_stage,
                        )
                        for c in self._ensure_registry().all()
                    ]
                    raise RuntimeError(
                        f"{self.method}: async pipeline stalled with no pending "
                        f"feedback and no ready chains while work remains "
                        f"(chain_id, iterations, n_iterations, open_stage)={iters}"
                    )
                continue

            record = self.wait_for_any_feedback(
                timeout=timeout,
                queues=self._pending_feedback_queues(),
            )
            self.absorb_generation([record])
            self._sync_generation_from_chains()
            self._on_generation_completed()
            self._checkpoint_at_sampling_barrier(
                reason=f"{self.method}_async_step_{self._generation}"
            )

        return total_submitted

    def _run_barrier_generations(self, *, timeout: float) -> int:
        """Classical generation barrier path (ensemble half-steps / PT)."""
        total_submitted = 0
        # Resume mid-barrier: drain pending first.
        if self._pending_uuids:
            results = self.wait_for_generation(
                timeout=timeout, queues=self._pending_feedback_queues()
            )
            self.absorb_generation(results)
            self._on_generation_completed()

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

            self._sync_generation_from_chains()
            self._on_generation_completed()
            self._checkpoint_at_sampling_barrier(
                reason=f"{self.method}_sampling_step_{self._generation}"
            )
        return total_submitted

    def _finalize_run(self) -> None:
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
        # The normal loop uses a feedback-only checkpoint. The final snapshot
        # is the one place where MCMC asks for the durable archive barrier.
        if self._save_checkpoint_callback is not None:
            self.checkpoint_at_barrier(reason=f"{self.method}_finished")
        self._logger.info(
            "%s finished: proposed=%d accepted=%d accept_rate=%.4f async=%s",
            self.method,
            self._total_proposed,
            self._total_accepted,
            (
                self._total_accepted / self._total_proposed
                if self._total_proposed
                else 0.0
            ),
            self._uses_async_chain_pipeline(),
        )

    # ----------------------------------------------------------------- driver
    def run_adaptive(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Run the MCMC feedback loop (async template or generation barriers)."""
        if timeout is not None:
            generation_timeout = timeout
        timeout = generation_timeout
        self._require_redis(f"{self.method}.run_adaptive")
        self._ensure_seed_sequence()
        self._ensure_registry()
        self._on_run_started()
        if self._finished and self._all_finished() and not self._pending_uuids:
            return 0

        if self._uses_async_chain_pipeline():
            total_submitted = self._run_async_independent_chains(timeout=timeout)
        else:
            total_submitted = self._run_barrier_generations(timeout=timeout)

        self._finalize_run()
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
        payload["async_chains"] = bool(self._uses_async_chain_pipeline())
        payload["sharded_feedback"] = bool(self._uses_sharded_feedback())
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
            sample.chain_id = int(meta["chain_id"])
            sample.step = int(meta.get("step", 0))
            sample.stage = int(stage)
            if self._uses_sharded_feedback():
                sample.feedback_queue = self._feedback_queue_for_chain(
                    int(meta["chain_id"])
                )
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
                "proposal_scales": self._normalize_scales(),
                "finished": self._finished,
                "total_accepted": self._total_accepted,
                "total_proposed": self._total_proposed,
                "failed_uuids": list(self._failed_uuids),
                "uuid_to_meta": dict(self._uuid_to_meta),
                "summary": self._summary,
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
        state.update(self._export_method_state())
        native_blob = self._export_native_mcmc_state()
        if native_blob is not None:
            state.update(
                {
                    "native_mcmc_state_blob": native_blob,
                    "native_mcmc_state_format": _MCMC_NATIVE_STATE_FORMAT,
                    "native_mcmc_state_version": _MCMC_NATIVE_STATE_VERSION,
                }
            )
        return state

    def _validate_checkpoint_against_config(self, state: Mapping[str, Any]) -> None:
        """Refuse resume when current Bounds no longer match the checkpoint.

        Prevents YAML edits (num_chains / num_iters / proposal_scale / seed)
        from silently continuing under an old runtime shape while summary still
        reports the new card.
        """
        method = str(state.get("method") or "").strip()
        if method and method != str(self.method):
            raise ValueError(
                f"{self.method}: checkpoint method {method!r} does not match "
                f"current Method {self.method!r}"
            )

        ck_nchains = state.get("nchains")
        if ck_nchains is not None and int(ck_nchains) != int(self._nchains):
            raise ValueError(
                f"{self.method}: checkpoint num_chains={int(ck_nchains)} does not "
                f"match current Bounds.num_chains={int(self._nchains)}"
            )

        ck_niters = state.get("niters")
        if ck_niters is not None and int(ck_niters) != int(self._niters):
            raise ValueError(
                f"{self.method}: checkpoint num_iters={int(ck_niters)} does not "
                f"match current Bounds.num_iters={int(self._niters)}"
            )

        ck_seed = state.get("seed")
        if ck_seed is not None and int(ck_seed) != int(self._seed):
            raise ValueError(
                f"{self.method}: checkpoint seed={int(ck_seed)} does not match "
                f"current Bounds.seed={int(self._seed)}"
            )

        raw_chains = list(state.get("chains") or [])
        if raw_chains and len(raw_chains) != int(self._nchains):
            raise ValueError(
                f"{self.method}: checkpoint has {len(raw_chains)} chain(s) but "
                f"current Bounds.num_chains={int(self._nchains)}"
            )

        ck_scales = state.get("proposal_scales")
        if ck_scales is not None:
            current = [float(x) for x in self._normalize_scales()]
            try:
                saved = [float(x) for x in list(ck_scales)]
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{self.method}: checkpoint proposal_scales is not numeric"
                ) from exc
            if len(saved) != len(current) or any(
                abs(a - b) > 1e-12 for a, b in zip(saved, current)
            ):
                raise ValueError(
                    f"{self.method}: checkpoint proposal_scale {saved!r} does not "
                    f"match current Bounds.proposal_scale {current!r}"
                )

    def _assert_registry_matches_config(self) -> None:
        """Validate restored engines against the live Bounds shape."""
        if self._registry is None:
            return
        chains = self._registry.all()
        if len(chains) != int(self._nchains):
            raise ValueError(
                f"{self.method}: restored {len(chains)} chain(s) but "
                f"Bounds.num_chains={int(self._nchains)}"
            )
        for chain in chains:
            eng_niters = getattr(chain.engine, "n_iterations", None)
            if eng_niters is not None and int(eng_niters) != int(self._niters):
                raise ValueError(
                    f"{self.method}: restored chain {chain.chain_id} has "
                    f"n_iterations={int(eng_niters)} but "
                    f"Bounds.num_iters={int(self._niters)}"
                )

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        self._validate_checkpoint_against_config(state)
        try:
            self._feedback_import_state(state)
            self._on_runtime_state_imported()
            self._finished = bool(state.get("finished", False))
            self._total_accepted = int(state.get("total_accepted", 0) or 0)
            self._total_proposed = int(state.get("total_proposed", 0) or 0)
            self._failed_uuids = list(state.get("failed_uuids") or [])
            self._uuid_to_meta = {
                str(k): dict(v)
                for k, v in dict(state.get("uuid_to_meta") or {}).items()
            }
            self._summary = state.get("summary")
            self._import_method_state(state)
            native_blob = state.get("native_mcmc_state_blob")
            if isinstance(native_blob, (bytes, bytearray)) and self._import_native_mcmc_state(
                bytes(native_blob)
            ):
                self._assert_registry_matches_config()
                return
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
            self._assert_registry_matches_config()
        except Exception:
            # Leave a clean sampler rather than a half-restored registry that
            # still reports the new YAML Bounds in summary().
            self._registry = None
            self._finished = False
            self._pending_uuids = set()
            self._uuid_to_meta = {}
            self._summary = None
            raise


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


# Compatibility name for callers that imported the old monolithic class.
# New code should import ``MCMCBaseSampler`` or the concrete class from its
# method module (for example ``Sampling.toymcmc.ToyMCMCSampler``).
MCMCSampler = MCMCBaseSampler


def create_mcmc() -> MCMCBaseSampler:
    from jarvishep2.Sampling.mcmc import create_mcmc as factory

    return factory()


def create_toymcmc() -> MCMCBaseSampler:
    """Create the lightweight baseline MCMC profile used by toy/test cards."""
    from jarvishep2.Sampling.toymcmc import create_toymcmc as factory

    return factory()


def create_ammcmc() -> MCMCBaseSampler:
    from jarvishep2.Sampling.ammcmc import create_ammcmc as factory

    return factory()


def create_am() -> MCMCBaseSampler:
    from jarvishep2.Sampling.am import create_am as factory

    return factory()


def create_dram() -> MCMCBaseSampler:
    from jarvishep2.Sampling.dram import create_dram as factory

    return factory()


def create_ensemble() -> MCMCBaseSampler:
    from jarvishep2.Sampling.ensemble_mcmc import create_ensemble as factory

    return factory()


def create_ensemble_alias() -> MCMCBaseSampler:
    from jarvishep2.Sampling.ensemble_mcmc import create_ensemble_alias as factory

    return factory()


def create_demcmc() -> MCMCBaseSampler:
    from jarvishep2.Sampling.demcmc import create_demcmc as factory

    return factory()


def create_ptmcmc() -> MCMCBaseSampler:
    from jarvishep2.Sampling.ptmcmc import create_ptmcmc as factory

    return factory()


def create_pt() -> MCMCBaseSampler:
    from jarvishep2.Sampling.ptmcmc import create_pt as factory

    return factory()


def create_pt_ensemble() -> MCMCBaseSampler:
    from jarvishep2.Sampling.ptensemble import create_pt_ensemble as factory

    return factory()


__all__ = [
    "MCMCBaseSampler",
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
    "create_toymcmc",
    "mcmc_sample_uuid",
]
