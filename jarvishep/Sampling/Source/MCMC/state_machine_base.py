#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
from abc import ABCMeta, abstractmethod
from collections import deque
from enum import Enum
from typing import Any, Dict, Sequence

from jarvishep.log_kv import format_two_column_log
from jarvishep.Sampling.sampler import SamplingVirtial
from jarvishep.sample import Sample

from .chain_runtime import ChainRegistry, ChainRuntime
from .controller import (
    MCMCControlGuard,
    MCMCControlPatch,
    MCMCControllerProtocol,
    NoopMCMCController,
)
from .engine_contract import validate_engine_contract
from .metrics_bus import MCMCMetricsBus, MCMCMetricsFrame


class MCMCState(str, Enum):
    INIT = "INIT"
    PROPOSE = "PROPOSE"
    SUBMIT = "SUBMIT"
    WAIT = "WAIT"
    UPDATE = "UPDATE"
    EXCHANGE = "EXCHANGE"
    TERMINATE = "TERMINATE"
    FAILED = "FAILED"


class MCMCStateMachineBase(SamplingVirtial):
    __metaclass__ = ABCMeta

    def __init__(self) -> None:
        super().__init__()
        self._state: MCMCState = MCMCState.INIT
        self._chain_registry: ChainRegistry | None = None
        self._pending_futures: set = set()
        self._future_to_sample: Dict[Any, Sample] = {}
        self._future_to_chain_id: Dict[Any, int] = {}
        self._inflight_chain_ids: set[int] = set()
        self._ready_queue = deque()
        self._ready_set: set[int] = set()
        self._state_machine_ready = False
        self._controller: MCMCControllerProtocol = NoopMCMCController()
        self._metrics_bus = MCMCMetricsBus()
        self._metrics_step = 0

    def set_controller(self, controller: MCMCControllerProtocol | None) -> None:
        self._controller = controller if controller is not None else NoopMCMCController()

    def _transition(self, new_state: MCMCState, reason: str = "") -> None:
        old = self._state
        self._state = new_state
        if getattr(self, "logger", None) is None:
            return
        if reason:
            self.logger.info(f"MCMC State -> {old.value} -> {new_state.value} ({reason})")
        else:
            self.logger.info(f"MCMC State -> {old.value} -> {new_state.value}")

    @property
    def state(self) -> MCMCState:
        return self._state

    def map_point_into_distribution(self, row) -> Dict[str, float]:
        result = {}
        for ii in range(len(row)):
            result[self.vars[ii].name] = self.vars[ii].map_standard_random_to_distribution(row[ii])
        return result

    def history_all(self, chain_id: int):
        chain = self._must_registry().get(chain_id)
        return chain.history.all()

    def history_tail(self, chain_id: int, n: int):
        chain = self._must_registry().get(chain_id)
        return chain.history.tail(n)

    def is_cold_chain(self, chain_id: int) -> bool:
        return self._must_registry().is_cold(chain_id)

    def chain_snapshot(self, chain_id: int) -> Dict[str, Any]:
        return self._must_registry().get(chain_id).snapshot()

    def metrics_all(self):
        return self._metrics_bus.all()

    def metrics_tail(self, n: int):
        return self._metrics_bus.tail(n)

    def metrics_latest(self):
        return self._metrics_bus.latest()

    def _must_registry(self) -> ChainRegistry:
        if self._chain_registry is None:
            raise RuntimeError("MCMC state machine is not initialized.")
        return self._chain_registry

    def _initialize_state_machine(self) -> None:
        self._transition(MCMCState.INIT, "initialize state machine")
        self._chain_registry = self._create_chain_registry()
        self._validate_chain_engines()
        self._pending_futures = set()
        self._future_to_sample = {}
        self._future_to_chain_id = {}
        self._inflight_chain_ids = set()
        self._ready_queue.clear()
        self._ready_set.clear()

        for cid in self._must_registry().ids():
            self._enqueue_chain(cid)

        self._state_machine_ready = True
        self._publish_metrics(event="run_start")
        patch = self._call_controller_hook("on_run_start", self._controller_context())
        self._apply_control_patch(patch, hook_name="on_run_start")

    def _validate_chain_engines(self) -> None:
        registry = self._must_registry()
        method = str(getattr(self, "method", "MCMC"))
        for chain in registry.all():
            validate_engine_contract(
                chain.engine,
                sampler_method=method,
                chain_id=int(chain.chain_id),
            )

    def _controller_context(self) -> Dict[str, Any]:
        registry = self._must_registry()
        return {
            "state": self._state.value,
            "method": getattr(self, "method", "MCMC"),
            "nchains": int(getattr(self, "_nchains", len(registry))),
            "niters": int(getattr(self, "_niters", 0)),
            "pending_futures": int(len(self._pending_futures)),
            "ready_queue": list(self._ready_queue),
            "cold_chains": registry.cold_chain_ids(),
            "chains": [chain.snapshot() for chain in registry.all()],
            "metrics_latest": (
                None if self.metrics_latest() is None else self.metrics_latest().__dict__.copy()
            ),
        }

    def _call_controller_hook(self, hook_name: str, *args, **kwargs):
        hook = getattr(self._controller, hook_name, None)
        if not callable(hook):
            return None
        try:
            return hook(*args, **kwargs)
        except Exception as exc:
            if getattr(self, "logger", None) is not None:
                self.logger.error(f"MCMC controller hook failed -> {hook_name} -> {exc}")
            return None

    def _apply_control_patch(self, patch: MCMCControlPatch | None, hook_name: str = "") -> bool:
        if patch is None:
            return False

        registry = self._must_registry()
        ok, reason = MCMCControlGuard.validate_patch(patch, nchains=len(registry))
        if not ok:
            if getattr(self, "logger", None) is not None:
                self.logger.warning(
                    format_two_column_log(
                        "MCMC controller patch rejected",
                        [("hook", hook_name or "unknown"), ("reason", reason)],
                    )
                )
            return False

        if patch.temperature_ladder is not None:
            ladder = [float(x) for x in patch.temperature_ladder]
            cold_idx = 0
            cold_val = min(ladder) if ladder else 1.0
            for ii, temp in enumerate(ladder):
                if temp == cold_val:
                    cold_idx = ii
                    break
            for ii, temp in enumerate(ladder):
                registry.set_temperature(ii, temp)
                registry.mark_cold(ii, ii == cold_idx)

        if patch.exchange_interval is not None and hasattr(self, "_exchange_interval"):
            self._exchange_interval = int(patch.exchange_interval)

        if patch.exchange_rule is not None and hasattr(self, "_exchange_rule"):
            self._exchange_rule = str(patch.exchange_rule)

        if patch.swap_pairing_mode is not None and hasattr(self, "_swap_pairing_mode"):
            self._swap_pairing_mode = str(patch.swap_pairing_mode)

        if patch.proposal_scales is not None:
            self._apply_proposal_scales(patch.proposal_scales)

        if patch.chain_priority_weight is not None:
            self._chain_priority_weight = [float(x) for x in patch.chain_priority_weight]

        if getattr(self, "logger", None) is not None:
            self.logger.info(
                format_two_column_log(
                    "MCMC controller patch applied",
                    [("hook", hook_name or "unknown")],
                )
            )
        return True

    def _enqueue_chain(self, chain_id: int) -> None:
        cid = int(chain_id)
        if cid in self._ready_set:
            return
        if cid in self._inflight_chain_ids:
            return
        chain = self._must_registry().get(cid)
        if chain.iter >= int(self._niters):
            return
        self._ready_queue.append(cid)
        self._ready_set.add(cid)

    def _pop_ready_chain(self) -> int | None:
        while self._ready_queue:
            cid = int(self._ready_queue.popleft())
            self._ready_set.discard(cid)
            if cid in self._inflight_chain_ids:
                continue
            chain = self._must_registry().get(cid)
            if chain.iter >= int(self._niters):
                continue
            return cid
        return None

    def _all_chains_finished(self) -> bool:
        for cid in self._must_registry().ids():
            if self._must_registry().get(cid).iter < int(self._niters):
                return False
        return True

    def _max_inflight(self) -> int:
        return int(getattr(self, "_nchains", 1) or 1)

    def _bounded_inflight(self, chain_count: int) -> int:
        chain_cap = max(1, int(chain_count))
        worker_cap = int(getattr(self, "max_workers", chain_cap) or chain_cap)
        if worker_cap <= 0:
            return chain_cap
        return max(1, min(chain_cap, worker_cap))

    def _next_save_dir(self) -> str | None:
        return None

    def _proposal_unit_point(self, chain: ChainRuntime):
        return next(chain.engine)

    def _propose_params_for_chain(self, chain: ChainRuntime):
        while True:
            proposal = self._proposal_unit_point(chain)
            params = self.map_point_into_distribution(proposal)
            if self._selectionexp:
                if self.evaluate_selection(self._selectionexp, params):
                    return params
                continue
            return params

    def _extract_logl(self, sample_info: Dict[str, Any]) -> float:
        observables = sample_info.get("observables", {})
        if "LogL" not in observables:
            raise KeyError("LogL not found in sample observables.")
        return float(observables["LogL"])

    def _update_chain_from_result(self, chain: ChainRuntime, sample_info: Dict[str, Any]) -> bool:
        logl_new = self._extract_logl(sample_info)
        beta = 1.0
        if chain.temperature > 0:
            beta = 1.0 / float(chain.temperature)

        accepted = None
        try:
            accepted = chain.engine.update(logl_new, beta=beta)
        except TypeError:
            accepted = chain.engine.update(logl_new)

        if accepted is None:
            accepted = True

        chain.iter = int(getattr(chain.engine, "iterations", chain.iter + 1))
        chain.last_logl = float(getattr(chain.engine, "last_loglikelihood", logl_new))
        if bool(accepted):
            chain.accepted += 1
        else:
            chain.rejected += 1
        chain.window_iter += 1

        chain.history.append_from_values(
            iter=chain.iter,
            state="accepted" if bool(accepted) else "rejected",
            proposal=getattr(chain.engine, "proposed_param", None),
            logl=chain.last_logl,
            accepted=bool(accepted),
            temperature=chain.temperature,
            meta={"chain_id": chain.chain_id},
        )
        return bool(accepted)

    def _on_before_submit(self, chain: ChainRuntime, sample: Sample) -> None:
        sample.info["chain_id"] = int(chain.chain_id)
        sample.info["chain_temperature"] = float(chain.temperature)
        sample.info["chain_is_cold"] = bool(chain.is_cold)

    def _on_after_update(self, chain: ChainRuntime, sample_info: Dict[str, Any], accepted: bool) -> None:
        _ = chain
        _ = sample_info
        _ = accepted

    def _should_exchange(self) -> bool:
        return False

    def _do_exchange(self) -> Dict[str, Any]:
        return {"attempted": 0, "accepted": 0}

    def _apply_proposal_scales(self, proposal_scales: Sequence[float] | float) -> None:
        _ = proposal_scales

    def _maybe_exchange_if_needed(self) -> None:
        if not self._should_exchange():
            return
        self._transition(MCMCState.EXCHANGE, "exchange gate")
        pre_patch = self._call_controller_hook("on_pre_exchange", self._controller_context())
        self._apply_control_patch(pre_patch, hook_name="on_pre_exchange")
        metrics = self._do_exchange()
        self._publish_metrics(event="post_exchange", extra=dict(metrics))
        post_patch = self._call_controller_hook(
            "on_post_exchange", self._controller_context(), dict(metrics)
        )
        self._apply_control_patch(post_patch, hook_name="on_post_exchange")

    @abstractmethod
    def _create_chain_registry(self) -> ChainRegistry:
        raise NotImplementedError

    def run_nested(self):
        if not self._state_machine_ready:
            self._initialize_state_machine()

        base_sample_cfg = self.info["sample"]
        self._transition(MCMCState.PROPOSE, "start run")

        while True:
            if self._all_chains_finished() and not self._pending_futures:
                break

            self._maybe_exchange_if_needed()
            pre_patch = self._call_controller_hook("on_pre_step", self._controller_context())
            self._apply_control_patch(pre_patch, hook_name="on_pre_step")

            self._transition(MCMCState.SUBMIT, "dispatch new proposals")
            while len(self._pending_futures) < self._max_inflight():
                cid = self._pop_ready_chain()
                if cid is None:
                    break
                chain = self._must_registry().get(cid)
                try:
                    params = self._propose_params_for_chain(chain)
                except StopIteration:
                    chain.iter = int(self._niters)
                    continue

                sample = Sample(params)
                save_dir = self._next_save_dir()
                sample_cfg = self.build_sample_config(base_sample_cfg, save_dir=save_dir)
                sample.set_config(sample_cfg)
                self._on_before_submit(chain, sample)
                try:
                    future = self.factory.submit_task(sample.info)
                except Exception:
                    self._on_sample_completed(sample.info)
                    sample.close()
                    raise
                self._pending_futures.add(future)
                self._future_to_sample[future] = sample
                self._future_to_chain_id[future] = cid
                self._inflight_chain_ids.add(cid)

            if not self._pending_futures:
                if self._all_chains_finished():
                    break
                if not self._ready_queue:
                    for chain in self._must_registry().all():
                        if chain.iter < int(self._niters):
                            self._enqueue_chain(chain.chain_id)
                    if not self._ready_queue and not self._pending_futures:
                        self._transition(MCMCState.FAILED, "deadlock detection")
                        raise RuntimeError(
                            "MCMC state machine deadlock: no pending futures and no ready chains."
                        )
                continue

            self._transition(MCMCState.WAIT, "wait for futures")
            done, _ = concurrent.futures.wait(
                self._pending_futures,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )

            self._transition(MCMCState.UPDATE, "apply completed results")
            for future in done:
                self._pending_futures.discard(future)
                sample = self._future_to_sample.pop(future, None)
                cid = self._future_to_chain_id.pop(future, None)
                if cid is not None:
                    self._inflight_chain_ids.discard(int(cid))

                try:
                    future.result()
                    if sample is None or cid is None:
                        continue
                    chain = self._must_registry().get(int(cid))
                    accepted = self._update_chain_from_result(chain, sample.info)
                    self._on_after_update(chain, sample.info, accepted)

                    outcome = {
                        "chain_id": int(cid),
                        "accepted": bool(accepted),
                        "iter": int(chain.iter),
                        "logl": chain.last_logl,
                    }
                    post_patch = self._call_controller_hook(
                        "on_post_step", self._controller_context(), outcome
                    )
                    self._apply_control_patch(post_patch, hook_name="on_post_step")

                    if chain.iter < int(self._niters):
                        self._enqueue_chain(cid)
                    self._publish_metrics(
                        event="post_step",
                        extra={
                            "chain_id": int(cid),
                            "accepted": bool(accepted),
                            "iter": int(chain.iter),
                        },
                    )
                except Exception as exc:
                    suuid = sample.uuid if sample is not None else "UNKNOWN"
                    self._transition(MCMCState.FAILED, "future failed")
                    self.logger.error(
                        format_two_column_log(
                            "[WorkerFactory] future exception consumed",
                            [("uuid", suuid), ("error", exc)],
                        )
                    )
                    raise
                finally:
                    if sample is not None:
                        self._on_sample_completed(sample.info)
                        sample.close()

        self._transition(MCMCState.TERMINATE, "all chains completed")
        self._publish_metrics(event="run_end")

    def _publish_metrics(self, event: str, extra: Dict[str, Any] | None = None) -> None:
        if self._chain_registry is None:
            return
        chains = self._chain_registry.all()
        acc_rates = []
        ess_proxies = []
        autocorr_proxies = []
        for chain in chains:
            denom = chain.accepted + chain.rejected
            if denom > 0:
                acc_rates.append(float(chain.accepted) / float(denom))
            ess_proxy, autocorr_proxy = self._estimate_chain_proxies(chain)
            if ess_proxy is not None:
                ess_proxies.append(float(ess_proxy))
            if autocorr_proxy is not None:
                autocorr_proxies.append(float(autocorr_proxy))
        acc_mean = float(sum(acc_rates) / len(acc_rates)) if acc_rates else 0.0
        ess_mean = float(sum(ess_proxies) / len(ess_proxies)) if ess_proxies else None
        autocorr_mean = (
            float(sum(autocorr_proxies) / len(autocorr_proxies)) if autocorr_proxies else None
        )

        swap_accept = None
        if isinstance(extra, dict) and "attempted" in extra and "accepted" in extra:
            attempted = float(extra.get("attempted", 0) or 0)
            accepted = float(extra.get("accepted", 0) or 0)
            swap_accept = (accepted / attempted) if attempted > 0 else 0.0

        frame = MCMCMetricsFrame(
            step=self._metrics_step,
            event=str(event),
            state=self._state.value,
            pending_futures=len(self._pending_futures),
            ready_queue=len(self._ready_queue),
            acceptance_rate_mean=acc_mean,
            ess_proxy_mean=ess_mean,
            autocorr_lag1_proxy_mean=autocorr_mean,
            swap_acceptance_rate=swap_accept,
            queue_depth_hint=len(self._pending_futures),
            meta=dict(extra or {}),
        )
        self._metrics_bus.publish(frame)
        self._metrics_step += 1

    def _estimate_chain_proxies(self, chain: ChainRuntime) -> tuple[float | None, float | None]:
        tail = chain.history.tail(32)
        logl_seq = [
            float(event.logl)
            for event in tail
            if getattr(event, "logl", None) is not None
        ]
        if len(logl_seq) < 3:
            return None, None

        n = len(logl_seq)
        mean = sum(logl_seq) / float(n)
        centered = [value - mean for value in logl_seq]
        var = sum(value * value for value in centered)
        if var <= 1e-16:
            return float(n), 0.0

        cov_lag1 = sum(centered[i] * centered[i - 1] for i in range(1, n))
        rho1 = cov_lag1 / var
        if rho1 > 0.99:
            rho1 = 0.99
        if rho1 < -0.99:
            rho1 = -0.99

        # Low-cost proxy for RL observation: ESS ~= N / (1 + 2 * max(rho1, 0)).
        ess = float(n) / float(1.0 + 2.0 * max(0.0, rho1))
        return ess, float(rho1)

    def _emit_chain_summary(self, label: str | None = None) -> None:
        if self._chain_registry is None or getattr(self, "logger", None) is None:
            return
        summary = self.chain_summary()
        if summary is None:
            return

        sampler_name = str(label or getattr(self, "method", "MCMC"))
        self.logger.warning(
            "{} Summary ->\n\tchains      -> {}\n\taccepted    -> {}\n\trejected    -> {}\n\tacc_rate    -> {:.4f}\n\tbest_logl   -> {}".format(
                sampler_name,
                int(summary["chains"]),
                int(summary["accepted"]),
                int(summary["rejected"]),
                float(summary["acc_rate"]),
                "N/A" if summary["best_logl"] is None else f"{float(summary['best_logl']):.6f}",
            )
        )

    def chain_summary(self) -> Dict[str, Any] | None:
        if self._chain_registry is None:
            return None
        chains = self._chain_registry.all()
        if not chains:
            return None

        accepted = sum(int(chain.accepted) for chain in chains)
        rejected = sum(int(chain.rejected) for chain in chains)
        total = accepted + rejected
        acc_rate = (float(accepted) / float(total)) if total > 0 else 0.0

        best_logl = None
        for chain in chains:
            logl = chain.last_logl
            if logl is None:
                continue
            if best_logl is None or float(logl) > float(best_logl):
                best_logl = float(logl)

        return {
            "method": str(getattr(self, "method", "MCMC")),
            "chains": int(len(chains)),
            "accepted": int(accepted),
            "rejected": int(rejected),
            "acc_rate": float(acc_rate),
            "best_logl": best_logl,
        }
