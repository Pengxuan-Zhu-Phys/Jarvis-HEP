#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
from typing import Any, Dict

from jarvishep.log_kv import format_two_column_log
from jarvishep.sample import Sample

from .chain_runtime import ChainRuntime
from .state_machine_base import MCMCState, MCMCStateMachineBase


class MCMCMultiStageStateMachineBase(MCMCStateMachineBase):
    """Parallel state-machine runtime for samplers with delayed-rejection stages."""

    def __init__(self) -> None:
        MCMCStateMachineBase.__init__(self)
        self._future_to_stage_idx: Dict[Any, int] = {}
        self._chain_stage_idx: Dict[int, int] = {}

    def _initialize_state_machine(self) -> None:
        MCMCStateMachineBase._initialize_state_machine(self)
        self._future_to_stage_idx = {}
        self._chain_stage_idx = {}

    def _initial_stage_index(self, chain: ChainRuntime) -> int:
        _ = chain
        return 0

    def _proposal_unit_point_for_stage(self, chain: ChainRuntime, stage_index: int):
        _ = stage_index
        return self._proposal_unit_point(chain)

    def _propose_params_for_stage(self, chain: ChainRuntime, stage_index: int):
        while True:
            proposal = self._proposal_unit_point_for_stage(chain, stage_index)
            params = self.map_point_into_distribution(proposal)
            if self._selectionexp:
                if self.evaluate_selection(self._selectionexp, params):
                    return params
                continue
            return params

    def _consume_stage_result(
        self,
        chain: ChainRuntime,
        sample_info: Dict[str, Any],
        stage_index: int,
    ) -> Dict[str, Any]:
        accepted = self._update_chain_from_result(chain, sample_info)
        return {
            "iteration_done": True,
            "accepted": bool(accepted),
            "next_stage": None,
            "logl": chain.last_logl,
            "stage_attempts": int(stage_index) + 1,
        }

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
                stage_idx = int(self._chain_stage_idx.get(cid, self._initial_stage_index(chain)))
                try:
                    params = self._propose_params_for_stage(chain, stage_idx)
                except StopIteration:
                    chain.iter = int(self._niters)
                    self._chain_stage_idx.pop(cid, None)
                    continue

                sample = Sample(params)
                save_dir = self._next_save_dir()
                sample_cfg = self.build_sample_config(base_sample_cfg, save_dir=save_dir)
                sample.set_config(sample_cfg)
                self._on_before_submit(chain, sample)
                sample.info["chain_stage"] = int(stage_idx)
                try:
                    future = self.factory.submit_task(sample.info)
                except Exception:
                    self._on_sample_completed(sample.info)
                    sample.close()
                    raise
                self._pending_futures.add(future)
                self._future_to_sample[future] = sample
                self._future_to_chain_id[future] = cid
                self._future_to_stage_idx[future] = int(stage_idx)
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
                            "MCMC multi-stage deadlock: no pending futures and no ready chains."
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
                stage_idx = int(self._future_to_stage_idx.pop(future, 0))
                if cid is not None:
                    self._inflight_chain_ids.discard(int(cid))

                try:
                    future.result()
                    if sample is None or cid is None:
                        continue

                    chain = self._must_registry().get(int(cid))
                    outcome = self._consume_stage_result(chain, sample.info, stage_idx)
                    iteration_done = bool(outcome.get("iteration_done", True))

                    if iteration_done:
                        accepted = bool(outcome.get("accepted", False))
                        chain.iter = int(getattr(chain.engine, "iterations", chain.iter + 1))
                        if getattr(chain.engine, "last_loglikelihood", None) is not None:
                            chain.last_logl = float(chain.engine.last_loglikelihood)
                        else:
                            chain.last_logl = float(
                                outcome.get("logl", self._extract_logl(sample.info))
                            )
                        if accepted:
                            chain.accepted += 1
                        else:
                            chain.rejected += 1
                        chain.window_iter += 1

                        stage_attempts = int(outcome.get("stage_attempts", stage_idx + 1))
                        chain.history.append_from_values(
                            iter=chain.iter,
                            state="accepted" if accepted else "rejected",
                            proposal=getattr(chain.engine, "proposed_param", None),
                            logl=chain.last_logl,
                            accepted=accepted,
                            temperature=chain.temperature,
                            meta={
                                "chain_id": int(chain.chain_id),
                                "stage_attempts": stage_attempts,
                                "final_stage": int(stage_idx),
                            },
                        )

                        self._on_after_update(chain, sample.info, accepted)
                        self._chain_stage_idx.pop(int(cid), None)

                        post_outcome = {
                            "chain_id": int(cid),
                            "accepted": bool(accepted),
                            "iter": int(chain.iter),
                            "logl": chain.last_logl,
                            "stage_attempts": stage_attempts,
                        }
                        post_patch = self._call_controller_hook(
                            "on_post_step", self._controller_context(), post_outcome
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
                                "stage_attempts": stage_attempts,
                            },
                        )
                    else:
                        next_stage = outcome.get("next_stage", None)
                        if next_stage is None:
                            raise RuntimeError(
                                f"Multi-stage sampler requires next_stage after rejection, chain={cid} stage={stage_idx}"
                            )
                        self._chain_stage_idx[int(cid)] = int(next_stage)
                        self._enqueue_chain(int(cid))
                        self._publish_metrics(
                            event="post_stage_reject",
                            extra={
                                "chain_id": int(cid),
                                "iter": int(chain.iter),
                                "stage": int(stage_idx),
                                "next_stage": int(next_stage),
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
