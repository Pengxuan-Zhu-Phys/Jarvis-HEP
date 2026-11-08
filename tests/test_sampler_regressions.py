"""Sampler audit: prior correctness, bounded recovery state and nested runtime."""
import os
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from jarvishep2.distributor import Distributor
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.sampling.Source.MCMC.mcmc_chain import MCMCChain
from jarvishep2.sampling.Source.MCMC.ammcmc_chain import AMMCMCChain
from jarvishep2.sampling.Source.MCMC.dram_chain import DRAMChain
from jarvishep2.sampling.Source.MCMC.engine_ensemble import EnsembleChain
from jarvishep2.sampling.Source.MCMC.engine_demcmc import DEMCMCChain
from jarvishep2.sampling.dynesty_sampler import DynestySampler, NestedRedisLogL
from jarvishep2.sampling.multinest_sampler import MultiNestSampler
from jarvishep2.sampling.redis_evaluation_pool import RedisEvaluationPool


def config(tmp_path, method, **bounds):
    return {
        "task_root": str(tmp_path), "task_result_dir": str(tmp_path),
        "Runtime": {"mode": "redis", "workers": 4, "batch_size": 1,
                    "checkpoint": {"enabled": False}},
        "Sampling": {"Method": method,
                     "Bounds": {"num_chains": 4, "num_iters": 12, "seed": 7, **bounds},
                     "Variables": [{"name": "x", "distribution": {
                         "type": "Flat", "parameters": {"min": 0., "max": 1.}}}]},
    }


def engine(cls, seed=10):
    kw = {"rng": np.random.default_rng(seed)}
    if cls in (EnsembleChain, DEMCMCChain):
        kw.update(chain_id=0, population_getter=lambda _: [np.array([.2]), np.array([.8])])
    if cls is DEMCMCChain:
        kw.update(de_gamma=.5, de_noise=.1)
    if cls in (AMMCMCChain, DRAMChain):
        kw["adapt_enabled"] = False
    return cls([.5], .25, 60000, **kw)


def test_invalid_initial_likelihood_never_becomes_chain_state(subtests):
    for cls in [MCMCChain, AMMCMCChain, EnsembleChain, DEMCMCChain, DRAMChain]:
        with subtests.test(cls=cls):
            c = engine(cls)
            initial = c.param.copy()
            next(c)
            assert not c.update(-np.inf)
            np.testing.assert_array_equal(c.param, initial)
            assert c.last_loglikelihood is None
            next(c)
            assert c.update(0.)
            assert c.last_loglikelihood == 0.


def test_uniform_target_keeps_boundary_probability(subtests):
    for cls in [MCMCChain, AMMCMCChain, EnsembleChain, DEMCMCChain]:
        with subtests.test(cls=cls):
            c = engine(cls)
            values = []
            for i in range(40000):
                next(c)
                c.update(0.)  # Engine must reject outside points, even with finite logL.
                if i >= 1000:
                    values.append(c.param[0])
            values = np.asarray(values)
            edge_fraction = np.mean((values < .1) | (values > .9))
            assert .18 < edge_fraction < .22, (cls.__name__, edge_fraction)
            assert abs(np.var(values) - 1 / 12) < .006


def test_selection_rejection_stays_local_and_finishes(subtests, tmp_path_factory):
    for method in ["MCMC", "AMMCMC", "DRAM", "EnsembleMCMC", "DEMCMC", "PTMCMC", "PTEnsemble"]:
        with subtests.test(method=method):
            tmp_path = tmp_path_factory.mktemp("case")
            cfg = config(tmp_path, method)
            cfg["Sampling"]["selection"] = "x < 0"
            sampler = Distributor.set_method(method)
            sampler.set_config(cfg)
            sampler.set_redis(make_fakeredis_queue())
            sampler._submit_group = Mock(side_effect=AssertionError("invalid point reached Redis"))
            sampler.run_adaptive(timeout=1)
            assert sampler._all_finished()
            assert sampler._total_accepted == 0
            assert sampler._total_failed == 0
            assert not sampler._uuid_to_meta
            sampler._submit_group.assert_not_called()


def test_failed_uuid_history_is_bounded_and_count_survives_resume(tmp_path):
    cfg = config(tmp_path, "MCMC", num_iters=800)
    sampler = Distributor.set_method("MCMC")
    sampler.set_config(cfg)
    for _ in range(150):
        batch = sampler.propose_generation()
        sampler.absorb_generation([{"uuid": s.uuid, "logL": -np.inf} for s in batch])
    assert sampler._total_failed == 600
    assert len(sampler._failed_uuids) == 256
    state = sampler.export_runtime_state()
    restored = Distributor.set_method("MCMC")
    restored.set_config(cfg)
    restored.import_runtime_state(state)
    assert restored.summary()["failed_samples"] == 600
    assert len(restored._failed_uuids) == 256
    restored.assert_checkpoint_attribute_contract()
    sampler._close_history_csv()
    restored._close_history_csv()


def test_adaptive_history_includes_repeated_rejected_states():
    c = engine(AMMCMCChain)
    next(c)
    c.update(0.)
    for _ in range(10):
        next(c)
        c.update(-np.inf)
    assert len(c._history) == 11
    np.testing.assert_array_equal(c._history[0], c._history[-1])


def test_dram_refuses_unimplemented_higher_order_acceptance():
    with pytest.raises(ValueError, match="dr_steps"):
        DRAMChain([.5], .1, 10, dr_steps=3)


def test_nested_queue_uses_workers_not_transport_chunk(subtests, tmp_path_factory):
    for cls in [DynestySampler, MultiNestSampler]:
        with subtests.test(cls=cls):
            tmp_path = tmp_path_factory.mktemp("case")
            sampler = cls()
            sampler.set_config(config(tmp_path, sampler.method))
            kwargs = sampler._build_constructor_kwargs(pool=Mock(), rstate=np.random.default_rng(0))
            assert kwargs["queue_size"] == 4


def test_completed_nested_checkpoint_does_not_restart(subtests, tmp_path_factory):
    for cls in [DynestySampler, MultiNestSampler]:
        with subtests.test(cls=cls):
            tmp_path = tmp_path_factory.mktemp("case")
            sampler = cls()
            sampler.set_config(config(tmp_path, sampler.method))
            sampler.import_runtime_state({"finished": True, "summary": {"ncall": 123}})
            sampler.set_redis(make_fakeredis_queue())
            assert sampler.run_adaptive() == 123
            assert sampler._sampler is None


def test_nested_selection_applies_to_batch_and_direct_calls(tmp_path):
    sampler = DynestySampler()
    cfg = config(tmp_path, "Dynesty")
    cfg["Sampling"]["selection"] = "x < 0.5"
    sampler.set_config(cfg)
    queue = Mock()
    pool = RedisEvaluationPool(queue, build_sample=sampler._build_sample_for_pool,
                               accepts_payload=sampler._accepts_pool_payload)
    def wait(uuids, futs):
        pool._forget_waiters(uuids)
        return [42.] * len(uuids)
    pool._wait_feedback = wait
    values = pool._redis_batch_logl([np.array([.8]), np.array([.2]), np.array([.9])])
    assert [v.val for v in values] == [-1e300, 42., -1e300]
    pushed = queue.push_many_tasks.call_args.args[0]
    # Selection rejects take no sample_index, so DATABASE indices stay dense.
    assert [task["sample_index"] for task in pushed] == [0]
    assert pool._call_index == 1
    queue.reset_mock()
    assert NestedRedisLogL(pool)(np.array([.8])) == -1e300
    queue.push_many_tasks.assert_not_called()


def test_waiter_registration_is_atomic_on_duplicate():
    pool = RedisEvaluationPool(Mock(), build_sample=Mock())
    pool._register_waiters(["existing"])
    with pytest.raises(ValueError, match="duplicate"):
        pool._register_waiters(["new", "existing"])
    assert set(pool._waiters) == {"existing"}


class SavedEngine:
    def __init__(self):
        self.saves = 0

    def save(self, filename):
        self.saves += 1


def test_nested_timed_save_submits_when_the_writer_is_idle(tmp_path, monkeypatch):
    sampler = DynestySampler()
    sampler.set_config(config(tmp_path, "Dynesty"))
    sampler._sampler = SavedEngine()
    sampler._pool_call_index = 8
    pickle_native = Mock(return_value=b"engine")
    monkeypatch.setattr(sampler, "_pickle_native_sampler", pickle_native)

    class Writer:
        def __init__(self):
            self.idle = False
            self.payloads = []

        def available(self):
            return self.idle

        def wait(self, timeout=10.0):
            return True

        def submit(self, payload):
            self.payloads.append(payload)
            return True

    writer = Writer()
    sampler._pickle_writer = writer
    sampler._install_dynesty_save_hook()
    sampler._sampler.save("unused")
    # A busy checkpoint process returns before the engine is pickled.
    assert sampler._sampler.saves == 0
    pickle_native.assert_not_called()
    writer.idle = True
    sampler._sampler.save("unused")
    pickle_native.assert_called_once()
    assert writer.payloads[0]["sampler_state"]["checkpoint_index"] == 8
    assert writer.payloads[0]["sampler_state"]["native_sampler_blob"] == b"engine"
    assert not os.path.exists(sampler._engine_path())


def test_nested_disabled_checkpoint_never_writes_engine(tmp_path, monkeypatch):
    sampler = DynestySampler()
    sampler.set_config(config(tmp_path, "Dynesty"))
    sampler.config.setdefault("EnvReqs", {}).setdefault("V2", {}).setdefault("checkpoint", {})["enabled"] = False
    sampler._sampler = SavedEngine()
    monkeypatch.setattr(sampler, "_pickle_native_sampler", lambda *a: b"engine")
    sampler.configure_checkpoint(checkpoint_file=str(tmp_path / "runtime.pickle"))
    assert sampler.save_checkpoint(force=True) is False
    sampler.export_runtime_state()
    assert not os.path.exists(sampler._engine_path())


def test_nested_resume_reattaches_current_parallelism(tmp_path):
    sampler = DynestySampler()
    sampler.set_config(config(tmp_path, "Dynesty"))
    sampler.import_runtime_state({"constructor_kwargs": {"queue_size": 1}})
    layer = SimpleNamespace(queue_size=1, use_pool_evolve=True, loglikelihood=None)
    pool = Mock()
    sampler._reattach_native_runtime(layer, pool=pool, bridge=NestedRedisLogL(pool))
    assert layer.queue_size == 4
    assert layer.mapper == pool.map


def test_native_legacy_adaptive_history_is_bounded_on_restore():
    import pickle
    from jarvishep2.sampling.mcmc_sampler import _unpickle_mcmc_engine
    c = engine(AMMCMCChain)
    c._history = [np.array([.5]) for _ in range(10000)]
    restored = _unpickle_mcmc_engine(pickle.dumps(c), Mock())
    assert restored._history.maxlen == restored._history_maxlen()
    assert len(restored._history) == restored._history_maxlen()


def test_dram_uniform_target_with_delayed_boundary_rejection():
    c = engine(DRAMChain)
    values = []
    followups = 0
    for i in range(20000):
        c.propose_stage(0)
        out = c.consume_stage_result(0, 0.)
        if not out["iteration_done"]:
            followups += 1
            c.propose_stage(1)
            c.consume_stage_result(1, 0.)
        if i >= 1000:
            values.append(c.param[0])
    values = np.asarray(values)
    assert followups > 500
    assert .18 < np.mean((values < .1) | (values > .9)) < .22


def test_dram_reverse_acceptance_one_forces_stage_two_rejection():
    c = engine(DRAMChain)
    c.param = np.array([.5])
    c.last_loglikelihood = 0.
    c._stage_proposals = {0: np.array([.6]), 1: np.array([.4])}
    c._stage_logl = {0: -1.}
    c._stage_alpha = {0: np.exp(-1.)}
    assert c._accept_prob_stage2(-2., 1.) == 0.


def test_nested_pool_threads_are_independent_of_batch_size():
    import threading
    from collections import namedtuple
    SamplerArgument = namedtuple("SamplerArgument", "value")
    barrier = threading.Barrier(4, timeout=2)
    pool = RedisEvaluationPool(Mock(), build_sample=Mock(), batch_size=1, njobs=4)
    def walk(arg):
        barrier.wait()
        return arg.value
    assert pool.map(walk, [SamplerArgument(i) for i in range(4)]) == list(range(4))


def test_mcmc_checkpoint_only_captures_when_requested_and_idle(tmp_path, monkeypatch):
    sampler = Distributor.set_method("EnsembleMCMC")
    sampler.set_config(config(tmp_path, sampler.method))
    writer = Mock()
    sampler.configure_checkpoint(checkpoint_file=str(tmp_path / "runtime.pickle"), writer=writer)
    snapshot = Mock(wraps=sampler.checkpoint_state)
    monkeypatch.setattr(sampler, "checkpoint_state", snapshot)
    writer.reset_mock()
    writer.requested.return_value = False
    assert not sampler.save_checkpoint()
    writer.requested.return_value = True
    writer.available.return_value = False
    assert not sampler.save_checkpoint()
    snapshot.assert_not_called()
    writer.wait.assert_not_called()
    writer.available.return_value = True
    assert sampler.save_checkpoint()
    snapshot.assert_called_once()
    writer.submit.assert_called_once()


def mcmc_with_feedback(cfg, score=None):
    """Run the real sampler loops with immediate deterministic Worker feedback."""
    sampler = Distributor.set_method(cfg["Sampling"]["Method"])
    sampler.set_config(cfg)
    queue = make_fakeredis_queue()
    sampler.set_redis(queue)
    proposals = {}

    def submit(task):
        proposals[task["uuid"]] = task["u_coords"]
        logl = score(task) if score else -float(np.sum((np.asarray(task["u_coords"]) - 0.3) ** 2))
        queue.publish_feedback({"uuid": task["uuid"], "logL": logl},
                               queue=task.get("feedback_queue"))

    queue.push_task = submit
    queue.push_many_tasks = lambda tasks: [submit(task) for task in tasks]
    return sampler, proposals


@pytest.mark.parametrize("method", [
    "MCMC", "ToyMCMC", "AMMCMC", "DRAM", "DEMCMC", "EnsembleMCMC", "PTMCMC", "PTEnsemble",
])
def test_mcmc_checkpoint_preserves_chain_trajectories(tmp_path, method):
    extra = {"adapt_start_iter": 2, "adapt_window": 3} if method in {"AMMCMC", "DRAM"} else {}
    full, full_proposals = mcmc_with_feedback(config(tmp_path / "full", method, **extra))
    full.run_adaptive(timeout=2)

    cfg = config(tmp_path / "resumed", method, **extra)
    interrupted, before = mcmc_with_feedback(cfg)
    writer = Mock()
    interrupted.configure_checkpoint(checkpoint_file=str(tmp_path / "runtime.pickle"), writer=writer)

    def stop_at_checkpoint(payload):
        if payload["sampler_state"]["progress"]["generation"] >= 5:
            raise KeyboardInterrupt
        return True

    writer.submit.side_effect = stop_at_checkpoint
    with pytest.raises(KeyboardInterrupt):
        interrupted.run_adaptive(timeout=2)
    state = writer.submit.call_args.args[0]["sampler_state"]
    assert 5 <= state["progress"]["generation"] < 12
    assert bool(state["pending"]) == interrupted._uses_async_chain_pipeline()

    resumed, after = mcmc_with_feedback(cfg)
    resumed.restore_checkpoint(state)
    assert set(resumed.repropose_unfinished()) == set(state["pending"])
    assert after == {uuid: before[uuid] for uuid in state["pending"]}
    resumed.run_adaptive(timeout=2)
    assert before | after == full_proposals
    actual = resumed.checkpoint_state()
    expected = full.checkpoint_state()
    for checkpoint in (actual, expected):
        checkpoint["progress"]["summary"].pop("diagnostics", None)
        for chain in checkpoint["chains"]:
            for event in chain["history"]["_events"]:
                event.timestamp_utc = ""
    assert actual == expected


@pytest.mark.parametrize("already_proposed", [False, True])
def test_dram_checkpoint_between_delayed_rejection_stages(tmp_path, already_proposed):
    cfg = config(tmp_path, "DRAM", proposal_scale=0.01)

    def score(task):
        return -1000.0 if task["step"] == 1 and task["stage"] == 0 else 0.0

    sampler, before = mcmc_with_feedback(cfg, score)
    for _ in range(2):
        sampler._submit_sample_batch(sampler.propose_generation())
        sampler.absorb_generation(sampler.wait_for_generation(timeout=2, queues=sampler._pending_feedback_queues()))
    assert all(chain.open_stage == 1 for chain in sampler._ensure_registry().all())
    if already_proposed:
        sampler._submit_sample_batch(sampler._propose_for_stages(stages="followup"))
    state = sampler.checkpoint_state()
    resumed, after = mcmc_with_feedback(cfg, score)
    resumed.restore_checkpoint(state)
    replay = resumed.repropose_unfinished()
    if already_proposed:
        assert len(replay) == 4
        assert after == {uuid: before[uuid] for uuid in replay}
    else:
        assert not replay
        expected = {
            s.uuid: s.u_coords.tolist() for s in sampler._propose_for_stages(stages="followup")
        }
        assert {
            s.uuid: s.u_coords.tolist() for s in resumed._propose_for_stages(stages="followup")
        } == expected


def test_mcmc_legacy_checkpoint_import_and_runtime_precedence(tmp_path):
    from pathlib import Path

    from jarvishep2.runtime.core import Jarvis2Core
    from jarvishep2.sampling.runtime_checkpoint import build_payload, save_checkpoint
    from jarvishep2.sampling.sample_checkpoint import atomic_pickle

    cfg = config(tmp_path, "ToyMCMC")
    cfg["Runtime"]["checkpoint"]["enabled"] = True
    core = Jarvis2Core(cfg)
    sampler, proposals = mcmc_with_feedback(core.config)
    sampler._top_up_ready_chains()
    legacy_state = sampler.export_runtime_state()
    path = Path(core.checkpoint_file())
    legacy = path.with_name("state.pkl")
    legacy_payload = build_payload(
        run_spec={"normalized_config": core.config},
        sampler_state=legacy_state,
    )
    save_checkpoint(str(legacy), legacy_payload)
    original = legacy.read_bytes()
    core.config["_unit_execution"] = {"execution_identity": "test-units"}
    core._get_resume()._preload_sample_checkpoint(resume=False, fresh=False)
    assert core._resume_policy == "fresh"
    core.config.pop("_unit_execution")

    unit_legacy = dict(legacy_payload, integrity={**legacy_payload["integrity"], "unit_execution": "test-units"})
    save_checkpoint(str(legacy), unit_legacy)
    with pytest.raises(ValueError, match="unit-aware MCMC resume requires runtime.pickle"):
        core.prepare_resume(resume=True)
    legacy.write_bytes(original)

    core.prepare_resume(resume=True)
    resumed, replay = mcmc_with_feedback(core.config)
    core.set_sampler(resumed)
    assert set(resumed.repropose_unfinished()) == set(proposals)
    assert replay == proposals

    writer = Mock()
    resumed.configure_checkpoint(checkpoint_file=str(path), writer=writer)
    payload = writer.submit.call_args.args[0]
    atomic_pickle(str(path), payload)
    assert legacy.read_bytes() == original

    invalid = dict(payload)
    invalid["sampler_state"] = {"pending": {}}
    atomic_pickle(str(path), invalid)
    core.prepare_resume(resume=True)
    fresh, _ = mcmc_with_feedback(core.config)
    with pytest.raises(ValueError, match="incomplete checkpoint"):
        core.set_sampler(fresh)
    assert path.exists()
    assert legacy.read_bytes() == original


def test_nested_timer_restarts_after_slow_save(monkeypatch):
    from jarvishep2.sampling.Source.Dynesty.py.dynesty import utils
    now = [100.]
    monkeypatch.setattr(utils.time, "time", lambda: now[0])
    timer = utils.DelayTimer(30.)
    now[0] += 31.
    assert timer.is_time()
    now[0] += 60.
    timer.reset()
    now[0] += 1.
    assert not timer.is_time()
