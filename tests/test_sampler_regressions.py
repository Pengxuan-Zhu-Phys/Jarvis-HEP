"""Sampler audit: prior correctness, bounded recovery state and nested runtime."""
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


@pytest.mark.parametrize("cls", [MCMCChain, AMMCMCChain, EnsembleChain, DEMCMCChain, DRAMChain])
def test_invalid_initial_likelihood_never_becomes_chain_state(cls):
    c = engine(cls)
    initial = c.param.copy()
    next(c)
    assert not c.update(-np.inf)
    np.testing.assert_array_equal(c.param, initial)
    assert c.last_loglikelihood is None
    next(c)
    assert c.update(0.)
    assert c.last_loglikelihood == 0.


@pytest.mark.parametrize("cls", [MCMCChain, AMMCMCChain, EnsembleChain, DEMCMCChain])
def test_uniform_target_keeps_boundary_probability(cls):
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


@pytest.mark.parametrize("method", ["MCMC", "AMMCMC", "DRAM", "EnsembleMCMC", "DEMCMC", "PTMCMC", "PTEnsemble"])
def test_selection_rejection_stays_local_and_finishes(tmp_path, method):
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


@pytest.mark.parametrize("cls", [DynestySampler, MultiNestSampler])
def test_nested_queue_uses_workers_not_transport_chunk(tmp_path, cls):
    sampler = cls()
    sampler.set_config(config(tmp_path, sampler.method))
    kwargs = sampler._build_constructor_kwargs(pool=Mock(), rstate=np.random.default_rng(0))
    assert kwargs["queue_size"] == 4


@pytest.mark.parametrize("cls", [DynestySampler, MultiNestSampler])
def test_completed_nested_checkpoint_does_not_restart(tmp_path, cls):
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
    assert len(queue.push_many_tasks.call_args.args[0]) == 1
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


def test_nested_timed_save_does_not_reserialize_in_callback(tmp_path, monkeypatch):
    sampler = DynestySampler()
    sampler.set_config(config(tmp_path, "Dynesty"))
    sampler._sampler = SavedEngine()
    pickle_native = Mock(return_value=b"engine")
    monkeypatch.setattr(sampler, "_pickle_native_sampler", pickle_native)
    side_save = Mock()
    monkeypatch.setattr(sampler, "_save_engine_side_file", side_save)
    states = []
    sampler._save_checkpoint_callback = lambda **kw: states.append(sampler.checkpoint_runtime_state(safe=True))
    sampler._install_dynesty_save_hook()
    sampler._sampler.save("unused")
    assert sampler._sampler.saves == 1
    pickle_native.assert_called_once()
    side_save.assert_not_called()
    assert states[0]["native_sampler_blob"] == b"engine"
    assert not sampler._native_snapshot_ready


def test_nested_disabled_checkpoint_never_writes_engine(tmp_path, monkeypatch):
    sampler = DynestySampler()
    sampler.set_config(config(tmp_path, "Dynesty"))
    sampler._sampler = SavedEngine()
    monkeypatch.setattr(sampler, "_pickle_native_sampler", lambda: b"engine")
    side_save = Mock()
    monkeypatch.setattr(sampler, "_save_engine_side_file", side_save)
    sampler.export_runtime_state()
    side_save.assert_not_called()


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


def test_slow_mcmc_checkpoint_does_not_trigger_save_storm(tmp_path, monkeypatch):
    import jarvishep2.sampling.mcmc_sampler as module
    sampler = Distributor.set_method("EnsembleMCMC")
    sampler.set_config(config(tmp_path, sampler.method))
    now = [100.]
    monkeypatch.setattr(module.time, "monotonic", lambda: now[0])
    saves = []
    def save(**kwargs):
        saves.append(now[0])
        now[0] += 60.  # Save takes longer than the configured 30s interval.
        return True
    monkeypatch.setattr(sampler, "persist_runtime_checkpoint", save)
    sampler._sampling_checkpoint_interval_sec = 30.
    assert sampler._checkpoint_at_sampling_barrier(reason="test")
    now[0] += 1.
    assert not sampler._checkpoint_at_sampling_barrier(reason="test")
    now[0] += 30.
    assert sampler._checkpoint_at_sampling_barrier(reason="test")
    assert len(saves) == 2


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
