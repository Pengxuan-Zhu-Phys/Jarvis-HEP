"""Host feedback protocol, chain durability, bounded recovery and selection."""

from __future__ import annotations

import copy
import csv
import math
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from jarvis_hep_sampler import SamplerRegistry
from jarvis_hep_sampler.errors import PluginError, PluginResumeMismatch

from jarvishep2.feedback_return import (
    build_feedback_record,
    normalize_feedback_record,
    resolve_feedback_return,
)
from jarvishep2.likelihood import LogLikelihoodEvaluator
from jarvishep2.sampling.generational_plugin_sampler import GenerationalPluginSampler
from jarvishep2.sampling.plugin_sampler import preflight_plugin, create_plugin_sampler
from jarvishep2.runtime.worker_config import build_worker_config

EXAMPLES = Path(__file__).resolve().parents[2] / "Jarvis-HEP-Sampler" / "examples"


@pytest.fixture
def config(tmp_path, monkeypatch):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    SamplerRegistry().load(EXAMPLES / "SliceMCMC.py")
    data = yaml.safe_load((EXAMPLES / "SliceMCMC.eggbox.task.yaml").read_text())
    data.update(
        task_root=str(tmp_path),
        scan_name="slice",
        task_result_dir=str(tmp_path / "outputs/slice"),
    )
    data["Sampling"]["Bounds"].update(max_evaluations=13, burn_in=0)
    return data


def adapter(config):
    prepared = preflight_plugin(
        config, task_root=config["task_root"], scan_name=config["scan_name"]
    )
    value = create_plugin_sampler(prepared.spec)
    assert isinstance(value, GenerationalPluginSampler)
    value.attach_preflight(prepared)
    value.set_config(config)
    return value


def feedback(samples, *, failed=False, score=0.0):
    return [
        build_feedback_record(
            SimpleNamespace(
                uuid=s.uuid,
                status="Failed" if failed else "Completed",
                observables={"LogL": score},
                params={},
            ),
            terminal_version=1,
        )
        for s in samples
    ]


def absorb(value, *, score=0.0):
    samples = value.propose_generation()
    if samples is None:
        return None
    rows = [normalize_feedback_record(r) for r in feedback(samples, score=score)]
    value.absorb_generation(rows[::-1])
    value._generation += 1
    return samples


def test_global_budget_chain_export_and_restore(config):
    s = adapter(config)
    assert (
        build_worker_config(config, task_result_dir=config["task_result_dir"])[
            "plugin_feedback_version"
        ]
        == 1
    )
    absorb(s)
    checkpoint = s.export_runtime_state()
    expected_batch = s.propose_generation()
    # Replay from the saved barrier and remove chain suffix left by newer generations.
    s.absorb_generation(
        [normalize_feedback_record(r) for r in feedback(expected_batch)]
    )
    suffix = s._trajectory_path.read_bytes()
    restored = GenerationalPluginSampler(s._spec)
    restored.attach_preflight(
        preflight_plugin(config, resume=True, checkpoint_state=checkpoint)
    )
    restored.set_config(config)
    restored.import_runtime_state(checkpoint)
    assert len(restored._trajectory_path.read_bytes()) < len(suffix)
    replay = restored.propose_generation()
    assert [sample.uuid for sample in replay] == [
        sample.uuid for sample in expected_batch
    ]
    restored.absorb_generation(
        [normalize_feedback_record(r) for r in feedback(replay)][::-1]
    )
    while absorb(restored) is not None:
        pass
    rows = list(csv.DictReader(restored._trajectory_path.open()))
    assert restored._index == restored._accepted_index == 13
    assert len(rows) == restored._session.retained == 9
    assert len({(r["chain"], r["step"]) for r in rows}) == 9
    restored.assert_checkpoint_attribute_contract()


def test_selection_rejections_consume_budget_without_worker_tasks(config):
    config["Sampling"]["selection"] = "x < 0"
    s = adapter(config)
    while (batch := s.propose_generation()) is not None:
        assert batch == []
        s.absorb_generation([])
    assert s.summary()["proposals"] == s.summary()["selection_rejected"] == 13
    assert s.summary()["worker_evaluations"] == s.summary()["retained_draws"] == 0


@pytest.mark.parametrize("field,value", [("LogL", -math.inf), ("LogL", 0.0)])
def test_valid_terminal_including_negative_infinity(field, value):
    sample = SimpleNamespace(
        uuid="u", status="Completed", observables={field: value}, params={}
    )
    record = normalize_feedback_record(
        build_feedback_record(sample, terminal_version=1)
    )
    assert record["evaluation_outcome"] == "ok" and record["score_present"] is True
    assert record["logL"] == value


@pytest.mark.parametrize(
    "status,score",
    [
        ("Failed", 0),
        ("Completed", None),
        ("Completed", True),
        ("Completed", math.nan),
        ("Completed", math.inf),
    ],
)
def test_failed_or_invalid_score_gets_explicit_terminal_failure(status, score):
    record = build_feedback_record(
        SimpleNamespace(
            uuid="u", status=status, observables={"LogL": score}, params={}
        ),
        terminal_version=1,
    )
    assert record["evaluation_outcome"] == "failed"
    assert record["score_present"] is False and record["logL"] is None


def test_failed_worker_halts_and_keeps_previous_barrier(config):
    s = adapter(config)
    original = s.export_runtime_state()
    samples = s.propose_generation()
    with pytest.raises(PluginError, match="evaluation failed"):
        s.absorb_generation(feedback(samples, failed=True))
    assert s.checkpoint_runtime_state(safe=True) == original
    assert s._session.candidates == 0


@pytest.mark.parametrize(
    "expression,observables",
    [("nan", {}), ("oo", {}), ("z", {"z": None}), ("z", {"z": True})],
)
def test_strict_likelihood_never_converts_invalid_values_to_zero_density(
    expression, observables
):
    evaluator = LogLikelihoodEvaluator(
        [{"name": "LogL", "expression": expression}], strict=True
    )
    with pytest.raises((ValueError, KeyError)):
        evaluator.calculate({"observables": observables})


def test_strict_likelihood_accepts_explicit_minus_infinity():
    evaluator = LogLikelihoodEvaluator(
        [{"name": "LogL", "expression": "z"}], strict=True
    )
    assert evaluator.calculate({"observables": {"z": -math.inf}}) == -math.inf


def test_unsafe_checkpoint_cannot_advance_past_archive(config):
    s = adapter(config)
    original = s.checkpoint_runtime_state(safe=True)
    samples = s.propose_generation()
    s._submitted_uuids = [p.uuid for p in samples]
    s.absorb_generation([normalize_feedback_record(r) for r in feedback(samples)])
    assert s.checkpoint_runtime_state(safe=True) == original


def test_corrupt_chain_is_refused_in_preflight(config):
    s = adapter(config)
    absorb(s)
    absorb(s)
    checkpoint = s.export_runtime_state()
    s._trajectory_path.write_text("broken")
    with pytest.raises(PluginResumeMismatch, match="chain.csv"):
        preflight_plugin(config, resume=True, checkpoint_state=checkpoint)


@pytest.mark.parametrize("drift", ["budget", "likelihood", "opera", "selection"])
def test_scientific_drift_rejected_before_runtime(config, drift):
    s = adapter(config)
    checkpoint = s.export_runtime_state()
    changed = copy.deepcopy(config)
    if drift == "budget":
        changed["Sampling"]["Bounds"]["max_evaluations"] += 1
    elif drift == "likelihood":
        changed["Sampling"]["LogLikelihood"][0]["expression"] = "z/2"
    elif drift == "opera":
        changed["Operas"]["Modules"][0]["input"][0]["expression"] = "x/2"
    else:
        changed["Sampling"]["selection"] = "x > 0"
    with pytest.raises(PluginResumeMismatch, match="science settings"):
        preflight_plugin(changed, resume=True, checkpoint_state=checkpoint)


def test_feedback_requirements_are_merged_and_cannot_disable_logl(config):
    s = adapter(config)
    config["Sampling"]["FeedbackReturn"] = {"mode": "fields", "fields": ["z"]}
    assert resolve_feedback_return(config, sampler=s)["fields"] == ["z"]
    config["Sampling"]["FeedbackReturn"]["include_logl"] = False
    with pytest.raises(PluginError, match="include_logl"):
        preflight_plugin(config)


def test_equivalent_feedback_duplicates_are_idempotent_and_conflicts_fail(
    config, monkeypatch
):
    from jarvishep2.sampling.feedback_sampler import FeedbackSampler

    s = adapter(config)
    samples = s.propose_generation()
    record = normalize_feedback_record(feedback(samples)[0])
    monkeypatch.setattr(
        FeedbackSampler, "_await_feedback_round", lambda *a, **kw: (record, "", None)
    )
    assert s._await_feedback_round()[0] == record
    assert s._await_feedback_round()[0] == record
    record = dict(record, logL=1.0)
    with pytest.raises(PluginError, match="conflicting duplicate"):
        s._await_feedback_round()


def test_missing_worker_capability_is_rejected(config):
    s = adapter(config)
    samples = s.propose_generation()
    rows = feedback(samples)
    rows[0].pop("plugin_feedback_version")
    with pytest.raises(PluginError, match="terminal feedback"):
        s.absorb_generation(rows)
    assert s._session.candidates == 0


def test_archive_lag_cannot_publish_new_checkpoint(config, monkeypatch):
    from jarvishep2.sampling import generational_plugin_sampler as implementation

    s = adapter(config)
    original = s.export_runtime_state()
    samples = s.propose_generation()
    s._submitted_uuids = [sample.uuid for sample in samples]
    s.absorb_generation([normalize_feedback_record(row) for row in feedback(samples)])
    s.redis = SimpleNamespace(get_archived_index_prefix=lambda scan: 0)
    times = iter([0.0, 31.0])
    monkeypatch.setattr(implementation.time, "monotonic", lambda: next(times))
    with pytest.raises(TimeoutError, match="archive did not confirm"):
        s.checkpoint_at_barrier(final=True)
    assert s.checkpoint_runtime_state(safe=True) == original


def test_chain_write_failure_preserves_previous_durable_state(config, monkeypatch):
    s = adapter(config)
    absorb(s)
    original = s.checkpoint_runtime_state(safe=True)
    samples = s.propose_generation()
    real_open = Path.open

    def broken_open(path, mode="r", *args, **kwargs):
        if path == s._trajectory_path and mode == "ab":
            raise OSError("chain disk full")
        return real_open(path, mode, *args, **kwargs)

    monkeypatch.setattr(Path, "open", broken_open)
    with pytest.raises(OSError, match="disk full"):
        s.absorb_generation(
            [normalize_feedback_record(row) for row in feedback(samples)]
        )
    assert s.checkpoint_runtime_state(safe=True) == original
