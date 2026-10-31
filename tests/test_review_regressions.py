"""Reproductions of runtime and configuration review findings."""
from __future__ import annotations

import json
from types import SimpleNamespace
from unittest import mock

import pytest

from jarvishep2.calculators.profile1d import Profile1DProfiler
from jarvishep2.command_parser import CommandParser
from jarvishep2.core import Jarvis2Core
from jarvishep2.factory import TaskFactory
from jarvishep2.sample import Sample, ExecutionStep
from jarvishep2.sampling.variables import Variable
from jarvishep2.task_validation import validate_task_config
from jarvishep2.workflow import resolve_module_layers


def test_nuisance_failure_remains_failed_after_pipeline():
    from jarvishep2.runtime.sample_executor import SampleExecutor
    sample = Sample(uuid="failed-nuisance", info={"store_samples": False})
    sample.execution_plan = [ExecutionStep(type="nuisance_optimize", name="nuisance", layer=0)]
    worker = SimpleNamespace(worker_config={}, _mapper=None, _scheduler=None,
                             _run_layer=lambda layer, point: point.set_status("Failed"))
    SampleExecutor(worker).process(sample)
    assert sample.status == "Failed"


@pytest.mark.parametrize("block", [{"method": "typo"}, {"target_mode": "typo"}])
def test_nuisance_unknown_modes_fail_before_evaluation(block):
    nuisance = {"variables": [{"name": "z", "distribution": {"parameters": {"min": 0, "max": 1}}}], **block}
    with pytest.raises(ValueError):
        Profile1DProfiler.from_config({"Sampling": {"Nuisance": nuisance}})
    issues = validate_task_config({"Sampling": {"Nuisance": nuisance}}).errors()
    assert any("Nuisance" in issue.path for issue in issues)


@pytest.mark.parametrize("lo", [0, -1])
def test_log_distribution_rejects_nonpositive_minimum(lo):
    with pytest.raises(ValueError, match="0 < min < max"):
        Variable("x", "", "Log", {"min": lo, "max": 10})
    report = validate_task_config({"Sampling": {"Method": "Random", "Bounds": {"point_number": 2},
        "Variables": [{"name": "x", "distribution": {"type": "Log", "parameters": {"min": lo, "max": 10}}}]}})
    assert any("minimum > 0" in item.message for item in report.errors())


def test_libdeps_installation_path_expands_base_before_decoding(tmp_path):
    parser = CommandParser.from_config({"LibDeps": {"path": str(tmp_path / "deps"), "Modules": [
        {"name": "Tool", "installation": {"path": "${LibDeps:path}/Tool"}}]}}, project_root=str(tmp_path))
    assert parser.resolve_static("${LibDeps:Tool}/bin/tool") == str(tmp_path / "deps/Tool/bin/tool")


def test_unknown_dependency_is_not_silently_assigned_a_layer():
    with pytest.raises(ValueError, match="unknown module.*typo"):
        resolve_module_layers([{"name": "A", "required_modules": ["typo"]}])


def test_retry_sends_do_not_inflate_factory_sample_counts():
    factory = TaskFactory.__new__(TaskFactory)
    factory.redis = mock.Mock()
    factory.redis.fetch_sample_stats.return_value = {"completed": 40, "failed": 0, "running": 0}
    factory.redis.get_queue_lengths.return_value = {"task_queue_length": 0}
    factory.redis.get_op_count.return_value = 41
    factory.workers = []
    factory._workers_total = 1
    factory._peak_workers_alive = 1
    factory._run_started_at = None
    assert factory.get_run_metrics()["submitted"] == 40


def test_outcome_submitted_matches_real_worker_tasks():
    core = Jarvis2Core()
    core._live_sample_counters = lambda: {"ok": 10, "failed": 2, "queued": 0, "running": 0}
    core._archiver_records_written = lambda: 12
    core.sampler = mock.Mock()
    core.sampler.summary.return_value = {"ncall": 5}
    outcome = core._capture_run_outcome(submitted=5)
    assert outcome.submitted == 12
    assert outcome.extras["ncall"] == 5


def test_refs_json_is_a_single_object(capsys):
    from jarvishep2.cli.client import main
    assert main(["--json", "--refs"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["ok"] and result["command"] == "refs"
    assert result["framework"] and result["samplers"]


def test_all_failed_exports_still_have_logl_column(tmp_path):
    from jarvishep2.database import SimpleHDF5Writer, convert_hdf5_to_csv
    path = tmp_path / "samples.hdf5"
    SimpleHDF5Writer(str(path)).add_record({"x": 1, "status": "Failed"})
    convert_hdf5_to_csv(str(path))
    assert "LogL" in path.with_suffix(".csv").read_text().splitlines()[0].split(",")


def test_ensemble_acceptance_uses_chain_temperature_and_stretch_jacobian():
    import numpy as np
    from jarvishep2.sampling.Source.MCMC.engine_ensemble import EnsembleChain
    engine = EnsembleChain([.5, .5], .1, 10, 0, lambda _: [], rng=np.random.default_rng(1))
    engine.last_loglikelihood = 0.0
    engine.proposed_param = np.array([.6, .6])
    engine._last_log_jacobian = np.log(2.)
    # At beta=1/4: -2/4 + log(2) > 0, so this move must be accepted.
    assert engine.update(-2., beta=.25)


@pytest.mark.parametrize("on_failure", ["typo", 3])
def test_unknown_failure_policy_is_rejected(on_failure):
    config = {"Scan": {"name": "s"}, "EnvReqs": {}, "Operas": {"Modules": []},
              "Sampling": {"Method": "ToyMCMC", "Bounds": {
                  "num_chains": 2, "num_iters": 4, "proposal_scale": .1,
                  "on_failure": on_failure}, "on_failure": on_failure,
                  "Variables": [{"name": "x", "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}}}]}}
    report = validate_task_config(config)
    assert any("on_failure" in item.path for item in report.errors())


def test_plugin_manual_generational_budget_matches_profile():
    from jarvishep2.cli.man import _plugin_method_page
    record = SimpleNamespace(options_schema=None, description="test", profile="generational",
                             version="1", source="user", source_path="sampler.py")
    with mock.patch("jarvishep2.sampler_catalog.plugin_record", return_value=record):
        page = _plugin_method_page("SliceMCMC")
    keys = {item["name"] for item in page["keys"]}
    assert "max_evaluations" in keys and "point_number" not in keys
    assert "max_evaluations:" in page["examples"][0]


def test_decimal_adaptation_scale_is_valid():
    from jarvishep2.sampling.Source.MCMC.config_contract import bounds_get_float
    assert bounds_get_float({"adapt_scale": 2.38}, "adapt_scale", minimum=1e-9) == 2.38


@pytest.mark.parametrize("field,value", [("radius", 0), ("adapt_scale", -2.38)])
def test_numeric_bound_diagnostic_does_not_call_a_number_the_wrong_type(field, value):
    from jarvishep2.cardload.task_schema import _schema_error_message
    error = SimpleNamespace(validator="oneOf", instance=value,
                            schema={"x-jarvis-numeric": True, "oneOf": [
                                {"type": "number", "exclusiveMinimum": 0}, {"type": "string"}]})
    assert _schema_error_message(error) == f"expected a number > 0, got {value!r}"


def test_project_module_is_importable_without_pythonpath_in_spawn(tmp_path, monkeypatch):
    import os
    import subprocess
    import sys
    import textwrap
    import yaml
    project = tmp_path / "project"
    (project / "bin").mkdir(parents=True)
    (project / ".jarvis-project.json").write_text("{}")
    (project / "local_physics.py").write_text("def evaluate(pipe):\n    pipe.send(42)\n    pipe.close()\n")
    card = project / "bin/task.yaml"
    card.write_text(yaml.safe_dump({"Scan": {"name": "s"}, "EnvReqs": {},
        "Sampling": {"Method": "Random", "Bounds": {"point_number": 1},
        "Variables": [{"name": "x", "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}}}]},
        "Operas": {"Modules": []}}))
    script = tmp_path / "runner.py"
    script.write_text(textwrap.dedent(f"""
        import multiprocessing as mp
        from jarvishep2.core import Jarvis2Core
        if __name__ == '__main__':
            Jarvis2Core().load_task_yaml({str(card)!r}, validate=False)
            from local_physics import evaluate
            ctx = mp.get_context('spawn')
            receive, send = ctx.Pipe(duplex=False)
            child = ctx.Process(target=evaluate, args=(send,))
            child.start()
            send.close()
            assert receive.poll(10)
            assert receive.recv() == 42
            child.join(10)
            assert child.exitcode == 0
    """))
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    result = subprocess.run([sys.executable, str(script)], cwd=tmp_path, env=env,
                            capture_output=True, text=True, timeout=25)
    assert result.returncode == 0, result.stderr
