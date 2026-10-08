"""External sampler S1: Jarvis-HEP-Sampler plug-ins behave like built-in samplers.

Covers the parity contract of ``docs/jarvis-hep-sampler/technical-design.md``
§1.2 (P1–P6, P9) for the shipped ``LatinHypercube`` and for user plug-ins
loaded with ``Jarvis sampler load``: catalog and Distributor, validation and
``Jarvis man`` without running plug-in code, the adapter against the SDK's
reference stream, batch invariance, checkpoints, preflight, logging,
telemetry, monitor, the CLI, and a full ``core.run``.
"""

from __future__ import annotations

import json
import logging
import os
import textwrap
import threading
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest import mock

import numpy as np
import pytest
import yaml
from jarvis_hep_sampler import SamplerRegistry
from jarvis_hep_sampler.errors import PluginConfigError, PluginOutputError
from jarvis_hep_sampler.testing import reference_stream

from jarvishep2 import sampler_catalog
from jarvishep2.distributor import STATELESS_METHODS, Distributor
from jarvishep2.logging.toplevel import JarvisContextFormatter
from jarvishep2.sampling.lifecycle_log import SamplerLifecycleLog
from jarvishep2.sampling.plugin_sampler import PluginSampler, preflight_plugin

TESTS_ROOT = Path(__file__).resolve().parent
RANDOM_YAML = TESTS_ROOT / "parity_project" / "random_opera.yaml"

HALF_CUBE = """
import os
from pathlib import Path

import numpy as np

_MARKER = os.environ.get("JHS_TEST_MARKER")
if _MARKER:  # proves which processes import the plug-in
    with open(_MARKER, "a", encoding="utf-8") as handle:
        handle.write(f"{os.getpid()}\\n")

PLUGIN = {
    "api_version": 1,
    "method": "HalfCube",
    "version": "1.0",
    "profile": "independent",
    "description": "Uniform points below an upper edge on every axis.",
}

OPTIONS_SCHEMA = {
    "type": "object",
    "properties": {
        "upper": {"type": "number", "exclusiveMinimum": 0, "maximum": 1, "default": 0.5},
        "explode_at": {"type": "integer"},
    },
    "additionalProperties": False,
}


def prepare(setup):
    setup.log.step("HalfCube Plan ->", [("upper", setup.options.get("upper", 0.5))])
    return {"upper": float(setup.options.get("upper", 0.5))}


def draw(context, rng, plan):
    if context.candidate_index == context.options.get("explode_at"):
        return [float("nan")] * context.dimension
    return rng.random(context.dimension) * plan["upper"]


def settings_rows(setup, plan):
    return [("upper", plan["upper"])]


def summary_rows(plan, stats):
    return [("kept fraction", round(stats.accepted / max(1, stats.candidates), 3))]


def status(plan, stats):
    return {"upper": plan["upper"]}
"""

NO_SCHEMA = """
PLUGIN = {"api_version": 1, "method": "FreeForm", "version": "1", "profile": "independent"}


def draw(context, rng, plan):
    return rng.random(context.dimension)
"""


# ------------------------------------------------------------------ fixtures
@pytest.fixture
def home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    path = tmp_path / "samplers-home"
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(path))
    return path


def _write(tmp_path: Path, name: str, source: str) -> Path:
    path = tmp_path / name
    path.write_text(textwrap.dedent(source), encoding="utf-8")
    return path


@pytest.fixture
def half_cube(tmp_path: Path, home: Path) -> Path:
    path = _write(tmp_path, "HalfCube.py", HALF_CUBE)
    SamplerRegistry().load(path, reserved_names=sampler_catalog.reserved_method_names())
    return path


def _config(tmp_path: Path, method: str, bounds: dict[str, Any], **sampling: Any) -> dict[str, Any]:
    variables = [
        {"name": name, "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}}}
        for name in ("x", "y")
    ]
    return {
        "Scan": {"name": "plugin-scan"},
        "scan_name": "plugin-scan",
        "task_root": str(tmp_path),
        "Runtime": {"mode": "redis", "workers": 2, "batch_size": 4},
        "Sampling": {"Method": method, "Variables": variables, "Bounds": dict(bounds), **sampling},
    }


def _sampler(config: dict[str, Any]) -> PluginSampler:
    sampler = Distributor.set_method(config["Sampling"]["Method"])
    assert isinstance(sampler, PluginSampler)
    sampler.set_config(config)
    return sampler


def _drain(sampler: PluginSampler) -> list[Any]:
    samples = []
    while (sample := sampler.propose_next()) is not None:
        samples.append(sample)
    return samples


class _Capture(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []
        self.setFormatter(JarvisContextFormatter(colorize=False))

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)

    def __enter__(self) -> _Capture:
        base = logging.getLogger("jarvis_hep")
        self._old = base.level
        base.addHandler(self)
        base.setLevel(logging.DEBUG)
        return self

    def __exit__(self, *exc: Any) -> None:
        base = logging.getLogger("jarvis_hep")
        base.removeHandler(self)
        base.setLevel(self._old)

    def messages(self) -> list[str]:
        return [record.getMessage() for record in self.records]


# ---------------------------------------------------------- catalog / Distributor
def test_shipped_plugin_flags_match_its_file(home: Path) -> None:
    for spec in sampler_catalog.builtin_plugin_specs():
        record = sampler_catalog.plugin_record(spec.name)
        assert record is not None and record.method == spec.name
        assert spec.stateless == (record.profile == "independent")
        assert spec.statistical == record.statistical
        assert spec.family == "plugin" and spec.source == "builtin-plugin"


def test_user_plugin_joins_the_catalog_and_distributor(half_cube: Path) -> None:
    spec = sampler_catalog.spec_for("HalfCube")
    assert spec is not None and spec.source == "user" and spec.family == "plugin"
    assert "HalfCube" in sampler_catalog.available_methods()
    assert "HalfCube" in Distributor.available_methods()
    assert "HalfCube" in STATELESS_METHODS
    assert Distributor.get_resume_status("HalfCube") == "implemented"
    assert sampler_catalog.canonical_name("halfcube") == "HalfCube"
    assert sampler_catalog.is_stateless("HalfCube")
    assert isinstance(Distributor.set_method("HalfCube"), PluginSampler)


def test_a_user_entry_never_shadows_a_shipped_method(tmp_path: Path, home: Path) -> None:
    clone = _write(tmp_path, "lhs.py", HALF_CUBE.replace('"HalfCube"', '"LatinHypercube"'))
    SamplerRegistry().load(clone)  # e.g. loaded before Jarvis-HEP shipped the name
    assert sampler_catalog.spec_for("LatinHypercube").source == "builtin-plugin"


def test_unknown_method_points_at_jarvis_sampler_load(home: Path) -> None:
    with pytest.raises(NotImplementedError, match="Jarvis sampler load"):
        Distributor.set_method("NotLoadedYet")


def test_unreadable_registry_keeps_builtins_and_is_reported(home: Path) -> None:
    home.mkdir(parents=True)
    (home / "sources.json").write_text("{broken", encoding="utf-8")
    assert "Random" in sampler_catalog.available_methods()
    assert "not valid JSON" in (sampler_catalog.user_registry_error() or "")


# ------------------------------------------------------------------ adapter
def test_points_match_the_sdk_reference_stream(tmp_path: Path, half_cube: Path) -> None:
    bounds = {"seed": 5, "point_number": 40, "upper": 0.8}
    samples = _drain(_sampler(_config(tmp_path, "HalfCube", bounds, selection="x < 0.4")))
    expected = list(
        reference_stream(
            str(sampler_catalog.plugin_record("HalfCube").snapshot_path),
            point_number=40,
            seed=5,
            variable_names=["x", "y"],
            options={"upper": 0.8},
            selection=lambda u: bool(u[0] < 0.4),
        )
    )
    assert [s.sample_index for s in samples] == [p.accepted_index for p in expected]
    assert all(np.array_equal(s.u_coords, p.u) for s, p in zip(samples, expected, strict=True))
    assert len({s.uuid for s in samples}) == len(samples)
    again = _drain(_sampler(_config(tmp_path, "HalfCube", bounds, selection="x < 0.4")))
    assert [s.uuid for s in again] == [s.uuid for s in samples]


def test_identity_drives_uuids(tmp_path: Path, half_cube: Path) -> None:
    base = {"seed": 5, "point_number": 3, "upper": 0.8}

    def first_uuid(bounds: dict[str, Any]) -> str:
        return _drain(_sampler(_config(tmp_path, "HalfCube", bounds)))[0].uuid

    reference = first_uuid(base)
    assert first_uuid(dict(base)) == reference
    assert first_uuid({**base, "upper": 0.7}) != reference
    assert first_uuid({**base, "seed": 6}) != reference


def test_submission_is_batch_invariant(tmp_path: Path) -> None:
    from jarvishep2.redis_queue import make_fakeredis_queue

    seen = []
    for batch_size in (1, 7, 64):
        # No Workers drain the fake queue, so the in-flight limit must cover the
        # whole budget (Runtime normalisation drops max_inflight; Bounds.max_worker
        # is the reserved key FixedSetSampler reads).
        config = _config(
            tmp_path, "LatinHypercube", {"seed": 3, "point_number": 50, "max_worker": 1000}
        )
        config["Runtime"]["batch_size"] = batch_size
        sampler = _sampler(config)
        queue = make_fakeredis_queue()
        queue.connect()
        sampler.set_redis(queue)
        assert sampler.run_distributed() == 50
        tasks = {}
        while (task := queue.pull_task(timeout=1)) is not None:
            tasks[task["sample_index"]] = (task["uuid"], tuple(task["u_coords"]))
        seen.append(tasks)
    assert len(seen[0]) == 50 and seen[0] == seen[1] == seen[2]


def test_bad_output_stops_at_the_candidate(tmp_path: Path, half_cube: Path) -> None:
    sampler = _sampler(_config(tmp_path, "HalfCube", {"point_number": 10, "explode_at": 3}))
    with pytest.raises(PluginOutputError) as info:
        _drain(sampler)
    assert info.value.candidate_index == 3


def test_non_finite_mapped_value_is_rejected(tmp_path: Path) -> None:
    config = _config(
        tmp_path,
        "LatinHypercube",
        {"seed": 1, "point_number": 2, "centered": True},
        Mapper=[{"name": "lx", "expression": "log(x - 0.25)"}],
    )
    with pytest.raises(PluginOutputError, match="lx = -inf"):
        _drain(_sampler(config))


def test_checkpoint_round_trip_and_prefix_replay(tmp_path: Path, half_cube: Path) -> None:
    config = _config(tmp_path, "HalfCube", {"seed": 2, "point_number": 30}, selection="y < 0.3")
    uninterrupted = _drain(_sampler(config))

    sampler = _sampler(config)
    sampler.assert_checkpoint_attribute_contract()
    head = [sampler.propose_next() for _ in range(4)]
    state = sampler.export_runtime_state()
    assert "numpy_random_state" not in state and state["plugin"]["identity"]
    resumed = _sampler(config)
    resumed.import_runtime_state(state)
    tail = _drain(resumed)
    assert [s.uuid for s in head + tail] == [s.uuid for s in uninterrupted]

    replayed = _sampler(config)
    replayed.advance_to_persisted_prefix(6)
    assert [s.uuid for s in _drain(replayed)] == [s.uuid for s in uninterrupted[6:]]


# ---------------------------------------------------------------- preflight
def test_preflight_is_a_no_op_for_builtins(tmp_path: Path) -> None:
    assert preflight_plugin(_config(tmp_path, "Random", {"point_number": 3})) is None


def test_preflight_keeps_a_scan_copy_that_resume_prefers(tmp_path: Path, half_cube: Path) -> None:
    config = _config(tmp_path, "HalfCube", {"point_number": 5})
    first = preflight_plugin(config, task_root=str(tmp_path), scan_name="s1")
    assert first.source_copy is not None and first.source_copy.is_file()
    assert first.source_copy.parent == tmp_path / "checkpoints" / "s1" / "HalfCube" / "plugin"

    half_cube.write_text(
        textwrap.dedent(HALF_CUBE).replace('"version": "1.0",', '"version": "1.0",  # edited'),
        encoding="utf-8",
    )
    SamplerRegistry().load(half_cube)  # the registry now holds a different file
    with _Capture() as capture:
        resumed = preflight_plugin(config, task_root=str(tmp_path), scan_name="s1", resume=True)
    assert resumed.resumed_from_copy
    assert resumed.plugin.source_sha256 == first.plugin.source_sha256
    assert any("own copy" in message for message in capture.messages())
    fresh = preflight_plugin(config, task_root=str(tmp_path), scan_name="s2")
    assert fresh.plugin.source_sha256 != first.plugin.source_sha256


def test_preflight_failure_is_logged_with_the_plugin_traceback(tmp_path: Path, half_cube: Path) -> None:
    from jarvishep2.runtime.core import Jarvis2Core

    config = _config(tmp_path, "HalfCube", {"point_number": 5, "upper": 2.0})
    core = SimpleNamespace(
        config=config,
        info={"task_root": str(tmp_path), "scan_name": "s"},
        _resume_policy="fresh",
        _logger=logging.getLogger("jarvis_hep.test"),
        _plugin_preflight=None,
    )
    with _Capture() as capture, pytest.raises(PluginConfigError) as info:
        Jarvis2Core.preflight_sampler_plugin(core)
    assert getattr(info.value, "_jarvis_logged", False)
    failure = next(r for r in capture.records if "meets error when preparing" in r.getMessage())
    assert failure.exc_info is not None
    assert "Jarvis-HEP.Sampler.HalfCube" in capture.format(failure)


def test_bootstrap_runs_the_preflight_before_redis(tmp_path: Path) -> None:
    from jarvishep2.runtime import _runtime_supervisor

    calls: list[str] = []

    class _Stop(Exception):
        pass

    core = mock.MagicMock()
    core._resume_policy = "fresh"
    core._resolve_database_dir.return_value = str(tmp_path)
    core.preflight_sampler_plugin.side_effect = lambda: calls.append("preflight")
    core.check_environment_requirements.side_effect = lambda: calls.append("environment")
    core.prepare_libraries.side_effect = lambda **_: calls.append("libraries")

    def _redis() -> None:
        calls.append("redis")
        raise _Stop

    core.init_redis.side_effect = _redis
    with (
        mock.patch("jarvishep2.process_cleanup.ensure_scan_name_available"),
        mock.patch.object(_runtime_supervisor, "build_command_parser"),
        pytest.raises(_Stop),
    ):
        _runtime_supervisor._RuntimeSupervisor(core).bootstrap_distributed_runtime()
    assert calls == ["environment", "preflight", "libraries", "redis"]


# ---------------------------------------------------- logging / telemetry / monitor
def test_lifecycle_log_and_plugin_log_facade(tmp_path: Path, half_cube: Path) -> None:
    with _Capture() as capture:
        sampler = _sampler(_config(tmp_path, "HalfCube", {"seed": 9, "point_number": 8}))
        _drain(sampler)
        log = SamplerLifecycleLog(sampler)
        log.start(resume=False, workers=2)
        log.result(SimpleNamespace(status="success", submitted=8, completed=8, failed=0))
    rendered = [capture.format(record) for record in capture.records]
    assert all("Jarvis-HEP.Sampler.HalfCube" in text for text in rendered)
    messages = capture.messages()
    plan = next(m for m in messages if m.startswith("HalfCube Plan ->"))
    assert "Field" in plan and "upper" in plan
    settings = next(m for m in messages if m.startswith("HalfCube Sampler Settings ->"))
    for key in ("seed", "point number", "plugin", "HalfCube 1.0", "upper"):
        assert key in settings
    summary = next(m for m in messages if m.startswith("HalfCube Sampler Summary ->"))
    assert "point_number reached" in summary and "kept fraction" in summary
    assert not [m for m in messages if "could not report" in m]


def test_telemetry_metadata_and_monitor(tmp_path: Path, half_cube: Path) -> None:
    from jarvishep2.monitor.queues import feedback_queue_plan
    from jarvishep2.monitor.samples import (
        SamplesBlockData,
        render_samples_block,
        sampler_display_from_sources,
    )
    from jarvishep2.runtime_metadata import sampler_metadata_from_config
    from jarvishep2.sampling.telemetry import sampler_status_json

    config = _config(tmp_path, "HalfCube", {"seed": 4, "point_number": 20, "upper": 0.6})
    sampler = _sampler(config)
    for _ in range(5):
        sampler.propose_next()
    status = json.loads(sampler_status_json(sampler))
    assert status["progress"] == {"kind": "finite", "current": 5, "target": 20}
    assert status["metrics"]["plugin_status"] == {"upper": 0.6}

    metadata = sampler_metadata_from_config(config)
    assert metadata["family"] == "plugin"
    assert metadata["config"]["version"] == "1.0"
    assert "upper" not in metadata["config"]  # options stay out of run metadata

    display = sampler_display_from_sources(metadata, {"sampler_status": json.dumps(status)})
    block = render_samples_block(60, SamplesBlockData(sampler=display))
    assert "CANDIDATES 5 / 20" in block.body
    assert block.border_subtitle.startswith("PLUGIN 1.0 · ACCEPTED 5 · SEED 4 · UPPER")
    assert feedback_queue_plan("HalfCube") == ("inactive", ())


def test_slow_status_cannot_block_control_lease(tmp_path: Path, half_cube: Path, monkeypatch) -> None:
    from jarvishep2.runtime import _runtime_supervisor as runtime
    from jarvishep2.sampling.telemetry import sampler_status_json

    sampler = _sampler(_config(tmp_path, "HalfCube", {"point_number": 10, "upper": .6}))
    entered, release, renewed, stop = (threading.Event() for _ in range(4))
    calls = []

    def slow_status(stats):
        calls.append(threading.current_thread().name)
        entered.set()
        assert release.wait(5), "test did not release the reporting hook"
        return {"upper": .6}

    monkeypatch.setattr(sampler._session, "status", slow_status)
    redis = mock.Mock(config={"port": 6379})
    redis.ping.return_value = True
    redis.refresh_control_lock.side_effect = lambda *a, **kw: (renewed.set() or True)
    core = SimpleNamespace(redis=redis, sampler=sampler, _shutdown_done=False,
                           _interrupt_requested=False, _managed_redis=None, archiver=None,
                           _control_lock_owner="review", _logger=mock.Mock())
    supervisor = runtime._RuntimeSupervisor(core)
    monkeypatch.setattr(supervisor, "_ensure_archiver_alive", lambda: None)
    monkeypatch.setattr(runtime, "CONTROL_LOCK_TTL_SEC", 3)
    monkeypatch.setattr(runtime, "_SUPERVISE_TICK_SEC", .01)
    sampling = threading.Thread(target=lambda: sampler._refresh_plugin_status(force=True),
                                name="TestSampling")
    lease = threading.Thread(target=supervisor._control_lease_loop, args=(stop, "review"),
                             name="TestLease")
    try:
        sampling.start()
        assert entered.wait(2)
        lease.start()
        assert renewed.wait(2), "a slow hook prevented the control lock refresh"
        status = json.loads(sampler_status_json(sampler))
        assert status["metrics"]["plugin_status"] == {"upper": .6}
        assert calls == ["TestSampling"]
        assert redis.publish_proc_board.called
    finally:
        release.set()
        stop.set()
        sampling.join(2)
        if lease.ident is not None:
            lease.join(2)
    assert not sampling.is_alive() and not lease.is_alive()


@pytest.mark.parametrize("mode", ["none", "valid", "broken"])
def test_explicit_options_are_logged_without_a_reporting_hook(tmp_path: Path, home: Path, mode) -> None:
    source = HALF_CUBE.split("def settings_rows")[0]
    if mode == "none":
        start, end = source.index("OPTIONS_SCHEMA ="), source.index("def prepare")
        source = source[:start] + source[end:]
    if mode == "broken":
        source += '\ndef settings_rows(setup, plan):\n    raise ValueError("reporting broke")\n'
    SamplerRegistry().load(_write(tmp_path, "HalfCube.py", source))
    sampler = _sampler(_config(tmp_path, "HalfCube", {"point_number": 4, "upper": .7}))
    assert dict(sampler.log_settings_rows())["upper"] == .7
    assert sampler.propose_next() is not None


def test_metadata_uses_the_prepared_source_after_reload(tmp_path: Path, half_cube: Path) -> None:
    from jarvishep2.runtime import _runtime_supervisor

    config = _config(tmp_path, "HalfCube", {"point_number": 8, "seed": 5})
    preflight = preflight_plugin(config, task_root=str(tmp_path), scan_name="plugin-scan")
    old_hash = preflight.plugin.source_sha256
    half_cube.write_text(HALF_CUBE.replace('"version": "1.0"', '"version": "2.0"'))
    SamplerRegistry().load(half_cube)
    resumed = preflight_plugin(config, task_root=str(tmp_path), scan_name="plugin-scan", resume=True)
    redis = mock.Mock()
    redis.connection_config.return_value = {"host": "localhost", "port": 1234, "db": 0}
    core = SimpleNamespace(config=config, info={"task_result_dir": str(tmp_path)}, redis=redis,
                           _plugin_preflight=resumed)
    _runtime_supervisor._RuntimeSupervisor(core)._publish_runtime_metadata()
    data = json.loads(Path(core.info["runtime_metadata_path"]).read_text())
    knobs = data["sampler"]["config"]
    assert knobs["version"] == "1.0"
    assert knobs["source_sha256"] == old_hash
    assert sampler_catalog.plugin_record("HalfCube").version == "2.0"


def test_stop_is_observed_between_rejected_candidates_and_in_replay(tmp_path: Path, half_cube: Path) -> None:
    for replay in (False, True):
        sampler = _sampler(_config(tmp_path, "HalfCube", {"point_number": 100000}, selection="x < 0"))
        sampler.set_interrupt_check(lambda: sampler._index >= 3)
        with pytest.raises(KeyboardInterrupt):
            sampler.advance_to_persisted_prefix(1) if replay else sampler.propose_next()
        assert sampler._index == 3 and sampler._accepted_index == 0


def test_interrupt_and_scan_mode_stop_backpressure_wait(tmp_path: Path, half_cube: Path) -> None:
    sampler = _sampler(_config(tmp_path, "HalfCube", {"point_number": 10}))
    sampler.redis = mock.Mock()
    sampler.redis.fetch_sample_stats.return_value = {"completed": 0, "failed": 0}
    sampler.redis.read_proc_board.return_value = {"scan_mode": "paused", "pause_reason": "workers down"}
    with pytest.raises(RuntimeError, match="scan_mode=paused"):
        sampler._wait_for_inflight_room(2, base_finished=0, max_inflight=2)
    sampler.set_interrupt_check(lambda: True)
    with pytest.raises(KeyboardInterrupt):
        sampler._wait_for_inflight_room(2, base_finished=0, max_inflight=2)


def test_indexed_block_adapter_against_a_known_sequence(tmp_path: Path, home: Path) -> None:
    source = '''
PLUGIN = {"api_version":1,"method":"KnownBlock","version":"1","profile":"independent","block_size":4}
def draw_block(context, rng, plan):
    return [[(context.block_index * context.block_size + row) / context.point_number]
            * context.dimension for row in range(context.block_size)]
'''
    SamplerRegistry().load(_write(tmp_path, "KnownBlock.py", source))
    config = _config(tmp_path, "KnownBlock", {"point_number": 11}, selection="x > 0.3")
    points = _drain(_sampler(config))
    assert [p.sample_index for p in points] == list(range(7))
    assert np.array_equal(np.array([p.u_coords[0] for p in points]), np.arange(4, 11) / 11)


# --------------------------------------------------------------- validation
def _card(tmp_path: Path, method: str, bounds: dict[str, Any]) -> Path:
    card = yaml.safe_load(RANDOM_YAML.read_text(encoding="utf-8"))
    card["Sampling"]["Method"] = method
    card["Sampling"]["Bounds"] = bounds
    path = tmp_path / f"{method}-{len(list(tmp_path.glob('*.yaml')))}.yaml"
    path.write_text(yaml.safe_dump(card), encoding="utf-8")
    return path


def _plugin_issues(report: Any) -> list[tuple[str, str]]:
    return [(i.level, i.code) for i in report.issues if i.code.startswith("JV2-PLG")]


def _validate(path: Path, *, strict: bool = False) -> Any:
    from jarvishep2.task_config import load_task_yaml
    from jarvishep2.task_validation import validate_task_config

    return validate_task_config(load_task_yaml(str(path)), strict=strict)


def test_validation_of_plugin_cards(tmp_path: Path, half_cube: Path) -> None:
    ok = _validate(_card(tmp_path, "HalfCube", {"point_number": 10, "upper": 0.3}))
    assert ok.ok and not _plugin_issues(ok)

    typo = _validate(_card(tmp_path, "HalfCube", {"point_number": 10, "uper": 0.3}))
    [problem] = typo.errors()
    assert (problem.code, problem.path) == ("JV2-PLG-004", "Sampling.Bounds.uper")
    assert "did you mean 'upper'" in problem.message

    wrong = _validate(_card(tmp_path, "HalfCube", {"point_number": 10, "upper": 2}))
    assert [(i.code, i.path) for i in wrong.errors()] == [("JV2-PLG-004", "Sampling.Bounds.upper")]

    reserved = _validate(_card(tmp_path, "HalfCube", {"point_number": 0}))
    assert [(i.code, i.path) for i in reserved.errors()] == [
        ("JV2-PLG-003", "Sampling.Bounds.point_number")
    ]

    unknown = _validate(_card(tmp_path, "NotLoaded", {"point_number": 1}))
    assert unknown.errors()[0].code == "JV2-MTH-003"
    assert "Jarvis sampler load" in unknown.errors()[0].suggestion


def test_options_without_a_schema_are_a_note_not_a_failure(tmp_path: Path, home: Path) -> None:
    from jarvishep2.task_validation import format_report

    SamplerRegistry().load(_write(tmp_path, "FreeForm.py", NO_SCHEMA))
    for strict in (False, True):
        report = _validate(_card(tmp_path, "FreeForm", {"point_number": 4, "knob": 1}), strict=strict)
        assert report.ok
        assert _plugin_issues(report) == [("info", "JV2-PLG-001")]
        assert report.notes() and not report.warnings()
    assert "successful (1 note(s))" in format_report(report)


def test_broken_snapshot_is_reported(tmp_path: Path, half_cube: Path) -> None:
    snapshot = Path(sampler_catalog.plugin_record("HalfCube").snapshot_path)
    snapshot.write_text(snapshot.read_text(encoding="utf-8") + "\n# tampered\n", encoding="utf-8")
    report = _validate(_card(tmp_path, "HalfCube", {"point_number": 10}))
    assert [i.code for i in report.errors()] == ["JV2-PLG-002"]


def test_validate_man_and_list_never_run_plugin_code(
    tmp_path: Path, home: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from jarvishep2.cli.man import resolve_man_request
    from jarvishep2.cli.sampler import dispatch_sampler

    marker = tmp_path / "imports.txt"
    monkeypatch.setenv("JHS_TEST_MARKER", str(marker))
    SamplerRegistry().load(_write(tmp_path, "HalfCube.py", HALF_CUBE))
    marker.unlink()
    _validate(_card(tmp_path, "HalfCube", {"point_number": 10}))
    resolve_man_request(["sampler.HalfCube"])
    resolve_man_request(["sampler"])
    assert dispatch_sampler(["list"]) == 0
    assert dispatch_sampler(["info", "HalfCube"]) == 0
    assert not marker.exists()


def test_man_pages_for_plugins(half_cube: Path) -> None:
    from jarvishep2.cli.man import resolve_man_request

    page = resolve_man_request(["sampler.halfcube"])
    assert page["context"] == {"method": "HalfCube"}
    assert [key["name"] for key in page["keys"]] == ["seed", "point_number", "upper", "explode_at"]
    assert page["capabilities"]["type"] == "plugin"
    assert "upper: 0.5" in page["examples"][0]
    index = {row["method"]: row["type"] for row in resolve_man_request(["sampler"])["methods"]}
    assert index["HalfCube"] == "plugin" and index["LatinHypercube"] == "plugin"
    assert resolve_man_request([], code="JV2-PLG-004")["diagnostics"][0] == "JV2-PLG-004"


def test_local_schema_refs_retain_type_description_and_default(tmp_path: Path, home: Path) -> None:
    from jarvishep2.cli.man import resolve_man_request

    schema = {"type": "object", "$defs": {"upper": {"type": "number", "maximum": 1,
              "default": .5, "description": "Upper edge of the unit cube."}},
              "properties": {"upper": {"$ref": "#/$defs/upper"}}}
    SamplerRegistry().load(_write(tmp_path, "FreeForm.py", NO_SCHEMA + '\nOPTIONS_SCHEMA = ' + repr(schema)))
    page = resolve_man_request(["sampler.FreeForm"])
    entry = next(k for k in page["keys"] if k["name"] == "upper")
    assert (entry["type"], entry["default"], entry["description"]) == ("number", .5, "Upper edge of the unit cube.")
    assert "upper: 0.5" in page["examples"][0]
    assert _validate(_card(tmp_path, "FreeForm", {"point_number": 10, "upper": .5})).ok


# --------------------------------------------------------------------- CLI
def test_cli_lifecycle(tmp_path: Path, home: Path, capsys: pytest.CaptureFixture[str]) -> None:
    from jarvishep2.cli.sampler import dispatch_sampler

    target = tmp_path / "MyScan.py"
    assert dispatch_sampler(["init", "MyScan", "--output", str(target)]) == 0
    assert dispatch_sampler(["init", "MyScan", "--output", str(target)]) != 0
    assert dispatch_sampler(["init", "Grid"]) != 0
    assert dispatch_sampler(["check", str(target), "--dimension", "3", "--points", "20"]) == 0
    assert dispatch_sampler(["load", str(target)]) == 0
    assert "Sampling.Method: MyScan" in capsys.readouterr().out
    assert dispatch_sampler(["check", "LatinHypercube", "--points", "20"]) == 0
    assert dispatch_sampler(["list"]) == 0
    listing = capsys.readouterr().out
    assert "MyScan" in listing and "LatinHypercube" in listing
    assert dispatch_sampler(["remove", "LatinHypercube", "--yes"]) != 0
    assert dispatch_sampler(["remove", "MyScan", "--yes"]) == 0
    assert sampler_catalog.spec_for("MyScan") is None
    assert target.exists()


def test_cli_check_reports_a_broken_plugin(tmp_path: Path, home: Path) -> None:
    from jarvishep2.cli.sampler import dispatch_sampler

    bad = _write(tmp_path, "Bad.py", NO_SCHEMA.replace("rng.random(context.dimension)", "[2.0] * context.dimension"))
    assert dispatch_sampler(["check", str(bad)]) != 0
    assert dispatch_sampler(["load", str(tmp_path / "missing.py")]) != 0


def test_cli_check_preserves_card_budget_and_reports_overrides(tmp_path: Path, half_cube: Path, capsys) -> None:
    from jarvishep2.cli.sampler import dispatch_sampler

    card = _card(tmp_path, "HalfCube", {"point_number": 128, "seed": 7, "upper": .5})
    with mock.patch("jarvis_hep_sampler.testing.check_plugin") as checker:
        checker.return_value = SimpleNamespace(method="HalfCube", ok=True, rows=lambda: [])
        assert dispatch_sampler(["check", str(half_cube), "--card", str(card)]) == 0
        assert checker.call_args.kwargs["point_number"] == 128
        assert checker.call_args.kwargs["seed"] == 7
        assert checker.call_args.kwargs["options"] == {"upper": .5}
        assert checker.call_args.kwargs["variable_names"] == ["x", "y"]
        assert "card overrides" in capsys.readouterr().out
        assert dispatch_sampler(["check", str(half_cube), "--card", str(card), "--points", "16", "--seed", "3"]) == 0
        assert checker.call_args.kwargs["point_number"] == 16
        assert checker.call_args.kwargs["seed"] == 3
        assert "changes the total candidate budget" in capsys.readouterr().out


# --------------------------------------------------------------- end to end
def test_core_run_with_a_user_plugin(tmp_path: Path, half_cube: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """P3: points go through Redis, Workers and the Archiver; Workers never import the plug-in."""
    from test_samplers_catalog import _stop_factory_workers
    from test_worker_mvp import _start_tcp_fakeredis

    from jarvishep2.database import SimpleHDF5Writer
    from jarvishep2.runtime.core import Jarvis2Core

    marker = tmp_path / "imports.txt"
    monkeypatch.setenv("JHS_TEST_MARKER", str(marker))
    card = _card(tmp_path, "HalfCube", {"point_number": 8, "seed": 13, "upper": 0.9})
    _stop_factory_workers()
    server, redis_config = _start_tcp_fakeredis()
    try:
        core = Jarvis2Core()
        core.load_task_yaml(str(card))
        result_dir = tmp_path / "result"
        core.config["task_result_dir"] = str(result_dir)
        core.config["Runtime"]["redis"] = redis_config
        core.config["Runtime"]["Watchdog"] = {"enabled": False}
        core.runtime = core.config["Runtime"]
        core._populate_info_from_config()
        outcome = core.run(write_run_summary=False)
    finally:
        _stop_factory_workers()
        server.shutdown()
        server.server_close()
    assert (outcome.completed, outcome.failed) == (8, 0)
    records = SimpleHDF5Writer(str(result_dir / "DATABASE" / "samples.hdf5")).read_records()
    expected = list(
        reference_stream(
            sampler_catalog.plugin_record("HalfCube").snapshot_path,
            point_number=8,
            seed=13,
            variable_names=["x", "y"],
            options={"upper": 0.9},
        )
    )
    got = sorted((int(r["sample_index"]), float(r["x"]), float(r["y"])) for r in records)
    assert [(i, round(x, 12), round(y, 12)) for i, x, y in got] == [
        (p.accepted_index, round(float(p.u[0]), 12), round(float(p.u[1]), 12)) for p in expected
    ]
    importers = set(marker.read_text(encoding="utf-8").split())
    assert importers == {str(os.getpid())}
    copies = list((tmp_path / "checkpoints").rglob("HalfCube-*.py"))
    assert len(copies) == 1
