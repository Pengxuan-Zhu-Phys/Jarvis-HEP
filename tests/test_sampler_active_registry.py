"""Real scans keep their bound source while separate CLI processes update the registry."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import signal
import socket
import subprocess
import sys
import time

import pytest
import yaml

from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.sampling.stateless_batch import deterministic_sampler_uuid

from test_sampler_cli_lifecycle import (
    CLI_BOOTSTRAP,
    SOURCE,
    SamplerCLI,
    _card,
    _free_port,
    _project,
)

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis is required"),
]

BUDGET = 37
GATE_INDEX = 16


def _active_source(path, version, coordinate):
    source = SOURCE.replace("VERSION", version).replace("COORDINATE", str(coordinate))
    source = source.replace("import os\n", "import os\nimport time\n")
    source = source.replace("    return [VALUE] * context.dimension", '''
    release = Path(os.environ["JHS_ACTIVE_RELEASE"])
    if context.candidate_index >= 16 and not release.exists():
        Path(os.environ["JHS_ACTIVE_ENTERED"]).write_text(str(context.candidate_index))
        deadline = time.monotonic() + 90
        while not release.exists():
            if time.monotonic() >= deadline:
                raise RuntimeError("test controller did not release the draw gate")
            time.sleep(.02)
    value = (context.candidate_index + .5) / context.point_number
    if VALUE == .75:
        value = 1 - value
    return [value] * context.dimension''')
    path.write_text(source, encoding="utf-8")
    return hashlib.sha256(source.encode()).hexdigest()


def _rows(path):
    if not path.exists():
        return []
    return sorted(SimpleHDF5Writer(str(path)).read_records(), key=lambda row: int(row["sample_index"]))


def _assert_rows(rows, binding, count):
    assert len(rows) == count
    prefix = f"plugin:MySampler:{binding['hard_identity'][:16]}"
    result = []
    for index, row in enumerate(rows):
        coordinate = (index + .5) / BUDGET
        assert int(row["sample_index"]) == index
        assert float(row["x"]) == float(row["y"]) == coordinate
        assert row["uuid"] == deterministic_sampler_uuid(prefix=prefix, seed=7, sample_index=index)
        result.append({"index": index, "x": float(row["x"]), "y": float(row["y"]),
                       "uuid": row["uuid"]})
    assert len({row["uuid"] for row in result}) == count
    return result


def _wait_for_durable_gate(process, database, entered, log):
    deadline = time.monotonic() + 75
    while time.monotonic() < deadline:
        assert process.poll() is None, log.read_text(encoding="utf-8")[-8000:]
        if entered.exists():
            try:
                rows = _rows(database)
            except OSError:
                rows = []
            assert len(rows) <= GATE_INDEX
            if len(rows) == GATE_INDEX:
                assert int(entered.read_text()) == GATE_INDEX
                assert [int(row["sample_index"]) for row in rows] == list(range(GATE_INDEX))
                return rows
        time.sleep(.05)
    pytest.fail(f"no live draw gate and durable prefix before deadline: {log}")


def _assert_port_closed(port):
    with socket.socket() as sock:
        sock.settimeout(1)
        assert sock.connect_ex(("127.0.0.1", port)) != 0, "managed Redis leaked after CLI exit"


@pytest.mark.parametrize(
    ("workers", "batch_size", "interrupt"), [(1, 1, False), (2, 4, True)],
    ids=["one-worker-continue", "two-workers-interrupt-resume"],
)
def test_active_scan_survives_registry_reload_and_remove(
    tmp_path, request, workers, batch_size, interrupt,
):
    cli = SamplerCLI(tmp_path)
    project = _project(tmp_path / "active-project")
    source = tmp_path / "MySampler.py"
    release = tmp_path / "release-draw"
    entered = tmp_path / "draw-entered"
    cli.env.update(JHS_ACTIVE_RELEASE=str(release), JHS_ACTIVE_ENTERED=str(entered))
    digest_v1 = _active_source(source, "1.0", .25)
    port = _free_port()
    card = _card(project, "active", port)
    config = yaml.safe_load(card.read_text(encoding="utf-8"))
    config["EnvReqs"]["V2"].update(workers=workers, batch_size=batch_size)
    config["Sampling"]["Bounds"]["point_number"] = BUDGET
    card.write_text(yaml.safe_dump(config), encoding="utf-8")
    cli("load-v1", "sampler", "load", source)
    database = project / "outputs/active/DATABASE/samples.hdf5"
    identity = database.parent / "sampler_identity.json"
    checkpoint = project / "checkpoints/active/MySampler/state.pkl"
    log = cli.logs / "active-run.log"
    command = [sys.executable, "-c", CLI_BOOTSTRAP, "run", str(card), "--skip-draw-flowchart"]
    evidence = {
        "pytest_node": request.node.nodeid,
        "invocation": list(request.config.invocation_params.args),
        "python": sys.executable, "workers": workers, "batch_size": batch_size,
        "interrupt_after_remove": interrupt, "budget": BUDGET,
        "v1_source_sha256": digest_v1, "events": [], "run_command": command,
        "environment_isolated": all(str(tmp_path) in cli.env[key] for key in (
            "JARVIS_SAMPLERS_HOME", "JARVIS_SCAN_REF_REGISTRY", "JARVIS_REDIS_NOTICE_STATE",
            "JHS_TEST_CATALOG_CACHE", "JHS_TEST_MARKER", "JARVIS_OPERAS_SOURCES_FILE",
        )),
        "home_unchanged": cli.env.get("HOME") == os.environ.get("HOME"),
        "pythonpath_injected": "PYTHONPATH" in cli.env,
    }
    assert evidence["environment_isolated"] and evidence["home_unchanged"]
    assert not evidence["pythonpath_injected"]
    evidence_path = tmp_path / "active-registry-evidence.json"
    started = time.monotonic()

    def event(name, **details):
        evidence["events"].append({"name": name, "elapsed_seconds": round(time.monotonic() - started, 3),
                                   **details})
        evidence_path.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")

    with log.open("w", encoding="utf-8") as handle:
        process = subprocess.Popen(command, cwd=cli.root, env=cli.env, stdout=handle,
                                   stderr=subprocess.STDOUT)
        try:
            prefix = _wait_for_durable_gate(process, database, entered, log)
            binding = json.loads(identity.read_text(encoding="utf-8"))
            assert binding["metadata"]["version"] == "1.0"
            assert binding["source_sha256"] == digest_v1
            bound_source = project / binding["source_file"]
            assert hashlib.sha256(bound_source.read_bytes()).hexdigest() == digest_v1
            original_identity = identity.read_bytes()
            prefix_rows = _assert_rows(prefix, binding, len(prefix))
            event("live-durable-prefix", pid=process.pid, count=len(prefix_rows))

            digest_v2 = _active_source(source, "2.0", .75)
            cli("reload-during-run", "sampler", "load", source)
            assert process.poll() is None
            assert "2.0" in cli("info-during-run", "sampler", "info", "MySampler")
            store = json.loads((cli.registry / "sources.json").read_text(encoding="utf-8"))
            assert store["entries"]["MySampler"]["source_sha256"] == digest_v2
            event("reloaded-while-active", pid_alive=True, v2_source_sha256=digest_v2)

            cli("remove-during-run", "sampler", "remove", "MySampler", "--yes")
            assert process.poll() is None
            assert source.exists() and hashlib.sha256(source.read_bytes()).hexdigest() == digest_v2
            store = json.loads((cli.registry / "sources.json").read_text(encoding="utf-8"))
            assert "MySampler" not in store["entries"]
            assert identity.read_bytes() == original_identity
            assert hashlib.sha256(bound_source.read_bytes()).hexdigest() == digest_v1
            event("removed-while-active", pid_alive=True, source_preserved=True)

            refused = _card(project, "unregistered", port)
            imports = cli.marker.read_bytes()
            cli("new-run-after-remove", "run", refused, "--skip-draw-flowchart", expected=2)
            assert cli.marker.read_bytes() == imports
            assert not (project / "outputs/unregistered").exists()
            assert process.poll() is None
            event("new-run-refused", exit_code=2, existing_scan_alive=True)

            if interrupt:
                process.send_signal(signal.SIGINT)
                # The plug-in's deliberate wait must finish before the runtime
                # can reach its normal interrupt safe point.
                release.touch()
                code = process.wait(timeout=60)
                assert code == 130, log.read_text(encoding="utf-8")[-8000:]
                _assert_port_closed(port)
                persisted = _rows(database)
                assert 0 < len(persisted) < BUDGET
                before_resume = _assert_rows(persisted, binding, len(persisted))
                assert before_resume[:len(prefix_rows)] == prefix_rows
                assert checkpoint.exists()
                event("interrupted-after-remove", exit_code=code, durable_rows=len(persisted),
                      redis_port_closed=True,
                      checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest())
                cli("resume-without-registration", "run", card, "--resume", "--skip-draw-flowchart")
            else:
                release.touch()
                code = process.wait(timeout=60)
                assert code == 0, log.read_text(encoding="utf-8")[-8000:]
                before_resume = prefix_rows
                _assert_port_closed(port)
                event("completed-after-remove", exit_code=code, redis_port_closed=True)

            _assert_port_closed(port)
            result = _assert_rows(_rows(database), binding, BUDGET)
            assert result[:len(before_resume)] == before_resume
            assert identity.read_bytes() == original_identity
            assert hashlib.sha256(bound_source.read_bytes()).hexdigest() == digest_v1
            assert hashlib.sha256(source.read_bytes()).hexdigest() == digest_v2
            evidence["archived_records"] = result
            evidence["artifact_sha256"] = {
                str(path.relative_to(tmp_path)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (source, card, identity, bound_source, database, checkpoint)
                if path.exists()
            }
            event("verified", archived=BUDGET, unique_uuids=BUDGET,
                  original_durable_prefix_preserved=True, source_version="1.0",
                  redis_port_closed=True)
        finally:
            release.touch()
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=10)
            evidence["initial_run_exit_code"] = process.returncode
            evidence["duration_seconds"] = round(time.monotonic() - started, 3)
            evidence["cli_logs"] = sorted(str(path.relative_to(tmp_path)) for path in cli.logs.glob("*.log"))
            evidence_path.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
