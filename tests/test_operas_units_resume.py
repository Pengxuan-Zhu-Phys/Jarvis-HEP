"""Refusal must precede recovery, resource allocation and any DATABASE mutation."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from jarvishep2.io.database import SimpleHDF5Writer
from jarvishep2.operas_units import UnitBridgeError
from jarvishep2.unit_runtime import preflight_units, publish_unit_manifest


@pytest.fixture
def core(tmp_path):
    from jarvis_operas import get_global_operas_registry, load_user_ops

    load_user_ops(
        str(Path(__file__).parent / "fixtures/unit_ops.py"),
        get_global_operas_registry(),
        namespace="hep_unit_fixture",
        refresh_sympy=False,
    )
    config = dict(
        task_result_dir=str(tmp_path),
        Sampling={"Method": "Random"},
        Operas={
            "Modules": [
                dict(
                    name="Limit",
                    operator="hep_unit_fixture.limit",
                    input=[dict(name="mass_GeV", entry="mass", unit="TeV")],
                    output=[dict(name="rate", entry="value", unit="fb")],
                )
            ]
        },
    )
    return SimpleNamespace(
        config=config, info={}, checkpoint_file=lambda: str(tmp_path / "runtime.pickle")
    )


def _checkpoint(core, *, cursor=1):
    from jarvishep2.sampling.sample_checkpoint import (
        FORMAT, VERSION, atomic_pickle, config_identity,
    )

    descriptor = core.config.get("_unit_execution") or {}
    atomic_pickle(
        core.checkpoint_file(),
        dict(
            format=FORMAT, version=VERSION,
            method=core.config["Sampling"]["Method"].lower(),
            config_id=config_identity(core.config),
            sampler_state={"next_sample_index": cursor},
            integrity={"unit_execution": descriptor.get("execution_identity")},
        ),
    )


def _row(core, index=0):
    descriptor = core.config["_unit_execution"]
    return dict(
        uuid=f"u{index}",
        sample_index=index,
        status="Completed",
        rate=38.332,
        _jarvis_units=dict(
            schema_version=1,
            execution_identity=descriptor["execution_identity"],
            observable_units={"/rate": "fb"},
            calls=[],
            failures=[],
        ),
    )


def _database(core, rows):
    path = Path(core.config["task_result_dir"]) / "DATABASE/samples.hdf5"
    SimpleHDF5Writer(str(path)).add_records(rows)
    return path


@pytest.mark.parametrize(
    "change", ["unit", "policy", "binding", "removed", "contract", "jo_version"]
)
def test_semantic_changes_refuse_without_writing(core, change, monkeypatch):
    preflight_units(core)
    publish_unit_manifest(core)
    _checkpoint(core)
    database = _database(core, [_row(core)])
    original = database.read_bytes()
    mod = core.config["Operas"]["Modules"][0]
    if change == "unit":
        mod["output"][0]["unit"] = "pb"
    if change == "policy":
        mod["unit_policy"] = "strict"
    if change == "binding":
        mod["output"][0]["name"] = "other"
    if change == "removed":
        core.config.pop("Operas")
    if change == "contract":
        from jarvis_operas import get_global_operas_registry, OperaFunction

        declaration = get_global_operas_registry().get("hep_unit_fixture.limit")
        metadata = {
            "unit_contract": get_global_operas_registry().info(
                "hep_unit_fixture.limit"
            )["unit_contract"]
        }
        metadata["unit_contract"]["inputs"]["mass_GeV"]["unit"] = "MeV"
        get_global_operas_registry().register(
            OperaFunction(
                namespace="hep_unit_fixture",
                name="limit",
                arity=None,
                return_dtype=None,
                numpy_impl=declaration.numpy_impl,
                metadata=metadata,
            ),
            replace=True,
        )
    if change == "jo_version":
        from importlib import metadata
        from packaging.version import Version

        original_version = metadata.version
        changed_version = f"{Version(original_version('Jarvis-Operas')).major + 1}.0.0"
        monkeypatch.setattr(
            metadata,
            "version",
            lambda name: changed_version if name == "Jarvis-Operas" else original_version(name),
        )
    with pytest.raises(UnitBridgeError, match="JV2-UNIT-004"):
        preflight_units(core, resume=True)
    assert database.read_bytes() == original


@pytest.mark.parametrize(
    "damage",
    [
        "missing_manifest",
        "bad_manifest",
        "missing_checkpoint",
        "missing_checkpoint_without_rows",
        "mixed_row",
        "bad_schema",
    ],
)
def test_missing_or_mixed_identity_refuses(core, damage):
    preflight_units(core)
    publish_unit_manifest(core)
    _checkpoint(core)
    row = _row(core)
    manifest = Path(core._unit_manifest_path)
    if damage == "missing_manifest":
        manifest.unlink()
        core.config.pop("Operas")
    if damage == "bad_manifest":
        manifest.write_text("{}")
    if damage in {"missing_checkpoint", "missing_checkpoint_without_rows"}:
        Path(core.checkpoint_file()).unlink()
    if damage == "mixed_row":
        row["_jarvis_units"]["execution_identity"] = "other"
    if damage == "bad_schema":
        row["_jarvis_units"]["schema_version"] = 2
    has_rows = damage != "missing_checkpoint_without_rows"
    path = _database(core, [row]) if has_rows else None
    before = path.read_bytes() if path else None
    with pytest.raises(UnitBridgeError, match="JV2-UNIT-004"):
        preflight_units(core, resume=True)
    if path:
        assert path.read_bytes() == before


def test_unit_resume_identity_applies_to_every_builtin_sampler(core, subtests):
    from jarvishep2.sampler_catalog import builtin_specs
    from jarvishep2.sampling.sample_checkpoint import atomic_pickle, read_pickle

    preflight_units(core)
    publish_unit_manifest(core)
    _database(core, [_row(core)])
    for spec in builtin_specs():
        with subtests.test(method=spec.name):
            core.config["Sampling"]["Method"] = spec.name
            _checkpoint(core)
            preflight_units(core, resume=True)

            payload = read_pickle(core.checkpoint_file())
            payload["integrity"]["unit_execution"] = "other"
            atomic_pickle(core.checkpoint_file(), payload)
            with pytest.raises(UnitBridgeError, match="checkpoint and unit manifest disagree"):
                preflight_units(core, resume=True)


def test_resume_verifies_rows_above_the_checkpoint_cursor(core):
    # DATABASE is append-only: rows past the cursor stay and are checked too.
    preflight_units(core)
    publish_unit_manifest(core)
    _checkpoint(core, cursor=1)
    _database(
        core, [_row(core), dict(uuid="stale", sample_index=1, status="Completed")]
    )
    with pytest.raises(UnitBridgeError, match="unverified/mixed"):
        preflight_units(core, resume=True)


def test_explicit_fresh_can_replace_unit_policy_and_manifest(core):
    preflight_units(core)
    publish_unit_manifest(core)
    _checkpoint(core)
    core.config.pop("Operas")
    preflight_units(core, resume=False)
    publish_unit_manifest(core)
    assert json.loads(Path(core._unit_manifest_path).read_text())["mode"] == "legacy"


def test_unit_schema_validation_does_not_discover_user_functions(core):
    from unittest.mock import patch
    from jarvishep2.task_validation import validate_task_config

    with patch(
        "jarvishep2.cardload.task_validation._operas_constant_full_names",
        side_effect=AssertionError("discovery"),
    ):
        report = validate_task_config(core.config, check_modules=True)
    assert not any("discovery" in item.message for item in report.issues)
