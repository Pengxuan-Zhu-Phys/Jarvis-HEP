"""Host contract tests: real JO registry, binding, source lifetime and archive."""

from __future__ import annotations

import copy
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from jarvishep2.operas import OperasModule
from jarvishep2.operas_units import (
    DATABASE_KEY,
    UnitBridgeError,
    entry_pointer,
    execution_identity,
    static_issues,
    read_entry,
    OperaResult,
)
from jarvishep2.sample import Sample
from jarvishep2.sample_logger import SampleLogger, BufferedSampleLogger

FIXTURE = Path(__file__).parent / "fixtures/unit_ops.py"


@pytest.fixture(autouse=True)
def registry():
    from jarvis_operas import get_global_operas_registry, load_user_ops

    registry = get_global_operas_registry()
    load_user_ops(
        str(FIXTURE), registry, namespace="hep_unit_fixture", refresh_sympy=False
    )
    return registry


def module_config(**updates):
    cfg = dict(
        name="Limit",
        operator="hep_unit_fixture.limit",
        unit_policy="strict",
        input=[dict(name="mass_GeV", entry="mass", unit="TeV")],
        output=[
            dict(name="limit", entry="value", unit="fb"),
            dict(name="raw", entry=""),
        ],
    )
    cfg.update(updates)
    return cfg


def sample():
    s = Sample.from_params({"mass": 0.2})
    s.enable_units("unit-id")
    return s


@pytest.mark.parametrize(
    "mode,operator", [("call", "limit"), ("acall", "limit"), ("acall", "alimit")]
)
@pytest.mark.parametrize(
    "value,unit", [(200.0, "GeV"), (0.2, "TeV"), (200000.0, "MeV")]
)
def test_equivalent_units_and_original_records(mode, operator, value, unit):
    config = module_config(
        call_mode=mode,
        operator="hep_unit_fixture." + operator,
        input=[dict(name="mass_GeV", entry="mass", unit=unit)],
    )
    original = {"mass": value}
    snapshot = copy.deepcopy(original)
    log = BufferedSampleLogger(extra={"module": "units"})
    result = OperasModule("Limit", config).execute_result(
        original, {}, sample_logger=log
    )
    assert result.values["limit"] == pytest.approx(38.332, rel=1e-13)
    assert result.values["raw"]["mass_received"] == 200.0
    assert result.values["raw"]["value"] == 0.038332
    assert result.values["raw"]["unit"] == "pb"
    assert result.values["raw"]["record_id"] == "canonical-record"
    assert result.values["raw"]["provenance"] == {"dependencies": ["mass-record"]}
    assert original == snapshot and result.units["limit"] == "fb"
    assert log.event_count >= 1


@pytest.mark.parametrize(
    "value,unit,code",
    [
        (0.2, "pb", "dimension_mismatch"),
        (0.2, "gev", "unknown_unit"),
        (True, "GeV", "type_error"),
        (complex(1, 2), "GeV", "type_error"),
        (float("nan"), "GeV", "non_finite"),
        ([1.0, 2.0], "GeV", "shape_mismatch"),
    ],
)
def test_bad_input_never_calls_physics(registry, value, unit, code):
    func = registry.get("hep_unit_fixture.limit").numpy_impl
    before = len(func.__globals__["CALLS"])
    with pytest.raises(UnitBridgeError) as error:
        OperasModule(
            "Limit",
            module_config(input=[dict(name="mass_GeV", entry="mass", unit=unit)]),
        ).execute({"mass": value}, {})
    assert error.value.jo_code == code
    assert len(func.__globals__["CALLS"]) == before


def test_strict_assumption_default_and_missing_entry():
    no_assertion = [dict(name="mass_GeV", entry="mass")]
    with pytest.raises(UnitBridgeError, match="missing_unit"):
        OperasModule("L", module_config(input=no_assertion)).execute(
            {"mass": 200.0}, {}
        )
    result = OperasModule(
        "L", module_config(input=no_assertion, unit_policy="canonical")
    ).execute_result({"mass": 200.0}, {})
    assert result.call["audit"][0]["origin"] == "assumed_canonical"
    defaults = OperasModule("L", module_config(input=[])).execute_result({}, {})
    assert defaults.call["audit"][0]["origin"] == "default_canonical"
    with pytest.raises(UnitBridgeError, match="missing entry"):
        OperasModule("L", module_config()).execute({}, {})


def test_file_logger_and_physical_source_units(tmp_path):
    log_path = tmp_path / "point.log"
    log = SampleLogger.open(str(log_path), module="units")
    try:
        assert copy.deepcopy(log) is log
        cfg = module_config(kwargs={"source_units": "physical-string"})
        OperasModule("L", cfg).execute({"mass": 0.2}, {}, sample_logger=log)
    finally:
        log.close()
    assert "canonical mass 200.0" in log_path.read_text()


def test_sources_alias_expression_conflict_and_overwrite():
    s = sample()
    first = OperasModule(
        "L", module_config(output=[dict(name="stored_rate", entry="value", unit="fb")])
    ).execute_result(s.observables, s.info)
    s.commit_opera_result(first)
    consumer = dict(
        operator="hep_unit_fixture.consume",
        unit_policy="strict",
        input=[dict(name="rate_pb", entry="stored_rate")],
        output=[dict(name="rate", entry="value")],
    )
    out = OperasModule("C", consumer).execute_result(s.observables, s.info)
    assert out.values["rate"] == pytest.approx(0.038332)
    assert out.call["audit"][0]["origin"] == "source_record"
    consumer["input"][0]["unit"] = "pb"
    with pytest.raises(UnitBridgeError, match="source_conflict"):
        OperasModule("C", consumer).execute(s.observables, s.info)
    consumer["input"] = [dict(name="rate_pb", expression="stored_rate")]
    with pytest.raises(UnitBridgeError, match="missing_unit"):
        OperasModule("C", consumer).execute(s.observables, s.info)
    s.observables["stored_rate"] = s.observables["stored_rate"]
    assert "/stored_rate" not in s.info["unit_execution"]["observable_units"]
    assert len(s.info["unit_execution"]["calls"]) == 1


def test_nullable_missing_projection_and_timeout_do_not_commit(registry):
    s = sample()
    cfg = module_config(operator="hep_unit_fixture.nullable")
    assert OperasModule("N", cfg).execute(s.observables, s.info)["limit"] is None
    cfg["operator"] = "hep_unit_fixture.missing"
    with pytest.raises(UnitBridgeError, match="output_path_missing"):
        s.commit_opera_result(
            OperasModule("M", cfg).execute_result(s.observables, s.info)
        )
    cfg = module_config(timeout_sec=0.001, kwargs={"delay": 0.05})
    with pytest.raises(UnitBridgeError, match="timed out"):
        s.commit_opera_result(
            OperasModule("T", cfg).execute_result(s.observables, s.info)
        )
    assert s.info["unit_execution"]["calls"] == []
    assert "limit" not in s.observables


def test_pointer_root_and_arrays():
    config = module_config(
        output=[
            dict(name="special", entry="/special/a~1b~0c.d/1"),
            dict(name="raw", entry=""),
        ]
    )
    out = OperasModule("L", config).execute({"mass": 0.2}, {})
    assert out["special"] == 2.0 and out["raw"]["unit"] == "pb"
    assert entry_pointer("quantities.xs.value") == "/quantities/xs/value"
    with pytest.raises(ValueError):
        entry_pointer("/bad~2escape")
    assert read_entry({"a": None}, "a") is None
    with pytest.raises(KeyError):
        read_entry({}, "a")


@pytest.mark.parametrize(
    "change", ["duplicate", "kwargs", "expression", "empty", "pointer", "reserved"]
)
def test_static_diagnostics(change):
    cfg = module_config()
    if change == "duplicate":
        cfg["output"].append(dict(name="limit", entry="value"))
    if change == "kwargs":
        cfg["kwargs"] = {"mass_GeV": 1.0}
    if change == "expression":
        cfg["input"][0]["expression"] = "mass"
    if change == "empty":
        cfg["input"][0]["entry"] = ""
    if change == "pointer":
        cfg["input"][0]["entry"] = "/bad~2x"
    if change == "reserved":
        cfg["output"][0]["name"] = DATABASE_KEY
    assert any(
        i.level == "error" for i in static_issues({"Operas": {"Modules": [cfg]}})
    )


def test_schema_accepts_only_opera_unit_fields():
    from jarvishep2.task_validation import validate_task_config

    card = dict(
        Scan={"name": "units"},
        EnvReqs={},
        Sampling={
            "Method": "Random",
            "Bounds": {"point_number": 1},
            "Variables": [
                {
                    "name": "mass",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0.1, "max": 0.3},
                    },
                }
            ],
        },
        Operas={"Modules": [module_config()]},
    )
    assert validate_task_config(card).ok
    card["Operas"]["Modules"][0]["input"][0]["unit"] = 1
    assert not validate_task_config(card).ok


@pytest.mark.parametrize("store_samples", [False, True])
def test_all_archive_branches_keep_units(tmp_path, store_samples):
    from jarvishep2.io.archiver import ArchiveProcessor
    from jarvishep2.io.database import (
        StreamingHDF5Writer,
        SimpleHDF5Writer,
        convert_hdf5_to_csv,
    )
    from jarvishep2.runtime_config import get_archiver_config

    path = tmp_path / "samples.hdf5"
    writer = StreamingHDF5Writer(str(path))
    proc = ArchiveProcessor.from_config(
        writer,
        sample_root=str(tmp_path / "SAMPLE"),
        archiver_config=get_archiver_config(
            {"EnvReqs": {"V2": {"store_samples": store_samples}}}
        ),
    )
    s = sample()
    s.commit_opera_result(
        OperasModule("L", module_config()).execute_result(s.observables, s.info)
    )
    try:
        row = s.to_info_dict()
        row["store_samples"] = store_samples
        proc.ingest(row)
        row = dict(
            uuid="failed",
            status="Failed",
            observables={},
            store_samples=store_samples,
            unit_execution=s.info["unit_execution"],
        )
        proc.ingest(row)
        proc.flush_batch(force=True)
    finally:
        writer.close()
    rows = SimpleHDF5Writer(str(path)).read_records()
    assert len(rows) == 2 and all(r[DATABASE_KEY]["schema_version"] == 1 for r in rows)
    assert rows[0][DATABASE_KEY]["observable_units"] == {"/limit": "fb"}
    convert_hdf5_to_csv(str(path))
    exported = list(csv.DictReader(path.with_suffix(".csv").open()))
    assert json.loads(exported[0][DATABASE_KEY]) == rows[0][DATABASE_KEY]


def test_identity_changes_only_with_unit_semantics():
    cfg = {"Operas": {"Modules": [module_config()]}}
    initial = execution_identity(cfg)
    cfg["EnvReqs"] = {"V2": {"workers": 8, "store_samples": False}}
    assert execution_identity(cfg) == initial
    cfg["Operas"]["Modules"][0]["output"][0]["unit"] = "pb"
    assert execution_identity(cfg) != initial


def test_bare_callable_cannot_bypass_units():
    with pytest.raises(UnitBridgeError):
        OperasModule(
            "bare", module_config(operator="jarvishep2.testing.eggbox.eggbox2d_numpy")
        ).preload()


def test_logical_opaque_and_named_mapping():
    with pytest.raises(UnitBridgeError, match="contract_missing"):
        OperasModule(
            "K",
            module_config(
                operator="hep_unit_fixture.logical", kwargs={"observables": {}}
            ),
        ).execute({"mass": 0.2}, {})
    assert (
        OperasModule("L", module_config()).execute(
            {"mass": 0.2, "report": {"a": 1.0}}, {}
        )["limit"]
        > 0
    )


@pytest.mark.parametrize("dtype", [np.float64, np.float32])
def test_array_shape_units_and_input_protection(dtype):
    values = np.array([[-0.2, 0.0, 0.2]], dtype=dtype)
    original = values.copy()
    cfg = module_config(
        operator="hep_unit_fixture.array",
        output=[dict(name="mass_back", entry="value", unit="TeV")],
    )
    result = OperasModule("A", cfg).execute({"mass": values}, {})
    np.testing.assert_allclose(
        result["mass_back"], original, rtol=1e-6 if dtype == np.float32 else 1e-13
    )
    np.testing.assert_array_equal(values, original)


def test_dynamic_output_conflicts_are_atomic():
    s = sample()
    with s.unit_producer("one"):
        s.commit_opera_result(OperaResult({"q": 1.0}, {"q": "GeV"}))
    with s.unit_producer("two"), pytest.raises(UnitBridgeError, match="same-layer"):
        s.commit_opera_result(OperaResult({"new": 2.0, "q": 3.0}, {"q": "TeV"}))
    assert "new" not in s.observables and s.observables["q"] == 1.0
    assert s.info["unit_execution"]["observable_units"] == {"/q": "GeV"}


def test_schema_bundle_and_manual_show_optional_unit():
    from jarvishep2.cli.agent_artifacts import schema_bundle
    from jarvishep2.cli.man import _operas_page

    definitions = schema_bundle()["schemas"][
        "https://jarvis-hep.org/schema/v2/core/common.json"
    ]["$defs"]
    for name in ("operasIO", "operasOutput"):
        assert "unit" in definitions[name]["properties"]
        assert "unit" not in definitions[name]["required"]
    assert any(key["name"] == "unit_policy" for key in _operas_page()["keys"])


def test_reserved_mapper_field_remains_an_archivable_failure():
    s = Sample.from_params({DATABASE_KEY: "forbidden"})
    with pytest.raises(UnitBridgeError) as failure:
        s.enable_units("unit-id")
    s.record_failure(failure.value)
    result = s.to_info_dict()
    assert result["status"] == "Failed" and DATABASE_KEY not in result["observables"]
    assert result["unit_execution"]["failures"][0]["code"] == "JV2-UNIT-005"


def test_whole_bag_rebind_and_mutable_ordinary_stage_clear_sources():
    s = sample()
    s.commit_opera_result(OperaResult({"q": np.array([1.0])}, {"q": "GeV"}))
    s.invalidate_mutable_unit_sources()
    assert not s.info["unit_execution"]["observable_units"]
    s.commit_opera_result(OperaResult({"q": 1.0}, {"q": "GeV"}))
    s.info["observables"] = dict(s.observables)
    s.pull_dual_truth_from_info()
    assert not s.info["unit_execution"]["observable_units"]


def test_direct_worker_blueprint_and_metadata_drift(tmp_path):
    from jarvishep2.runtime.worker_config import build_worker_config
    from jarvishep2.runtime.worker import Worker

    cfg = {"Operas": {"Modules": [module_config()]}}
    blueprint = build_worker_config(
        cfg,
        task_result_dir=str(tmp_path),
        extra={
            "unit_declarations": [{"operator": "wrong"}],
            "file_operation_mode": "inline",
        },
    )
    assert blueprint["unit_execution_identity"]
    worker = Worker(0, {}, blueprint)
    with pytest.raises(UnitBridgeError, match="differ from control"):
        worker._init_runtime()
