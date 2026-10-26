"""Unit execution identity at the scan/restore boundary (no unit arithmetic)."""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Mapping
from pathlib import Path

from jarvishep2.operas_units import (
    DATABASE_KEY,
    UnitBridgeError,
    digest,
    execution_identity,
)


def _refuse(reason):
    raise UnitBridgeError(
        reason + "; use a new scan name or explicitly start fresh",
        code="JV2-UNIT-004",
        path="$.Operas",
    )


def preflight_units(core, *, resume=False, check=False):
    """Must precede engine restore, HDF5 recovery and reuse of persisted rows."""
    config = core.config
    if not isinstance(config, Mapping):
        return
    descriptor = execution_identity(config)
    if descriptor is None:
        config.pop("_unit_execution", None)
    directory = (
        Path(
            core.info.get("task_result_dir")
            or config.get("task_result_dir")
            or os.getcwd()
        )
        / "DATABASE"
    )
    path = directory / "unit_execution.json"
    old = None
    if not check and path.exists():
        try:
            old = json.loads(path.read_text())
            if old["schema_version"] != 1 or old["mode"] not in {"legacy", "units"}:
                raise ValueError("unknown manifest schema/mode")
            if (
                old["mode"] == "units"
                and digest(old["identity_payload"]) != old["execution_identity"]
            ):
                raise ValueError("manifest digest mismatch")
        except (KeyError, ValueError, TypeError) as exc:
            if resume:
                _refuse(f"unit manifest cannot be verified: {exc}")
            old = {"mode": "unverifiable"}  # Explicit fresh replaces the marker.
    if descriptor is None and old is None:
        if resume and not check:
            # Removing a unit declaration together with a lost manifest must
            # not turn unit rows into ordinary completed results.
            from jarvishep2.io.database import _open_hdf5_readonly

            for shard in directory.glob("samples*.hdf5"):
                try:
                    with _open_hdf5_readonly(str(shard)) as handle:
                        if "records" in handle:
                            records = handle["records"]
                            for start in range(0, len(records), 1024):
                                if any(
                                    DATABASE_KEY in json.loads(raw)
                                    for raw in records[start : start + 1024]
                                ):
                                    _refuse(
                                        "DATABASE contains unit records but its unit manifest is missing"
                                    )
                except UnitBridgeError:
                    raise
                except OSError:
                    # Unknown, non-unit historical files retain their existing
                    # recovery path. The post-recovery check runs before reuse.
                    continue
        return
    core._unit_descriptor = descriptor
    if descriptor:
        config["_unit_execution"] = descriptor
    else:
        config.pop("_unit_execution", None)
    if check:
        return
    if resume:
        current = descriptor["execution_identity"] if descriptor else None
        previous = old.get("execution_identity") if old else None
        if current != previous:
            _refuse("unit policy, bindings, contract, JO version or activation changed")
        method = str((config.get("Sampling") or {}).get("Method", "")).lower()
        if current and method not in {"grid", "random", "bridson", "csv"}:
            _refuse(
                f"unit resume for {method!r} is not supported; fresh unit execution is supported"
            )
        checkpoint = core.checkpoint_file()
        shards = sorted(directory.glob("samples*.hdf5"))
        if current and shards and not os.path.isfile(checkpoint):
            _refuse("unit DATABASE exists without a verifiable checkpoint")
        prefix = None
        if os.path.isfile(checkpoint):
            try:
                from jarvishep2.sampling.runtime_checkpoint import load_checkpoint

                payload = load_checkpoint(checkpoint)
                if (payload.get("integrity") or {}).get("unit_execution") != current:
                    _refuse("checkpoint and unit manifest disagree")
                if (payload.get("persistence") or {}).get("replacing_existing"):
                    state = payload.get("sampler_state") or {}
                    prefix = int(
                        state.get("checkpoint_index", state.get("accepted_index", 0))
                        or 0
                    )
            except UnitBridgeError:
                raise
            except Exception as exc:
                _refuse(f"unit checkpoint cannot be verified: {exc}")
        # Read only. Do not run HDF5 recovery before deciding whether these
        # rows are even eligible for this scan's semantic identity.
        from jarvishep2.io.database import _open_hdf5_readonly

        for shard in shards:
            try:
                with _open_hdf5_readonly(str(shard)) as handle:
                    if "records" not in handle:
                        continue
                    records = handle["records"]
                    for start in range(0, len(records), 1024):
                        for raw in records[start : start + 1024]:
                            row = json.loads(raw)
                            if (
                                prefix is not None
                                and int(row.get("sample_index", prefix)) >= prefix
                            ):
                                continue
                            units = row.get(DATABASE_KEY)
                            if current:
                                if (
                                    not isinstance(units, dict)
                                    or units.get("schema_version") != 1
                                    or units.get("execution_identity") != current
                                ):
                                    _refuse(
                                        f"{shard.name} contains unverified/mixed unit records"
                                    )
                            elif units is not None:
                                _refuse(
                                    "unit records cannot be reused by a legacy scan"
                                )
            except UnitBridgeError:
                raise
            except Exception as exc:
                _refuse(f"unit DATABASE cannot be verified read-only: {exc}")
    core._unit_manifest_path = str(path)
    # An explicit fresh transition to legacy retains a marker so an interrupted
    # replacing_existing run does not mistake its historical tail for new data.
    core._unit_manifest_value = descriptor or dict(
        schema_version=1, mode="legacy", execution_identity=None
    )


def publish_unit_manifest(core):
    path = getattr(core, "_unit_manifest_path", None)
    if (
        not isinstance(path, str)
        or not path
        or core.config.get("_check_modules_sample_layout")
    ):
        return
    value = core._unit_manifest_value
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".units-", dir=os.path.dirname(path))
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(value, handle, sort_keys=True, allow_nan=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        if os.name == "posix":
            directory_fd = os.open(os.path.dirname(path), os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
