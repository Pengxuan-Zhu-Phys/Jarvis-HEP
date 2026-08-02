#!/usr/bin/env python3
"""Lightweight Sample model for Redis transport and Worker-side reconstruction."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import os
from typing import Any, Mapping, Protocol, runtime_checkable
from uuid import uuid4

import numpy as np

from jarvishep2.runtime_config import should_eager_materialize, should_materialize_on_failure
from jarvishep2.sample_logger import BufferedSampleLogger, SampleLogger


def _make_lazy_sample_logger(*, module: str) -> BufferedSampleLogger:
    return BufferedSampleLogger(
        extra={
            "module": module,
            "to_console": True,
            "Jarvis": True,
            "_log_domain": "jarvis_hep",
        }
    )


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


VALID_EXECUTION_STEP_TYPES = frozenset(
    {"calculator", "opera", "likelihood", "nuisance_optimize"}
)

_INFO_PROJECTION_KEYS = (
    "save_dir",
    "run_log",
    "logger_name",
    "pack_id",
    "NAttempt",
    "worker_id",
    "staging_path",
    "product_list",
    "bucket_id",
    "bucket_dir",
    "bucket_name",
    "error",
    "error_type",
    "failed_module",
)
_TIMING_KEYS = ("elapsed_s", "started_at", "finished_at")


@runtime_checkable
class UMapperProtocol(Protocol):
    """Worker-held u → x mapper contract (full UMapper implementation lands in D1)."""

    def map(self, u_coords: np.ndarray) -> Mapping[str, Any]:
        """Map normalized u_coords to physical parameter dict."""
        ...


def _validate_execution_step_type(step_type: str) -> str:
    normalized = str(step_type).strip()
    if normalized not in VALID_EXECUTION_STEP_TYPES:
        allowed = ", ".join(sorted(VALID_EXECUTION_STEP_TYPES))
        raise ValueError(f"invalid execution step type '{step_type}'; allowed: {allowed}")
    return normalized


def _as_float64_array(value: Any) -> np.ndarray:
    if value is None:
        return np.array([], dtype=np.float64)
    if isinstance(value, np.ndarray):
        return np.asarray(value, dtype=np.float64)
    if isinstance(value, (list, tuple)):
        return np.asarray(value, dtype=np.float64)
    return np.asarray([value], dtype=np.float64)


@dataclass
class ExecutionStep:
    """One workflow step in the execution plan (JSON-serializable)."""

    type: str
    name: str
    layer: int
    params: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "type": self.type,
            "name": self.name,
            "layer": self.layer,
            "params": dict(self.params),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ExecutionStep:
        return cls(
            type=_validate_execution_step_type(data.get("type", "")),
            name=str(data.get("name", "")),
            layer=int(data.get("layer", 0)),
            params=dict(data.get("params") or {}),
        )


@dataclass
class Sample:
    """Unit of work: identity + u_coords cross Redis; everything else is Worker-local."""

    uuid: str
    u_coords: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    execution_plan: list[ExecutionStep] = field(default_factory=list)
    opera_params: dict[str, Any] = field(default_factory=dict)
    sample_artifacts: str = "auto"
    priority: int = 0
    created_at: str | None = None
    # Monotonic logical proposal index.  Keep this after the historical wire
    # fields so positional ``Sample(uuid, u_coords, ...)`` callers retain their
    # established meaning.
    sample_index: int | None = None
    params: dict[str, Any] = field(default_factory=dict, repr=False)
    info: dict[str, Any] = field(default_factory=dict, repr=False)
    observables: dict[str, Any] = field(default_factory=dict, repr=False)
    status: str = "Created"
    _materialized: bool = field(default=False, repr=False)
    _logger: SampleLogger | BufferedSampleLogger | None = field(default=None, repr=False)
    _with_nuisance: bool = field(default=False, repr=False)
    _likelihood: Any = field(default=None, repr=False)
    _params_bound: bool = field(default=False, repr=False)

    @classmethod
    def from_params(cls, params: Mapping[str, Any]) -> Sample:
        """Build a control-side Sample from a parameter dict (V1-compatible seed)."""
        seed = dict(params)
        sample = cls(
            uuid=str(uuid4()),
            u_coords=np.array([], dtype=np.float64),
            params=seed,
            observables=dict(seed),
        )
        sample.observables["uuid"] = sample.uuid
        return sample

    @property
    def u(self) -> np.ndarray:
        return self.u_coords

    @property
    def likelihood(self) -> Any:
        return self._likelihood

    @likelihood.setter
    def likelihood(self, value: Any) -> None:
        self._likelihood = value

    @property
    def save_dir(self) -> str | None:
        save_dir = self.info.get("save_dir")
        return str(save_dir) if save_dir else None

    def update_uuid(self, new_uuid: str) -> None:
        self.uuid = str(new_uuid)
        self.observables["uuid"] = self.uuid
        self._mirror_fields_to_info()

    def to_task_dict(self) -> dict[str, Any]:
        """Light dict for Redis transport (invariant #7, #8)."""
        payload = {
            "uuid": self.uuid,
            "sample_index": self.sample_index,
            "u_coords": self.u_coords.tolist(),
            "execution_plan": [step.to_dict() for step in self.execution_plan],
            "opera_params": dict(self.opera_params),
            "sample_artifacts": self.sample_artifacts,
            "priority": int(self.priority),
            "created_at": self.created_at or _utc_now_iso(),
        }
        for forbidden in ("logger", "handlers", "params", "info", "observables"):
            if forbidden in payload:
                raise ValueError(f"forbidden key on wire: {forbidden}")
        return payload

    @classmethod
    def from_task_dict(cls, data: Mapping[str, Any]) -> Sample:
        """Reconstruct inside a Worker; logger and paths are not restored."""
        plan_raw = data.get("execution_plan") or []
        execution_plan = [
            step if isinstance(step, ExecutionStep) else ExecutionStep.from_dict(step)
            for step in plan_raw
        ]
        return cls(
            uuid=str(data["uuid"]),
            sample_index=(
                int(data["sample_index"])
                if data.get("sample_index") is not None
                else None
            ),
            u_coords=_as_float64_array(data.get("u_coords")),
            execution_plan=execution_plan,
            opera_params=dict(data.get("opera_params") or {}),
            sample_artifacts=str(data.get("sample_artifacts", "auto")),
            priority=int(data.get("priority", 0)),
            created_at=data.get("created_at"),
            _materialized=False,
            _logger=None,
        )

    def to_info_dict(self) -> dict[str, Any]:
        """Result/monitor projection without live handles (invariant #8)."""
        projection: dict[str, Any] = {
            "uuid": self.uuid,
            "sample_index": self.sample_index,
            "status": self.status,
            "params": dict(self.params),
            "observables": dict(self.observables),
            "likelihood": self._likelihood,
            "opera_params": dict(self.opera_params),
            "sample_artifacts": self.sample_artifacts,
            "execution_plan": [step.to_dict() for step in self.execution_plan],
            "priority": self.priority,
            "created_at": self.created_at,
        }
        if isinstance(self.info, dict):
            timings = {
                key: self.info[key] for key in _TIMING_KEYS if key in self.info
            }
            if timings:
                projection["timings"] = timings
            for key in _INFO_PROJECTION_KEYS:
                if key in self.info:
                    projection[key] = self.info[key]
        projection.pop("logger", None)
        return projection

    def bind_params(self, mapper: UMapperProtocol | None) -> None:
        """Map u_coords → params/observables via the Worker-held UMapper.

        Idempotent: safe to call multiple times; always re-maps from the current
        ``u_coords``. The full ``UMapper`` class lands in a later milestone; any object
        exposing ``map(u_coords) -> Mapping[str, Any]`` is accepted here.
        """
        if mapper is None:
            return
        if not hasattr(mapper, "map"):
            raise TypeError("mapper must provide a map(u_coords) method")

        mapped = mapper.map(self.u_coords)
        if not isinstance(mapped, Mapping):
            raise TypeError("mapper.map(u_coords) must return a mapping")

        self.params = dict(mapped)
        self.observables = dict(mapped)
        self.observables["uuid"] = self.uuid
        self._params_bound = True
        self._mirror_fields_to_info()

    def set_config(self, config: Mapping[str, Any]) -> None:
        self.info = dict(config)
        self.sample_artifacts = str(
            self.info.get("sample_artifacts", self.sample_artifacts)
        ).strip().lower()
        self.create_info()
        if self.info.get("nuisance", {}):
            self.combine_nuisance_card()
            self._with_nuisance = True
        if should_eager_materialize(self.info):
            self.materialize()

    def create_info(self) -> None:
        # Prefer an already-assigned SAMPLE/<bucket> parent (Worker Redis allocator).
        if not isinstance(self.info, dict):
            self.info = {}
        bucket_parent = self.info.get("_bucket_parent") or self.info.get("bucket_dir")
        if not bucket_parent:
            bucket_parent = self.info.get("save_dir")
        if not bucket_parent:
            bucket_parent = self.info.get("sample_dirs")
        sample_root = self._resolve_sample_root(self.info)
        logger_name = f"Sample@{self.uuid}"

        lazy_logger = _make_lazy_sample_logger(module=logger_name)
        self.status = "Init"
        self.info.update(
            {
                "sample_dirs": sample_root,
                "save_dir": None,
                "run_log": None,
                "logger_name": logger_name,
                "logger": lazy_logger,
                "handlers": {},
                "_materialized": False,
                "_bucket_parent": bucket_parent if isinstance(bucket_parent, str) else None,
            }
        )
        # D9.6: dual-truth keys share field objects (params/observables identity).
        self._mirror_fields_to_info()
        self._sync_logger_handles(lazy_logger)

    def combine_nuisance_card(self) -> None:
        card = self.info["nuisance"]
        attempt_uuid = "{}@{}".format(self.uuid, card["NAttempt"])
        active_params = dict(card["active"]["param"] or {})
        active_params.update({"uuid": attempt_uuid})
        self.params.update(active_params)
        self.observables = self.params
        self.info["NAttempt"] = card["NAttempt"]
        self._mirror_fields_to_info()

    def gather_nuisance(self) -> None:
        self.observables["uuid"] = self.uuid
        self._mirror_fields_to_info()

    @staticmethod
    def _resolve_sample_root(info: Mapping[str, Any]) -> str:
        bucket_parent = info.get("save_dir")
        if bucket_parent:
            return os.path.dirname(str(bucket_parent).rstrip(os.sep))

        sample_dirs = info.get("sample_dirs")
        if sample_dirs:
            return str(sample_dirs)

        task_root = info.get("task_result_dir", os.getcwd())
        return os.path.join(task_root, "SAMPLE")

    @staticmethod
    def _buffered_logger(sample_info: Mapping[str, Any] | None) -> BufferedSampleLogger | None:
        if not isinstance(sample_info, Mapping):
            return None
        logger = sample_info.get("logger")
        return logger if isinstance(logger, BufferedSampleLogger) else None

    def _active_logger(self) -> SampleLogger | BufferedSampleLogger | None:
        if isinstance(self.info, dict):
            return self.info.get("logger") or self._logger
        return self._logger

    def _sync_logger_handles(self, logger: SampleLogger | BufferedSampleLogger | None) -> None:
        self._logger = logger
        if isinstance(self.info, dict):
            self.info["logger"] = logger

    def child_logger(self, *, module: str | None = None) -> SampleLogger | BufferedSampleLogger | None:
        """Bind a child logger via logger_name without forcing materialization."""
        base = self._active_logger()
        if base is None:
            return None
        child_module = module or self.info.get("logger_name", f"Sample@{self.uuid}")
        return base.bind(module=child_module)

    def materialize(
        self,
        worker_id: str | None = None,
        *,
        bucket_parent: str | None = None,
        failure_message: str | None = None,
    ) -> str:
        if self._materialized:
            return str(self.info.get("save_dir"))

        if worker_id is not None:
            self.info["worker_id"] = worker_id

        buffered = self._buffered_logger(self.info)

        if bucket_parent is None:
            bucket_parent = self.info.get("_bucket_parent")
        if bucket_parent is None:
            bucket_parent = self._resolve_sample_root(self.info)
            os.makedirs(bucket_parent, exist_ok=True)

        save_dir = os.path.join(str(bucket_parent), self.uuid)
        run_log = os.path.join(save_dir, "Sample_running.log")

        os.makedirs(save_dir, exist_ok=True)
        self.info["save_dir"] = save_dir
        self.info["run_log"] = run_log
        self.info["_materialized"] = True
        self._materialized = True

        if buffered is not None and buffered.event_count > 0:
            buffered.replay_to(run_log)
            buffered.discard()
            self._open_sample_logger(announce_creation=False)
        else:
            self._open_sample_logger(announce_creation=True)

        if failure_message:
            logger = self._active_logger()
            if logger is not None:
                logger.error(failure_message)

        return save_dir

    def _open_sample_logger(self, *, announce_creation: bool = True) -> None:
        logger_name = f"Sample@{self.info['uuid']}"
        self.info["logger_name"] = logger_name

        logger = SampleLogger.open(
            self.info["run_log"],
            module=logger_name,
            extra={
                "to_console": True,
                "Jarvis": True,
                "_log_domain": "jarvis_hep",
            },
        )
        if announce_creation:
            logger.info("Sample created into the Disk")
        self._sync_logger_handles(logger)

    def materialize_failure_artifacts(self, error: Exception | str | None = None) -> str | None:
        if not should_materialize_on_failure(self.info):
            return None
        if self._materialized and self.info.get("save_dir"):
            return str(self.info["save_dir"])

        message = None
        if error is not None:
            message = f"Sample failed -> {error}"
        return self.materialize(failure_message=message)

    def _mirror_fields_to_info(self) -> None:
        """Project field SSOT into ``info`` (D9.6).

        ``params`` / ``observables`` are **shared object identity** with
        ``info["params"]`` / ``info["observables"]`` so in-place mutators used by
        Likelihood / Calculators stay consistent. ``status`` / ``uuid`` are
        mirrored as values.
        """
        if not isinstance(self.info, dict):
            self.info = {}
        self.info["uuid"] = self.uuid
        self.info["status"] = self.status
        self.info["params"] = self.params
        self.info["observables"] = self.observables

    def adopt_params(self, params: Mapping[str, Any], *, as_observables: bool = True) -> None:
        """Replace physical params (and optionally seed observables) then mirror."""
        self.params = dict(params)
        if as_observables:
            self.observables = dict(params)
            self.observables["uuid"] = self.uuid
        self._mirror_fields_to_info()

    def merge_observables(self, updates: Mapping[str, Any]) -> None:
        """Merge calculator/opera outputs into observables (single write entry)."""
        if not isinstance(updates, Mapping):
            raise TypeError("merge_observables requires a mapping")
        self.observables.update(dict(updates))
        self._mirror_fields_to_info()

    def set_status(self, status: str) -> None:
        """Update status field and its info projection."""
        self.status = str(status)
        self._mirror_fields_to_info()

    def record_pack_id(self, step_name: str, pack_id: str) -> None:
        """Record a calculator pack_id for the given execution-plan step."""
        if not isinstance(self.info, dict):
            self.info = {}
        pack_ids = dict(self.info.get("pack_ids") or {})
        pack_ids[str(step_name)] = str(pack_id)
        self.info["pack_ids"] = pack_ids
        self.info["pack_id"] = str(pack_id)

    def record_failure(
        self,
        error: BaseException | str,
        *,
        failed_module: str | None = None,
    ) -> None:
        """Stamp Failed status + machine-readable error fields (D11.1 / D9.6)."""
        self.set_status("Failed")
        if not isinstance(self.info, dict):
            self.info = {}
        self.info["error"] = str(error)
        if isinstance(error, BaseException):
            self.info["error_type"] = type(error).__name__
        else:
            self.info["error_type"] = "Error"
        if failed_module:
            self.info["failed_module"] = str(failed_module)
        elif not self.info.get("failed_module") and self.execution_plan:
            last = self.execution_plan[-1]
            self.info["failed_module"] = str(getattr(last, "name", "") or "")

    def record_handoff(
        self,
        *,
        staging_path: str,
        product_list: list[Any] | None = None,
        clear_save_dir: bool = True,
    ) -> None:
        """Record staging handoff paths after SAMPLE → staging move."""
        if not isinstance(self.info, dict):
            self.info = {}
        self.info["staging_path"] = str(staging_path)
        self.info["product_list"] = list(product_list or [])
        if clear_save_dir:
            self.info["save_dir"] = None

    def pull_dual_truth_from_info(self) -> None:
        """Adopt info dual-truth keys into fields after external in-place mutators.

        Used after Likelihood/Calculator paths that receive ``sample.info`` and may
        rebind ``info['observables']`` / ``info['params']`` to new dicts.
        """
        if not isinstance(self.info, dict):
            return
        obs = self.info.get("observables")
        if isinstance(obs, dict):
            self.observables = obs
        params = self.info.get("params")
        if isinstance(params, dict):
            self.params = params
        status = self.info.get("status")
        if status is not None and str(status):
            self.status = str(status)
        uuid_val = self.info.get("uuid")
        if uuid_val is not None and str(uuid_val):
            self.uuid = str(uuid_val)
        self._mirror_fields_to_info()

    def start(self) -> None:
        logger = self._active_logger()
        if logger is not None:
            logger.info("Sample -> {} is ready for submittion".format(self.uuid))
        self.set_status("Running")
        if self._with_nuisance and isinstance(self.info, dict):
            self.info.setdefault("nuisance", {})
            if isinstance(self.info["nuisance"], dict):
                self.info["nuisance"]["status"] = "Running"
            if logger is not None:
                logger.info(
                    "{}\nSample start {}-th nuisance attempt".format(
                        ">" * 60, self.info["nuisance"]["NAttempt"]
                    )
                )

    def close(self) -> None:
        if self._materialized:
            logger = self._active_logger()
            if logger is not None:
                logger.info(self._build_close_message())
        self.close_logger()

    def _build_close_message(self) -> str:
        observables = self.observables if isinstance(self.observables, dict) else {}
        if not observables:
            return "Sample closed"
        summary = self.format_summary(observables)
        return (
            "Sample SUMMARY\n"
            "============================================================================\n"
            f"{summary}\n"
            "============================================================================\n"
            "Sample closed"
        )

    def close_logger(self) -> None:
        logger = self._active_logger()
        if logger is not None:
            if isinstance(logger, BufferedSampleLogger):
                logger.discard()
            else:
                try:
                    logger.close()
                except Exception:
                    pass

        if isinstance(self.info, dict):
            self.info["handlers"] = {}
            self.info["logger"] = None
        self._logger = None

    def record(self) -> dict[str, Any]:
        """Final Archiver record: observables plus result metadata."""
        info = self.to_info_dict()
        metadata = {
            key: info[key]
            for key in (
                "status",
                "params",
                "likelihood",
                "timings",
                "save_dir",
                "pack_id",
                "worker_id",
                "logger_name",
                "sample_artifacts",
                "execution_plan",
                "priority",
                "sample_index",
                "created_at",
            )
            if key in info
        }
        return {
            "uuid": info["uuid"],
            "observables": dict(info.get("observables", {})),
            "metadata": metadata,
        }

    @staticmethod
    def format_summary(
        values: Mapping[str, Any],
        *,
        key_width: int | None = None,
        value_width: int = 60,
        float_precision: int = 6,
    ) -> str:
        def _format_value(value: Any) -> str:
            if value is None:
                rendered = "None"
            elif isinstance(value, np.generic):
                rendered = str(value.item())
            elif isinstance(value, bool):
                rendered = "True" if value else "False"
            elif isinstance(value, int):
                rendered = str(value)
            elif isinstance(value, float):
                if np.isnan(value):
                    rendered = "nan"
                elif np.isposinf(value):
                    rendered = "inf"
                elif np.isneginf(value):
                    rendered = "-inf"
                else:
                    av = abs(value)
                    if (av != 0.0 and av < 1e-4) or av >= 1e6:
                        rendered = f"{value:.{float_precision}e}"
                    else:
                        rendered = f"{value:.{float_precision}f}".rstrip("0").rstrip(".")
            else:
                rendered = str(value)

            if len(rendered) > value_width:
                if value_width <= 3:
                    rendered = rendered[:value_width]
                else:
                    rendered = rendered[: value_width - 3] + "..."
            return rendered

        keys = list(values.keys())
        if key_width is None:
            key_width = max((len(str(key)) for key in keys), default=1)
            key_width = min(key_width, 60)

        lines = []
        for key in keys:
            key_text = str(key)
            if len(key_text) > key_width:
                if key_width <= 3:
                    key_text = key_text[:key_width]
                else:
                    key_text = key_text[: key_width - 3] + "..."
            value_text = _format_value(values[key])
            lines.append(f"{key_text:<{key_width}}  {value_text:>{value_width}}")
        return "\n".join(lines)


def _attach_sample_for_materialize(sample_info: dict[str, Any]) -> Sample:
    seed_params = sample_info.get("params")
    if not isinstance(seed_params, dict) or not seed_params:
        seed_params = sample_info.get("observables", {})
    sample = Sample.from_params(seed_params if isinstance(seed_params, dict) else {})
    sample.info = sample_info
    # Re-bind dual-truth keys so shell Sample and info bag share field identity.
    if isinstance(sample_info.get("params"), dict):
        sample.params = sample_info["params"]
    if isinstance(sample_info.get("observables"), dict):
        sample.observables = sample_info["observables"]
    if sample_info.get("uuid"):
        sample.update_uuid(str(sample_info["uuid"]))
    else:
        sample._mirror_fields_to_info()
    if sample_info.get("status"):
        sample.set_status(str(sample_info["status"]))
    sample._materialized = bool(sample_info.get("_materialized"))
    sample._logger = sample_info.get("logger")
    sample.sample_artifacts = str(sample_info.get("sample_artifacts", "auto"))
    return sample


def ensure_sample_materialized(
    sample_info: dict[str, Any],
    *,
    failure_message: str | None = None,
) -> str | None:
    """Materialize sample artifacts on demand (e.g. @Sdir resolution)."""
    if not isinstance(sample_info, dict):
        return None
    if str(sample_info.get("sample_artifacts", "auto")).strip().lower() == "never":
        return None
    if sample_info.get("_materialized") and sample_info.get("save_dir"):
        return str(sample_info["save_dir"])

    sample = _attach_sample_for_materialize(sample_info)
    return sample.materialize(failure_message=failure_message)


def materialize_failure_artifacts(
    sample_info: dict[str, Any],
    *,
    error: Exception | str | None = None,
) -> str | None:
    """Materialize failed-sample artifacts and replay buffered sample logs."""
    if not isinstance(sample_info, dict):
        return None
    if not should_materialize_on_failure(sample_info):
        return None
    if sample_info.get("_materialized") and sample_info.get("save_dir"):
        return str(sample_info["save_dir"])

    message = None
    if error is not None:
        message = f"Sample failed -> {error}"

    return ensure_sample_materialized(sample_info, failure_message=message)


__all__ = [
    "ExecutionStep",
    "Sample",
    "UMapperProtocol",
    "VALID_EXECUTION_STEP_TYPES",
    "ensure_sample_materialized",
    "materialize_failure_artifacts",
]
