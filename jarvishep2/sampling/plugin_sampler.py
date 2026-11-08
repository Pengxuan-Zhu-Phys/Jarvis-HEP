#!/usr/bin/env python3
"""Run Jarvis-HEP-Sampler plug-ins as ordinary samplers (external sampler S1).

A plug-in only produces unit-cube points (``IndependentSession`` in the SDK).
Everything after that is the same as for ``Random``: mapping, a finite-value
check, ``selection``, accepted indices, deterministic UUIDs, Redis submission
with backpressure, checkpoints and the lifecycle log.

The run preflight (:func:`preflight_plugin`) loads the plug-in and runs its
``prepare`` before Redis, Workers or the Archiver start, and keeps a copy of
the plug-in source next to the scan's checkpoint. Workers never import the
plug-in. Design: ``Jarvis-Books Docs/Jarvis-HEP/specs/sampler-integration.md`` §5.4, §6.2, §7.
"""

from __future__ import annotations

import contextlib
import hashlib
import math
import os
import tempfile
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from jarvishep2.log_kv import PermilleProgress, format_two_column_log
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.sample import Sample
from jarvishep2.sampler_catalog import (
    PluginRecord,
    SamplerSpec,
    plugin_record,
    spec_for,
)
from jarvishep2.sampling.completion_sampler import CompletionSampler
from jarvishep2.sampling.runtime_checkpoint import operas_constants_fingerprint, stable_json_hash
from jarvishep2.sampling.sample_checkpoint import FORMAT, VERSION
from jarvishep2.sampling.stateless_batch import deterministic_sampler_uuid
from jarvishep2.sampling.sampling_utils import (
    BoolConversionError,
    evaluate_selection,
    physical_from_u,
)

IDENTITY_VERSION = 3


def plugin_logger(method: str) -> Any:
    """The ``Jarvis-HEP.Sampler.<Method>`` logger (same as ``lifecycle_log.sampler_logger``)."""
    return get_jarvis_logger(f"sampler.{method.lower()}", module=f"Sampler.{method}")


class JarvisLogBackend:
    """Backend for the SDK ``log`` facade: the sampler logger, Jarvis table layout."""

    def __init__(self, method: str) -> None:
        self.logger = plugin_logger(method)

    def emit(self, level: int, message: str) -> None:
        self.logger.log(level, message)

    def emit_table(self, level: int, title: str, rows: Sequence[tuple[str, Any]]) -> None:
        self.logger.log(level, format_two_column_log(title, rows))


def _sampling(config: Mapping[str, Any]) -> Mapping[str, Any]:
    sampling = config.get("Sampling")
    return sampling if isinstance(sampling, Mapping) else {}


def plugin_identity(
    config: Mapping[str, Any],
    *,
    metadata: Any,
    seed: int,
    point_number: int,
    options: Mapping[str, Any],
    generation_mode: str = "draw",
) -> tuple[str, dict[str, Any]]:
    """Hard scan identity (design §7.3): SHA-256 of the science-defining fields.

    UUIDs derive from it, so changing any of these fields starts a different
    set of samples. Source edits under the same plug-in ``version`` do not
    change it.
    """
    from jarvis_hep_sampler.options import thaw_options
    from jarvis_hep_sampler.streams import STREAM_PROTOCOL

    sampling = _sampling(config)
    fields = {
        "identity_version": IDENTITY_VERSION,
        "method": metadata.method,
        "plugin_version": metadata.version,
        "profile": metadata.profile,
        "api_version": metadata.api_version,
        "state_version": metadata.state_version or 1,
        "stream_protocol": STREAM_PROTOCOL,
        "generation_mode": generation_mode,
        "block_size": metadata.block_size if generation_mode == "draw_block" else None,
        "variables": sampling.get("Variables") or [],
        "mapper": sampling.get("Mapper") or [],
        "selection": sampling.get("selection"),
        "seed": int(seed),
        "point_number": int(point_number),
        "options": thaw_options(options),
    }
    fields.update(
        loglikelihood=sampling.get("LogLikelihood"),
        operas=config.get("Operas"),
        calculators=config.get("Calculators"),
        nuisance=sampling.get("Nuisance"),
        feedback_return=sampling.get("FeedbackReturn"),
    )
    if metadata.profile == "generational":
        fields["feedback_protocol"] = 1
    from jarvishep2.sampling.plugin_identity import canonical_hash

    return canonical_hash(fields), fields


def plugin_source_dir(
    task_root: str, scan_name: str, method: str, *, check: bool = False,
) -> Path:
    """Source beside this run's checkpoint, with a separate ``check`` directory."""
    from jarvishep2.sampling.runtime_checkpoint import checkpoint_path

    state = checkpoint_path(task_root=task_root, scan_name=scan_name,
                            sampler_name=method, check=check)
    return Path(state).parent / "plugin"


def _write_copy(directory: Path, name: str, payload: bytes) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / name
    if target.is_file() and target.read_bytes() == payload:
        return target
    handle = tempfile.NamedTemporaryFile(dir=directory, prefix=f".{name}.", delete=False)
    try:
        with handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(handle.name, target)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(handle.name)
        raise
    return target


@dataclass
class PluginPreflight:
    """Everything the run prepared for one plug-in scan before Redis started."""

    spec: SamplerSpec
    record: PluginRecord
    plugin: Any
    session: Any
    identity: str
    identity_fields: dict[str, Any] = field(repr=False)
    source_copy: Path | None = None
    resumed_from_copy: bool = False
    restored_state: Mapping[str, Any] | None = field(default=None, repr=False)
    trajectory_digest: Any = field(default=None, repr=False)
    checkpoint_integrity: dict[str, Any] = field(default_factory=dict, repr=False)

    @property
    def method(self) -> str:
        return self.spec.name

    def checkpoint_fields(self, index: int, accepted: int) -> dict[str, Any]:
        return {
            "plugin_checkpoint_version": 1,
            "plugin": {
                "method": self.method,
                "profile": self.plugin.metadata.profile,
                "version": self.plugin.metadata.version,
                "source_sha256": self.plugin.source_sha256,
                "identity": self.identity,
            },
            "maxp": self.session.point_number,
            "seed": self.session.seed,
            "selectionexp": self.identity_fields["selection"],
            "index": index,
            "accepted_index": accepted,
            "next_sample_index": accepted,
        }

    def runtime_metadata(self) -> dict[str, Any]:
        """Allowlisted provenance for the source and settings actually prepared."""
        result = {
            "method": self.plugin.metadata.method,
            "version": self.plugin.metadata.version,
            "source_sha256": self.plugin.source_sha256,
            "seed": self.session.seed,
            "point_number": self.session.point_number,
            "supports_checkpoint": self.session.supports_checkpoint,
        }
        if self.plugin.metadata.profile == "generational":
            result.update(profile="generational", max_evaluations=self.session.max_evaluations,
                          feedback_protocol=1)
        return result


def _variable_names(config: Mapping[str, Any]) -> list[str]:
    from jarvishep2.sampling.variables import load_variables

    return [variable.name for variable in load_variables(dict(config))]


def preflight_plugin(
    config: Mapping[str, Any], *, task_root: str | None = None,
    scan_name: str | None = None, resume: bool = False,
    checkpoint_state: Mapping[str, Any] | None = None,
) -> PluginPreflight | None:
    """Resolve exact scan source and reject drift before any run resources start."""
    from jarvis_hep_sampler import GenerationalSession, IndependentSession, SamplerLog, SamplerRegistry
    from jarvis_hep_sampler.errors import PluginLoadError, PluginResumeMismatch
    from jarvis_hep_sampler.loader import load_plugin_source
    from jarvishep2.cardload.task_config import sampling_method
    from jarvishep2.sampling.plugin_identity import (
        BINDING_VERSION, binding_record, check_identity, environment_versions,
        identity_path, read_binding, source_path,
        validate_checkpoint, write_binding,
    )

    # Programmatic callers may provide the layout separately.
    cfg = dict(config)
    if task_root:
        cfg["task_root"] = task_root
    if scan_name:
        cfg["scan_name"] = scan_name
    method = sampling_method(cfg)
    check = bool(cfg.get("_check_modules_sample_layout"))
    # Calculator smoke is a new session on every invocation, even when its
    # DATABASE/test rows remain. It never resumes or rebinds the formal scan.
    binding = read_binding(cfg) if resume and not check else None
    record = binding_record(cfg) if binding else plugin_record(method)
    spec = (SamplerSpec(name=method, source=record.source, family="plugin",
                        stateless=record.profile == "independent", statistical=record.statistical)
            if record else spec_for(method))
    if spec is None or not spec.is_plugin:
        return None
    if record is None:
        raise PluginLoadError(f"{method} is not a loaded sampler", method=method)
    path = identity_path(cfg)
    source_dir = (plugin_source_dir(task_root, scan_name, method, check=check)
                  if task_root and scan_name else None)
    database = path.parent
    if source_dir is not None and not check:
        has_data = any(p.is_file() and p.stat().st_size for pattern in ("*.hdf5", "*.h5", "*.csv")
                       for p in database.glob(pattern))
        checkpoint = source_dir.parent / "runtime.pickle"
        legacy = source_dir.parent / "state.pkl"
        if binding is None and (has_data or path.exists() or checkpoint.exists() or legacy.exists()):
            guidance = (
                "use --resume with the bound source or choose a new scan name"
                if checkpoint.exists() and path.exists()
                else "a valid runtime.pickle and source binding are required for resume. Choose a new scan name"
            )
            raise PluginResumeMismatch(f"existing scan cannot restart in place. {guidance}", method=method)
    if binding:
        from jarvis_hep_sampler.metadata import parse_metadata
        from jarvis_hep_sampler.options import freeze_options, parse_reserved, split_bounds
        meta = parse_metadata(binding["metadata"])
        reserved, options = split_bounds(_sampling(cfg).get("Bounds") or {}, meta.profile)
        settings = parse_reserved(reserved, meta.profile)
        identity, fields = plugin_identity(
            cfg, metadata=meta, seed=settings.seed, point_number=settings.point_number or settings.max_evaluations,
            options=freeze_options(options), generation_mode=binding["identity_fields"]["generation_mode"],
        )
        check_identity(binding["identity_fields"], fields, method=method)
        if checkpoint_state is not None:
            validate_checkpoint(checkpoint_state, identity=identity, source=binding["source_sha256"],
                                fields=fields, method=method)
        origin = source_path(cfg, binding)
        if origin.exists():
            payload = origin.read_bytes()
        else:
            # Only a byte-identical registered snapshot may repair a missing copy.
            from jarvishep2.sampler_catalog import registered_plugin_record
            current = SamplerRegistry().get(method)
            live = registered_plugin_record(method)
            candidate = (SamplerRegistry().snapshot_path(current) if current is not None
                         else Path(live.snapshot_path) if live and live.source != "scan" else None)
            if candidate is None or not candidate.is_file():
                raise PluginResumeMismatch("bound source is missing and no identical snapshot is available",
                                           method=method, source=str(origin))
            payload = candidate.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if digest != binding["source_sha256"]:
            raise PluginResumeMismatch("bound source SHA-256 mismatch; refusing resume",
                                       method=method, source=str(origin))
    else:
        origin = Path(record.snapshot_path)
        payload = origin.read_bytes()
        if hashlib.sha256(payload).hexdigest() != record.source_sha256:
            raise PluginLoadError("registered source SHA-256 mismatch", method=method, source=str(origin))
    plugin = load_plugin_source(payload, origin=origin)
    session_class = IndependentSession if plugin.metadata.profile == "independent" else GenerationalSession
    if plugin.metadata.profile == "generational":
        sampling = _sampling(cfg)
        if not sampling.get("LogLikelihood"):
            raise PluginLoadError("generational samplers require Sampling.LogLikelihood", method=method)
        if (sampling.get("FeedbackReturn") or {}).get("include_logl") is False:
            raise PluginLoadError("generational samplers require FeedbackReturn.include_logl: true", method=method)
    session = session_class.from_bounds(
        plugin, _sampling(cfg).get("Bounds") or {}, variable_names=_variable_names(cfg),
        log=SamplerLog(method, JarvisLogBackend(method)),
    )
    identity, fields = plugin_identity(
        cfg, metadata=plugin.metadata, seed=session.seed, point_number=session.point_number,
        options=session.options, generation_mode=("ask_tell" if plugin.metadata.profile == "generational" else
                         "draw_block" if plugin.function("draw_block") else "draw"),
    )
    if binding:
        check_identity(binding["identity_fields"], fields, method=method)
    if checkpoint_state is not None:
        validate_checkpoint(checkpoint_state, identity=identity, source=plugin.source_sha256,
                            fields=fields, method=method)
    trajectory_digest = None
    if checkpoint_state is None:
        session.prepare()
    else:
        saved_session = checkpoint_state.get("session")
        if not isinstance(saved_session, dict):
            raise PluginResumeMismatch("checkpoint has no SDK session", method=method)
        if plugin.metadata.profile == "generational":
            from jarvishep2.sampling.generational_plugin_sampler import validate_trajectory

            trajectory_digest = validate_trajectory(path.parent / "chain.csv", checkpoint_state, method=method)
        else:
            try:
                CompletionSampler.validate_completion_state(checkpoint_state)
                if checkpoint_state.get("dataset_id") is not None or checkpoint_state.get("source_state") != {}:
                    raise ValueError("invalid plug-in point source")
                for index, entry in checkpoint_state["unfinished"].items():
                    coords = entry.get("u_coords")
                    expected = deterministic_sampler_uuid(
                        prefix=f"plugin:{method}:{identity[:16]}",
                        seed=session.seed, sample_index=index,
                    )
                    if (not isinstance(coords, np.ndarray) or coords.shape != (session.dimension,)
                            or not np.all(np.isfinite(coords) & (coords >= 0) & (coords <= 1))
                            or not isinstance(entry.get("opera_params"), dict)
                            or entry["uuid"] != expected):
                        raise ValueError("invalid unfinished plug-in proposal")
            except (KeyError, TypeError, ValueError) as exc:
                raise PluginResumeMismatch(f"invalid checkpoint completion state: {exc}", method=method) from exc
        session.restore(saved_session)
    constants = operas_constants_fingerprint() if session.supports_checkpoint else None
    integrity = {
        "unit_execution": (cfg.get("_unit_execution") or {}).get("execution_identity"),
        "operas_constants_hash": stable_json_hash(constants) if constants is not None else None,
    }
    logger = plugin_logger(method)
    environment = environment_versions()
    if binding and environment != binding["environment"]:
        logger.warning("%s plug-in: environment changed (%s -> %s); exact cross-version reproduction is not promised",
                       method, binding["environment"], environment)
    copy = None
    if source_dir is not None:
        if binding:
            copy = source_path(cfg, binding)
            if not copy.exists():
                _write_copy(copy.parent, copy.name, payload)
            logger.warning("%s plug-in: resuming with the scan's own copy %s", method, plugin.source_sha256[:12])
        else:
            copy = _write_copy(source_dir, f"{method}-{plugin.source_sha256}.py", payload)
            write_binding(path, {
                "binding_version": BINDING_VERSION, "hard_identity": identity, "identity_fields": fields,
                "source_sha256": plugin.source_sha256,
                "source_file": str(copy.relative_to(Path(task_root).resolve())),
                "metadata": plugin.metadata.to_dict(), "options_schema": plugin.options_schema,
                "hooks": list(plugin.hooks), "original_source": record.source_path,
                "environment": environment,
            })
    return PluginPreflight(spec=spec, record=record, plugin=plugin, session=session,
                           identity=identity, identity_fields=fields, source_copy=copy,
                           resumed_from_copy=bool(binding), restored_state=checkpoint_state,
                           trajectory_digest=trajectory_digest, checkpoint_integrity=integrity)


class PluginCheckpoint:
    """Publish a snapshot from either external sampler adapter."""

    @property
    def supports_checkpoint(self) -> bool:
        return self._session is not None and self._session.supports_checkpoint

    def save_checkpoint(self, *, force: bool = False) -> bool:
        from jarvis_hep_sampler.errors import PluginOutputError

        writer = self._pickle_writer
        if writer is None or not self._session.can_checkpoint:
            return False
        if not force and not writer.requested():
            return False
        if force and not writer.wait():
            return False
        if not writer.available():
            return False
        try:
            state = self.checkpoint_state()
        except (PluginOutputError, OSError) as exc:
            self._logger.warning("Checkpoint capture failed -> %s", exc)
            return False
        return writer.submit({
            "format": FORMAT,
            "version": VERSION,
            "method": self.method.lower(),
            "config_id": self._config_id,
            "integrity": self._preflight.checkpoint_integrity,
            "sampler_state": state,
        })

    def checkpoint_setting(self) -> str:
        from jarvishep2.runtime_config import get_checkpoint_config

        if not self.supports_checkpoint:
            return "unsupported"
        return "enabled" if get_checkpoint_config(self.config)["enabled"] else "disabled"


class PluginSampler(PluginCheckpoint, CompletionSampler):
    """Independent plug-in points with terminal feedback and optional resume."""

    def __init__(self, spec: SamplerSpec) -> None:
        super().__init__()
        self.method = spec.name
        self.uuid_prefix = f"plugin:{spec.name}"
        self._spec = spec
        self._preflight: PluginPreflight | None = None
        self._session: Any = None
        self._identity = ""
        self._maxp = 0
        self._dimensions = 0
        self._submit_progress: PermilleProgress | None = None
        self._logger = plugin_logger(spec.name)
        self._plugin_status_snapshot: dict[str, Any] = {}
        self._plugin_status_next_at = 0.0

    # ------------------------------------------------------------- set-up
    def attach_preflight(self, preflight: PluginPreflight) -> None:
        """Use the plug-in the run preflight already loaded and prepared."""
        if preflight.method != self.method:
            raise ValueError(
                f"preflight is for {preflight.method}, not {self.method}"
            )
        self._preflight = preflight

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        if self._preflight is None:
            # Direct Python use and tests: no run preflight, no scan copy.
            self._preflight = preflight_plugin(self.config)
        preflight = self._preflight
        assert preflight is not None
        self._session = preflight.session
        self._identity = preflight.identity
        self.uuid_prefix = f"plugin:{self.method}:{self._identity[:16]}"
        self._seed = int(self._session.seed)
        self._maxp = int(self._session.point_number)
        self._dimensions = int(self._session.dimension)
        self._index = 0
        self._accepted_index = 0
        self._submit_progress = None
        self._refresh_plugin_status(force=True)

    @property
    def plugin_version(self) -> str:
        return self._preflight.plugin.metadata.version if self._preflight else ""

    @property
    def plugin_source_sha256(self) -> str:
        return self._preflight.plugin.source_sha256 if self._preflight else ""

    # ------------------------------------------------------------ sampling
    def propose_next(self) -> Sample | None:
        from jarvis_hep_sampler.errors import PluginOutputError

        while self._index < self._maxp:
            self._sampling_tick()
            candidate = self._index
            u_coords = self._session.point(candidate)
            self.check_runtime_control()
            self._index += 1
            self._emit_progress()
            physical = physical_from_u(u_coords, self.vars, self._mapper_pipeline)
            for name, value in physical.items():
                if isinstance(value, (int, float, np.number)) and not math.isfinite(float(value)):
                    raise PluginOutputError(
                        f"mapped value {name} = {float(value)!r} is not finite",
                        method=self.method,
                        candidate_index=candidate,
                    )
            if self._selectionexp:
                try:
                    keep = evaluate_selection(
                        self._selectionexp, physical, context=self._expression_context
                    )
                except BoolConversionError as exc:
                    raise ValueError(f"Invalid selection expression: {self._selectionexp}") from exc
                if not keep:
                    continue
            accepted_index = self._accepted_index
            self._accepted_index += 1
            sample = self._build_sample(u_coords)
            sample.sample_index = accepted_index
            sample.uuid = self._uuid_for_accepted_index(accepted_index)
            return self._remember(sample)
        return None

    def _emit_progress(self) -> None:
        """S4: log every new 1‰ of the candidate budget (not during prefix replay)."""
        if bool(getattr(self, "_suppress_submit_progress", False)):
            return
        if self._submit_progress is None:
            self._submit_progress = PermilleProgress(
                self._logger,
                total=self._maxp,
                label=f"{self.method} candidates sampled",
            )
        self._submit_progress.update(self._index)

    def _stream_exhausted(self) -> bool:
        return self._index >= self._maxp

    # --------------------------------------------------------- checkpoints
    def _source_state(self) -> dict[str, Any]:
        return {}

    def _restore_source_state(self, state: Mapping[str, Any]) -> None:
        if state != {}:
            raise ValueError("invalid plug-in point source")

    def checkpoint_state(self) -> dict[str, Any]:
        return {
            **super().checkpoint_state(),
            **self._preflight.checkpoint_fields(self._index, self._accepted_index),
            "session": self._session.checkpoint(),
        }

    def restore_checkpoint(self, state: Mapping[str, Any]) -> None:
        from jarvishep2.sampling.plugin_identity import validate_checkpoint

        validate_checkpoint(state, identity=self._identity, source=self.plugin_source_sha256,
                            fields=self._preflight.identity_fields, method=self.method)
        self._index = state["index"]
        super().restore_checkpoint(state)
        if self._preflight.restored_state is not state:
            self._session.restore(state["session"])

    # ------------------------------------------------------------- logging
    def _stats(self) -> Any:
        return self._session.stats(candidates=self._index, accepted=self._accepted_index)

    def log_settings_rows(self) -> list[tuple[str, Any]]:
        rows: list[tuple[str, Any]] = [
            ("point number", self._maxp),
            ("checkpoint", self.checkpoint_setting()),
            ("plugin", f"{self.method} {self.plugin_version} ({self.plugin_source_sha256[:12]})"),
            *super().log_settings_rows(),
        ]
        from jarvis_hep_sampler.options import thaw_options

        # Explicit card options are always visible. A reporting hook may
        # enrich/replace their display, or add defaults the host cannot infer.
        options = thaw_options(self._session.options)
        extra = dict(self._session.settings_rows())
        for key, value in list(options.items())[:32]:
            rows.append((key, _setting_text(extra.pop(key, value))))
        for key in options:
            extra.pop(key, None)
        room = max(0, 32 - min(len(options), 32))
        rows.extend((key, _setting_text(value)) for key, value in list(extra.items())[:room])
        omitted = len(options) + len(extra) - min(len(options), 32) - min(len(extra), room)
        if omitted:
            rows.append(("options omitted", omitted))
        return rows

    def log_summary_rows(self) -> list[tuple[str, Any]]:
        return [*super().log_summary_rows(), *self._session.summary_rows(self._stats())]

    def log_stop_reason(self) -> str | None:
        return self._session.stop_reason(self._stats())

    def plugin_status(self) -> dict[str, Any]:
        """Return cached metrics; telemetry never executes user code."""
        return dict(self._plugin_status_snapshot)

    def _sampling_tick(self) -> None:
        super()._sampling_tick()
        self._refresh_plugin_status()

    def _refresh_plugin_status(self, *, force: bool = False) -> None:
        """Run only on the sampling thread, never on a lease/monitor thread."""
        now = time.monotonic()
        if not force and now < self._plugin_status_next_at:
            return
        self._plugin_status_snapshot = self._session.status(self._stats())
        self._plugin_status_next_at = time.monotonic() + 1.0


def _setting_text(value: Any) -> Any:
    """Bound a table cell while keeping simple scalar values readable."""
    if value is None or isinstance(value, (bool, int, float)):
        return value
    text = str(value)
    return text if len(text) <= 256 else text[:253] + "..."


def create_plugin_sampler(spec: SamplerSpec) -> PluginSampler:
    """Distributor factory for plug-in methods."""
    if not spec.stateless:
        from jarvishep2.sampling.generational_plugin_sampler import GenerationalPluginSampler
        return GenerationalPluginSampler(spec)
    return PluginSampler(spec)


__all__ = [
    "JarvisLogBackend",
    "PluginPreflight",
    "PluginSampler",
    "create_plugin_sampler",
    "plugin_identity",
    "plugin_logger",
    "plugin_source_dir",
    "preflight_plugin",
]
