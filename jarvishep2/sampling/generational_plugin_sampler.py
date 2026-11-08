"""Host adapter for SDK ask/tell plugins, using real Worker terminal feedback.

Evaluation DATABASE and posterior chain.csv are separate. Checkpoints advance
only after complete tell and chain fsync. A partial generation is
replayed deterministically; any trajectory suffix after the checkpoint is cut.
"""

from __future__ import annotations

import csv
import hashlib
import io
import math
import os
import time
from pathlib import Path

import numpy as np

from jarvis_hep_sampler.errors import PluginOutputError, PluginResumeMismatch
from jarvishep2.feedback_return import (
    LIKELIHOOD_FEEDBACK_VERSION,
    WIRE_LOGL_KEY,
    terminal_log_density,
)
from jarvishep2.log_kv import PermilleProgress, format_two_column_log
from jarvishep2.mapper import MapperPipeline
from jarvishep2.sampling.feedback_sampler import FeedbackSampler
from jarvishep2.sampling.plugin_sampler import PluginCheckpoint, preflight_plugin, plugin_logger
from jarvishep2.sampling.sample_checkpoint import config_identity
from jarvishep2.sampling.sampling_utils import physical_from_u, evaluate_selection
from jarvishep2.sampling.stateless_batch import deterministic_sampler_uuid
from jarvishep2.sampling.variables import load_variables


class GenerationalPluginSampler(PluginCheckpoint, FeedbackSampler):
    """No plugin code executes in Workers or monitor threads."""

    is_generational_plugin = True
    uses_sample_checkpoint = True
    _checkpoint_saved_attributes = frozenset({
        "_accepted_index", "_evaluation_failures", "_generation", "_index",
        "_seed", "_selection_rejected", "_trajectory_offset",
    })
    _checkpoint_excluded_attributes = frozenset({
        "_batch_size", "_committing", "_config_id", "_current_feedback", "_finished",
        "_identity", "_last_feedback", "_last_step_log_at", "_logger", "_on_failure",
        "_pending_uuids", "_pickle_writer", "_plugin_status_next_at", "_plugin_status_snapshot",
        "_preflight", "_proposals", "_seed_seq", "_selectionexp", "_session", "_spec",
        "_submit_progress", "_synthetic", "_trajectory_hash", "_trajectory_initialized",
        "_trajectory_path", "_uuid_to_id", "method", "vars",
    })

    def __init__(self, spec):
        super().__init__()
        self.method, self._spec = spec.name, spec
        self._preflight = self._session = None
        self._identity = ""
        self._logger = plugin_logger(self.method)
        self.vars = []
        self._selectionexp = None
        self._index = self._accepted_index = 0
        self._selection_rejected = self._evaluation_failures = 0
        self._proposals, self._synthetic, self._uuid_to_id = [], [], {}
        self._trajectory_path = None
        self._trajectory_offset = 0
        self._trajectory_hash = hashlib.sha256()
        self._trajectory_initialized = False
        self._last_feedback, self._current_feedback = {}, {}
        self._finished = False
        self._committing = False
        self._plugin_status_snapshot = {}
        self._plugin_status_next_at = 0.0
        self._submit_progress = None
        self._last_step_log_at = None
        self._pickle_writer = None
        self._config_id = ""

    def attach_preflight(self, preflight):
        if preflight.method != self.method:
            raise ValueError("sampler preflight method differs")
        self._preflight = preflight

    def set_config(self, config_info):
        super().set_config(config_info)
        self._preflight = self._preflight or preflight_plugin(self.config)
        self._session = self._preflight.session
        self._identity = self._preflight.identity
        self._seed = self._session.seed
        self.vars = load_variables(self.config)
        self._mapper_pipeline = MapperPipeline.from_config(self.config)
        self._selectionexp = (self.config.get("Sampling") or {}).get("selection")
        self._config_id = config_identity(self.config)
        result_dir = self.config.get("task_result_dir")
        if not result_dir:
            root = Path(self.config.get("task_root") or os.getcwd())
            scan = (
                self.config.get("scan_name")
                or (self.config.get("Scan") or {}).get("name")
                or "scan"
            )
            result_dir = root / "outputs" / scan
        database = Path(result_dir) / "DATABASE"
        if self.config.get("_check_modules_sample_layout"):
            database /= "test"
        self._trajectory_path = database / "chain.csv"
        self._refresh_status()

    @property
    def plugin_version(self):
        return self._preflight.plugin.metadata.version

    @property
    def plugin_source_sha256(self):
        return self._preflight.plugin.source_sha256

    @property
    def _maxp(self):
        return self._session.max_evaluations if self._session else None

    def feedback_return_spec(self):
        return dict(
            self._preflight.plugin.metadata.feedback_return
            or {"mode": "minimal", "include_logl": True}
        )

    def _has_chain(self):
        return self._preflight.plugin.function("transitions") is not None

    def _ensure_trajectory(self):
        if self._trajectory_initialized:
            return
        self._trajectory_path.parent.mkdir(parents=True, exist_ok=True)
        if self._trajectory_offset == 0:
            columns = (
                ["generation", "chain", "step", "source_uuid"]
                + [f"u_{name}" for name in self._session.setup.variable_names]
                + list(self._session.setup.variable_names)
                + ["logL"]
            )
            stream = io.StringIO(newline="")
            csv.writer(stream).writerow(columns)
            payload = stream.getvalue().encode("utf-8")
            with self._trajectory_path.open("wb") as handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            if os.name == "posix":
                directory_fd = os.open(self._trajectory_path.parent, os.O_RDONLY)
                try:
                    os.fsync(directory_fd)
                finally:
                    os.close(directory_fd)
            self._trajectory_hash.update(payload)
            self._trajectory_offset = len(payload)
        self._trajectory_initialized = True

    def propose_generation(self):
        self._proposals = self._session.ask()
        if self._proposals is None:
            self._finished = True
            return None
        self._synthetic, self._uuid_to_id, self._current_feedback = [], {}, {}
        samples = []
        for offset, proposal in enumerate(self._proposals):
            physical = physical_from_u(proposal["u"], self.vars, self._mapper_pipeline)
            for name, value in physical.items():
                if isinstance(value, (int, float, np.number)) and not math.isfinite(
                    float(value)
                ):
                    raise PluginOutputError(
                        f"mapped {name} is not finite", method=self.method
                    )
            if self._selectionexp and not evaluate_selection(
                self._selectionexp, physical, context=self._expression_context
            ):
                self._synthetic.append(
                    {"id": proposal["id"], "status": "selection_rejected", "logL": None}
                )
                continue
            sample = self._build_sample(proposal["u"])
            sample.uuid = deterministic_sampler_uuid(
                prefix=f"plugin:{self.method}:{self._identity}",
                seed=self._seed,
                sample_index=self._session.candidates + offset,
            )
            sample.sample_index = self._accepted_index + len(samples)
            self._uuid_to_id[sample.uuid] = proposal["id"]
            samples.append(sample)
        return samples

    def _await_feedback_round(self, *args, **kwargs):
        record, mode, missing = super()._await_feedback_round(*args, **kwargs)
        if record is not None:
            uuid = str(record.get("uuid", ""))
            previous = self._current_feedback.get(uuid) or self._last_feedback.get(uuid)
            if previous is not None and previous != record:
                raise PluginOutputError(
                    "conflicting duplicate Worker feedback", method=self.method
                )
            if uuid in self._uuid_to_id:
                self._current_feedback[uuid] = dict(record)
        return record, mode, missing

    def absorb_generation(self, results):
        self._committing = True
        terminal = list(self._synthetic)
        failures = 0
        for record in results:
            if (
                record.get(LIKELIHOOD_FEEDBACK_VERSION) != 1
                or record.get("evaluation_outcome") not in {"ok", "failed"}
                or not isinstance(record.get("score_present"), bool)
            ):
                raise PluginOutputError(
                    "Worker lacks valid plugin terminal feedback v1", method=self.method
                )
            status = record["evaluation_outcome"]
            score = terminal_log_density(record)
            if (status == "ok" and score is None) or (
                status == "failed"
                and (record["score_present"] or record.get(WIRE_LOGL_KEY) is not None)
            ):
                raise PluginOutputError(
                    "inconsistent Worker terminal feedback", method=self.method
                )
            proposal_id = self._uuid_to_id.get(str(record["uuid"]))
            if proposal_id is None:
                raise PluginOutputError(
                    "unknown Worker feedback UUID", method=self.method
                )
            failures += status == "failed"
            fields = {
                key: value
                for key, value in record.items()
                if key
                not in {
                    "uuid",
                    WIRE_LOGL_KEY,
                    LIKELIHOOD_FEEDBACK_VERSION,
                    "evaluation_outcome",
                    "score_present",
                }
            }
            terminal.append(
                {"id": proposal_id, "status": status, "logL": score, "fields": fields}
            )
        events = self._session.tell(terminal)
        self._ensure_trajectory()
        stream = io.StringIO(newline="")
        writer = csv.writer(stream)
        id_to_uuid = {value: key for key, value in self._uuid_to_id.items()}
        for event in events:
            physical = self._mapper_pipeline.map(event["u"])
            writer.writerow(
                [
                    self._session.generation - 1,
                    event["chain"],
                    event["step"],
                    id_to_uuid.get(event["source_id"], ""),
                    *event["u"],
                    *(physical[name] for name in self._session.setup.variable_names),
                    event["logL"],
                ]
            )
        payload = stream.getvalue().encode("utf-8")
        if payload:
            with self._trajectory_path.open("ab") as handle:
                handle.write(payload)
                handle.flush()
            self._trajectory_hash.update(payload)
            self._trajectory_offset += len(payload)
        self._accepted_index += len(results)
        self._selection_rejected += len(self._synthetic)
        self._evaluation_failures += failures
        self._index = self._session.candidates
        self._committing = False
        self._refresh_status()
        self._last_feedback = self._current_feedback
        self._submitted_uuids.clear()
        self._proposals, self._synthetic, self._uuid_to_id = [], [], {}
        self._emit_progress()
        self._maybe_log_step()

    def at_safe_barrier(self):
        return self._session.can_checkpoint and not self._committing and not self._pending_uuids

    def configure_checkpoint(self, *, checkpoint_file, writer):
        self._pickle_writer = writer

    def shutdown_checkpointing(self):
        self._pickle_writer = None

    def checkpoint_at_barrier(self, *, reason="generation_barrier", final=False):
        return self.save_checkpoint(force=final)

    def run_adaptive(self, *, generation_timeout=3600.0, timeout=None):
        self._ensure_trajectory()
        self.save_checkpoint(force=True)
        self._emit_progress(force=True)
        self._maybe_log_step(force=True)
        try:
            count = super().run_adaptive(
                generation_timeout=generation_timeout, timeout=timeout
            )
            self.save_checkpoint(force=True)
        finally:
            # Read committed counters only: an unfinished ask/tell never looks
            # like completed progress on failure or interruption.
            self._emit_progress()
            self._maybe_log_step(force=True)
        if self._has_chain():
            self._logger.warning("Chain saved to %s", self._trajectory_path)
        return count

    def _emit_progress(self, *, force=False):
        if self._submit_progress is None:
            self._submit_progress = PermilleProgress(
                self._logger,
                total=self._session.max_evaluations,
                label=f"{self.method} proposals completed",
            )
        extra = f"worker evaluations={self._accepted_index}, generation={self._session.generation}"
        if self._has_chain():
            extra += f", retained={self._session.retained}"
        self._submit_progress.update(self._index, force=force, extra=extra)

    def _maybe_log_step(self, *, force=False):
        now = time.monotonic()
        if self._last_step_log_at is None:
            self._last_step_log_at = now
            if not force:
                return
        if not force and now - self._last_step_log_at < 60.0:
            return
        self._last_step_log_at = now
        rows = self._log_counter_rows()
        saved_index = getattr(self._pickle_writer, "last_saved_index", None)
        if saved_index is not None:
            rows.append(("checkpoint sample index", saved_index))
        rows.extend(self._algorithm_log_rows())
        self._logger.info(format_two_column_log(f"{self.method} Progress ->", rows))

    def _log_counter_rows(self):
        rows = [
            ("generation", self._session.generation),
            ("proposals / limit", f"{self._index} / {self._session.max_evaluations}"),
            ("worker evaluations", self._accepted_index),
            ("selection rejects", self._selection_rejected),
            ("evaluation failures", self._evaluation_failures),
        ]
        if self._has_chain():
            rows.append(("retained draws", self._session.retained))
        return rows

    def _algorithm_log_rows(self):
        from jarvishep2.sampling.plugin_sampler import _setting_text

        return [
            (key, _setting_text(value))
            for key, value in self._session.summary_rows(self._stats())[:32]
        ]

    def checkpoint_state(self):
        if not self.at_safe_barrier():
            raise PluginOutputError("checkpoint requires a complete tell/chain barrier", method=self.method)
        if self._trajectory_initialized:
            with self._trajectory_path.open("rb") as handle:
                os.fsync(handle.fileno())
        return {
            **self._preflight.checkpoint_fields(self._index, self._accepted_index),
            "generation": self._session.generation,
            "session": self._session.checkpoint(),
            "trajectory_offset": self._trajectory_offset,
            "trajectory_sha256": self._trajectory_hash.hexdigest(),
            "selection_rejected": self._selection_rejected,
            "completed_count": self._accepted_index - self._evaluation_failures,
            "failed_count": self._evaluation_failures,
        }

    def restore_checkpoint(self, state):
        from jarvishep2.sampling.plugin_identity import validate_checkpoint

        validate_checkpoint(
            state,
            identity=self._identity,
            source=self.plugin_source_sha256,
            fields=self._preflight.identity_fields,
            method=self.method,
        )
        if self._preflight.restored_state is state:
            self._trajectory_hash = self._preflight.trajectory_digest
        else:
            self._trajectory_hash = validate_trajectory(
                self._trajectory_path, state, method=self.method
            )
            self._session.restore(state["session"])
        self._trajectory_offset = state["trajectory_offset"]
        if self._trajectory_path.exists():
            with self._trajectory_path.open("r+b") as handle:
                handle.truncate(self._trajectory_offset)
                handle.flush()
                os.fsync(handle.fileno())
        self._trajectory_initialized = self._trajectory_offset > 0
        self._index, self._accepted_index = state["index"], state["accepted_index"]
        self._generation = self._session.generation
        self._selection_rejected = state["selection_rejected"]
        self._evaluation_failures = state["failed_count"]
        self._pending_uuids.clear()
        self._submitted_uuids.clear()
        self._submit_progress = None
        self._last_step_log_at = None
        self._refresh_status()
        rows = self._log_counter_rows()
        if self._has_chain():
            rows.append(("chain bytes", self._trajectory_offset))
        self._logger.info(format_two_column_log(f"{self.method} Resume Barrier ->", rows))

    def _stats(self):
        return self._session.stats(self._index, self._session.retained)

    def log_settings_rows(self):
        from jarvishep2.sampling.plugin_sampler import _setting_text

        extra = dict(self._session.settings_rows())
        options = dict(self._session.options)
        rows = [
            ("plugin", f"{self.method} {self.plugin_version} ({self.plugin_source_sha256[:12]})"),
            ("evaluation limit", self._session.max_evaluations),
            ("checkpoint", self.checkpoint_setting()),
        ]
        if self._selectionexp:
            rows.append(("selection", _setting_text(self._selectionexp)))
        merged = {key: extra.pop(key, value) for key, value in options.items()}
        merged.update(extra)
        rows.extend(
            (key, _setting_text(value)) for key, value in list(merged.items())[:32]
        )
        return rows

    def log_summary_rows(self):
        return [*self._log_counter_rows(), *self._algorithm_log_rows()]

    def log_stop_reason(self):
        return self._session.stop_reason(self._stats())

    def _refresh_status(self):
        if time.monotonic() >= self._plugin_status_next_at:
            self._plugin_status_snapshot = self._session.status(self._stats())
            self._plugin_status_next_at = time.monotonic() + 1.0

    def plugin_status(self):
        return {
            "evaluations": self._index,
            "retained": self._session.retained,
            "generation": self._session.generation,
            **self._plugin_status_snapshot,
        }

    def summary(self):
        return {
            "method": self.method,
            "profile": "generational",
            "proposals": self._index,
            "worker_evaluations": self._accepted_index,
            "selection_rejected": self._selection_rejected,
            "evaluation_failures": self._evaluation_failures,
            "retained_draws": self._session.retained,
            "max_evaluations": self._session.max_evaluations,
            "stop_reason": self.log_stop_reason(),
            "chain_file": str(self._trajectory_path),
        }


def validate_trajectory(path, state, *, method):
    """Read-only validation before resources start; returns a hash of the committed prefix."""
    offset = state.get("trajectory_offset")
    digest = hashlib.sha256()
    if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
        raise PluginResumeMismatch("invalid chain byte cursor", method=method)
    if offset:
        if not path.is_file() or path.stat().st_size < offset:
            raise PluginResumeMismatch(
                "chain.csv is missing or shorter than its checkpoint", method=method
            )
        with path.open("rb") as handle:
            remaining = offset
            while remaining:
                chunk = handle.read(min(1048576, remaining))
                digest.update(chunk)
                remaining -= len(chunk)
    if digest.hexdigest() != state.get("trajectory_sha256"):
        raise PluginResumeMismatch(
            "chain.csv committed prefix checksum differs", method=method
        )
    return digest
