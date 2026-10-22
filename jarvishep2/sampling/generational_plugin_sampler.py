"""Host adapter for SDK ask/tell plugins, using real Worker terminal feedback.

Evaluation DATABASE and posterior chain.csv are separate. Checkpoints advance
only after archive acknowledgements and chain fsync. A partial generation is
replayed deterministically; any trajectory suffix after the checkpoint is cut.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import io
import math
import os
import time
from pathlib import Path

import numpy as np

from jarvis_hep_sampler.errors import PluginOutputError, PluginResumeMismatch
from jarvishep2.mapper import MapperPipeline
from jarvishep2.sampling.feedback_sampler import FeedbackSampler
from jarvishep2.sampling.plugin_sampler import preflight_plugin, plugin_logger
from jarvishep2.sampling.sampling_utils import physical_from_u, evaluate_selection
from jarvishep2.sampling.stateless_batch import deterministic_sampler_uuid
from jarvishep2.sampling.variables import load_variables


class GenerationalPluginSampler(FeedbackSampler):
    """No plugin code executes in Workers or monitor threads."""

    is_generational_plugin = True
    _checkpoint_saved_attributes = frozenset(
        {
            "_seed",
            "_generation",
            "_accepted_index",
            "_index",
            "_trajectory_offset",
            "_selection_rejected",
            "_evaluation_failures",
            "_on_failure",
        }
    )
    _checkpoint_excluded_attributes = frozenset(
        {
            "method",
            "_spec",
            "_preflight",
            "_session",
            "_identity",
            "_logger",
            "vars",
            "_selectionexp",
            "_pending_uuids",
            "_seed_seq",
            "_batch_size",
            "_proposals",
            "_synthetic",
            "_uuid_to_id",
            "_trajectory_path",
            "_trajectory_hash",
            "_trajectory_initialized",
            "_last_feedback",
            "_current_feedback",
            "_finished",
            "_committing",
            "_plugin_status_snapshot",
            "_plugin_status_next_at",
        }
    )

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
        self._last_safe_state = copy.deepcopy(self.export_runtime_state())
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

    @property
    def archive_target_count(self):
        return self._accepted_index

    def feedback_return_spec(self):
        return dict(
            self._preflight.plugin.metadata.feedback_return
            or {"mode": "minimal", "include_logl": True}
        )

    def set_persisted_uuids(self, uuids):
        # Replay the last uncommitted generation to reconstruct feedback. The
        # archive deduplicates deterministic UUIDs/sample_index; no lifetime set.
        self._persisted_uuids = set()

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

    def _submit_sample_batch(self, samples):
        # Keep a bounded window across generations until the durable checkpoint.
        for start in range(0, len(samples), self._batch_size):
            batch = list(samples[start : start + self._batch_size])
            for sample in batch:
                self._register_pending(sample.uuid)
                self._submitted_uuids.append(str(sample.uuid))
            self._submit_group(batch)
        return len(samples)

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
                record.get("plugin_feedback_version") != 1
                or record.get("evaluation_outcome") not in {"ok", "failed"}
                or not isinstance(record.get("score_present"), bool)
            ):
                raise PluginOutputError(
                    "Worker lacks valid plugin terminal feedback v1", method=self.method
                )
            status = record["evaluation_outcome"]
            score = record.get("logL")
            if record["score_present"] != (status == "ok") or (
                status == "failed" and score is not None
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
                    "logL",
                    "plugin_feedback_version",
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
        self._last_feedback.update(self._current_feedback)
        self._proposals, self._synthetic, self._uuid_to_id = [], [], {}

    def at_safe_barrier(self):
        return (
            self._session is not None
            and self._session._pending is None
            and not self._submitted_uuids
            and not self._committing
        )

    def checkpoint_runtime_state(self, *, safe):
        # Interrupt callbacks request a live snapshot in older feedback engines.
        # This adapter explicitly refuses partial or archive-unconfirmed state.
        if safe and self.at_safe_barrier():
            self._sync_trajectory()
            self._last_safe_state = copy.deepcopy(self.export_runtime_state())
        return copy.deepcopy(self._last_safe_state)

    def _sync_trajectory(self):
        if self._trajectory_initialized:
            with self._trajectory_path.open("rb") as handle:
                os.fsync(handle.fileno())

    def checkpoint_at_barrier(self, *, reason="generation_barrier", final=False):
        # Bounded windows avoid a disk barrier per cheap Eggbox evaluation.
        if (
            not final
            and len(self._submitted_uuids) < 256
            and not self._checkpoint_requested
        ):
            return False
        scan = (
            self.config.get("scan_name")
            or (self.config.get("Scan") or {}).get("name")
            or "scan"
        )
        if self.redis is not None and self._submitted_uuids:
            deadline = time.monotonic() + 30
            # Tasks have dense sample_index even when selection rejects proposals.
            # The Archiver publishes the fsync-confirmed prefix for indexed rows.
            while self.redis.get_archived_index_prefix(scan) < self._accepted_index:
                self._fail_if_generation_blocked(what="plugin archive barrier")
                if time.monotonic() >= deadline:
                    raise TimeoutError(
                        "archive did not confirm plugin generation window; previous checkpoint preserved"
                    )
                time.sleep(0.05)
        self._sync_trajectory()
        self._submitted_uuids.clear()
        self._persisted_uuids.clear()
        self._last_feedback.clear()
        self._checkpoint_requested = False
        self._last_safe_state = copy.deepcopy(self.export_runtime_state())
        return self.persist_runtime_checkpoint(force=True, reason=reason)

    def run_adaptive(self, *, generation_timeout=3600.0, timeout=None):
        self._ensure_trajectory()
        self.checkpoint_at_barrier(reason="plugin_initial", final=True)
        count = super().run_adaptive(
            generation_timeout=generation_timeout, timeout=timeout
        )
        self.checkpoint_at_barrier(reason="plugin_finished", final=True)
        return count

    def export_runtime_state(self):
        return {
            "plugin": {
                "method": self.method,
                "version": self.plugin_version,
                "source_sha256": self.plugin_source_sha256,
                "identity": self._identity,
            },
            "maxp": self._session.max_evaluations,
            "seed": self._seed,
            "selectionexp": self._selectionexp,
            "index": self._index,
            "accepted_index": self._accepted_index,
            "generation": self._session.generation,
            "session": self._session.checkpoint(),
            "trajectory_offset": self._trajectory_offset,
            "trajectory_sha256": self._trajectory_hash.hexdigest(),
            "selection_rejected": self._selection_rejected,
            "evaluation_failures": self._evaluation_failures,
        }

    def import_runtime_state(self, state):
        from jarvishep2.sampling.plugin_identity import validate_checkpoint

        validate_checkpoint(
            state,
            identity=self._identity,
            source=self.plugin_source_sha256,
            fields=self._preflight.identity_fields,
            method=self.method,
        )
        self._session.restore(state["session"])
        self._trajectory_hash = validate_trajectory(
            self._trajectory_path, state, method=self.method
        )
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
        self._evaluation_failures = state["evaluation_failures"]
        self._pending_uuids.clear()
        self._submitted_uuids.clear()
        self._last_safe_state = copy.deepcopy(dict(state))

    def _stats(self):
        return self._session.stats(self._index, self._session.retained)

    def log_settings_rows(self):
        from jarvishep2.sampling.plugin_sampler import _setting_text

        extra = dict(self._session.settings_rows())
        options = dict(self._session.options)
        rows = [
            ("plugin", f"{self.method} {self.plugin_version}"),
            ("evaluation limit", self._session.max_evaluations),
        ]
        merged = {key: extra.pop(key, value) for key, value in options.items()}
        merged.update(extra)
        rows.extend(
            (key, _setting_text(value)) for key, value in list(merged.items())[:32]
        )
        return rows

    def log_summary_rows(self):
        return [*self.summary().items(), *self._session.summary_rows(self._stats())]

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
