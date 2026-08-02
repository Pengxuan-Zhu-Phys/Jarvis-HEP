#!/usr/bin/env python3
"""Calculator backend held long-term by Workers (sync scheduler path only)."""

from __future__ import annotations

import os
import time
from collections.abc import Mapping
from typing import Any
from jarvishep2.async_subprocess import AsyncSubprocessScheduler, SubprocessJob
from jarvishep2.command_parser import CommandParser
from jarvishep2.expression import ExpressionContext
from jarvishep2.io_portal import (
    apply_hep_io_save,
    build_io_context,
    evaluate_io_expression,
    read_io_output_sync,
    write_io_input_sync,
)
from jarvishep2.Module.calculator_spec import (
    CalculatorSpec,
    is_cwd_switch_command,
    next_cwd_from_command,
)
from jarvishep2.Module.runtime_preparer import RuntimePreparer
from jarvishep2.Sampling.sampling_utils import BoolConversionError, evaluate_selection


class CalculatorModule:
    """External calculator runner: Spec + RuntimePreparer + sync execute."""

    def __init__(self, name: str, config: Mapping[str, Any]) -> None:
        self.spec = CalculatorSpec.from_config(name, config)
        self.name = self.spec.name
        self.config = dict(self.spec.raw)
        self.clone_shadow = self.spec.clone_shadow
        self.installation = list(self.spec.installation)
        self.initialization = list(self.spec.initialization)
        self.execution = dict(self.config.get("execution") or {})
        self.commands_template = list(self.spec.commands)
        self.input_specs = list(self.spec.input_specs)
        self.output_specs = list(self.spec.output_specs)
        self.timeout = self.spec.timeout
        self.basepath = self.spec.basepath
        self.source = self.spec.source
        self.symlink_name = self.spec.symlink_name
        self.env_setup = list(self.spec.env_setup)
        self.selection = self.spec.selection
        self._subprocess_env: dict[str, str] | None = None
        self.sample_info: dict[str, Any] = {}
        self.PackID: str | None = None
        self.logger: Any | None = None
        self._templates_loaded = False
        self._template_parse_count = 0
        self._command_counter = 0
        self._scheduler: AsyncSubprocessScheduler | None = None
        self._command_parser: CommandParser | None = None
        self._expression_context: ExpressionContext | None = None
        self._file_ops: Any | None = None
        force = bool(self.config.get("force_reinstall", False))
        self._preparer = RuntimePreparer(
            self.spec,
            force_reinstall=force,
            install_epoch=int(self.config.get("_install_epoch", 0) or 0),
        )
        # V1: install stage writes to pack-local Installation_{name}-{PackID}.log
        # (under shadow runtime path), not Sample_running.log.
        self._install_logger: Any | None = None
        self._install_log_path: str | None = None

    @staticmethod
    def custom_format(record: Mapping[str, Any]) -> str:
        module = record.get("extra", {}).get("module", "No module")
        if "raw" in record.get("extra", {}):
            return "{message}\n"
        return f"{module} | {record.get('message', '')}"

    def preload_templates(self) -> None:
        if self._templates_loaded:
            return
        _ = list(self.commands_template)
        _ = list(self.input_specs)
        _ = list(self.output_specs)
        self._template_parse_count += 1
        self._templates_loaded = True

    def attach_scheduler(self, scheduler: AsyncSubprocessScheduler | None) -> None:
        self._scheduler = scheduler

    def attach_command_parser(self, parser: CommandParser | None) -> None:
        self._command_parser = parser

    def attach_expression_context(self, context: ExpressionContext | None) -> None:
        self._expression_context = context

    def attach_file_ops(self, file_ops: Any | None) -> None:
        """HEP FileOperation service (save/copy/delete) — not Portal."""
        self._file_ops = file_ops

    def _evaluate_io_expression(self, expression: str, values: Mapping[str, Any]) -> Any:
        return evaluate_io_expression(
            expression,
            values,
            context=self._expression_context,
            logger=self._logger(),
        )

    def bind_env(self, env: Mapping[str, str]) -> None:
        self._subprocess_env = {str(key): str(value) for key, value in env.items()}

    def env_setup_sources(self) -> list[str]:
        sources: list[str] = []
        for item in self.env_setup:
            if isinstance(item, Mapping):
                source = str(item.get("source", "")).strip()
                if source:
                    sources.append(source)
        return sources

    def acquire_pack_id(self, pack_id: str) -> None:
        """Bind a Redis-owned exclusive PackID (must be stable ``001``…, never a UUID)."""
        from jarvishep2.redis_queue import is_stable_calc_pack_id

        text = str(pack_id or "").strip()
        if not is_stable_calc_pack_id(text):
            raise ValueError(
                f"calculator '{self.name}' refused non-stable PackID {pack_id!r}; "
                f"expected Redis free-list slots like '001' (not 'ready' or UUIDs)"
            )
        self.PackID = text
        self._preparer.acquire_pack_id(text)

    def decode_shadow_path(self, path: str) -> str:
        return self._preparer.decode_shadow_path(path)

    def decode_shadow_commands(self, command: Mapping[str, Any]) -> dict[str, str]:
        return self._preparer.decode_shadow_commands(command)

    def shadow_runtime_path(self) -> str:
        return self._preparer.shadow_runtime_path()

    def install_log_path(self) -> str:
        """V1 path: ``{shadow_runtime}/Installation_{name}-{PackID}.log``."""
        pack = str(self.PackID or "NA")
        runtime = self.shadow_runtime_path()
        return os.path.join(runtime, f"Installation_{self.name}-{pack}.log")

    def ensure_install_logger(self) -> Any:
        """Open (or reuse) the pack-local installation log sink."""
        if self._install_logger is not None:
            return self._install_logger
        from jarvishep2.sample_logger import SampleLogger

        pack = str(self.PackID or "NA")
        path = self.install_log_path()
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        label = f"{self.name}-{pack}"
        self._install_logger = SampleLogger.open(path, module=label)
        self._install_log_path = path
        self._install_logger.info(f"Start install {label}")
        return self._install_logger

    def close_install_logger(self) -> None:
        logger = self._install_logger
        self._install_logger = None
        if logger is None:
            return
        closer = getattr(logger, "close", None)
        if callable(closer):
            try:
                closer()
            except Exception:
                pass

    def ensure_shadow_installed(self) -> None:
        try:
            self._preparer.ensure_shadow_installed(
                run_stage=self._run_stage_commands,
                logger=self._stream_logger_for_stage("install"),
            )
        finally:
            # Install is one-shot per pack; release the file handle after the stage.
            self.close_install_logger()

    def ensure_symlink_runtime(self, sample_info: Mapping[str, Any]) -> str | None:
        return self._preparer.ensure_symlink_runtime(sample_info)

    def prepare_runtime(self, sample_info: Mapping[str, Any]) -> None:
        self.sample_info = dict(sample_info)
        # Fresh command index sequence for install/init/execution of this sample.
        self._command_counter = 0
        if self.logger is None:
            self.update_sample_logger(sample_info)
        try:
            # Install stage: use pack Installation_*.log only when we may actually install.
            # Reuse hits still log a short line there if the stamp file already exists.
            install_log = self._logger()
            if self.clone_shadow and self.PackID:
                try:
                    if not self._preparer.can_reuse_install():
                        install_log = self.ensure_install_logger()
                    else:
                        # Append a one-liner to the pack install log without "Start install".
                        from jarvishep2.sample_logger import SampleLogger

                        path = self.install_log_path()
                        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
                        pack = str(self.PackID)
                        reuse_log = SampleLogger.open(
                            path, module=f"{self.name}-{pack}"
                        )
                        install_log = reuse_log
                        self._install_logger = reuse_log
                        self._install_log_path = path
                except Exception:
                    install_log = self._logger()
            self._preparer.prepare(
                sample_info,
                run_stage=self._run_stage_commands,
                logger=install_log,
            )
        finally:
            # Install log (if opened) is pack-local; never leave it bound into init/exec.
            self.close_install_logger()

    def _logger(self):
        if self.logger is not None:
            return self.logger
        if isinstance(self.sample_info, dict):
            return self.sample_info.get("logger")
        return None

    def _stream_logger_for_stage(self, stage: str) -> Any | None:
        """Install → pack Installation_*.log; init/exec → Sample_running.log."""
        if str(stage or "").strip().lower() == "install":
            return self.ensure_install_logger()
        return self._logger()

    def update_sample_logger(self, sample_info: Mapping[str, Any]) -> None:
        """Bind the per-Sample logger as ``Sample@… (Module-No.pack)`` (V1 contract)."""
        parent = None
        if isinstance(sample_info, Mapping):
            parent = sample_info.get("logger")
        if parent is None:
            self.logger = None
            return
        base_name = ""
        if isinstance(sample_info, Mapping):
            base_name = str(sample_info.get("logger_name") or "").strip()
        if not base_name:
            uuid = ""
            if isinstance(sample_info, Mapping):
                uuid = str(sample_info.get("uuid") or "").strip()
            base_name = f"Sample@{uuid or 'UNKNOWN'}"
        pack = str(self.PackID or "NA")
        module_label = f"{base_name} ({self.name}-No.{pack})"
        binder = getattr(parent, "bind", None)
        if callable(binder):
            self.logger = binder(module=module_label)
        else:
            self.logger = parent
        if self.logger is not None:
            self.logger.info("Module load instance and logger is correctly set!")

    def close_sample_logger(self) -> None:
        """Detach the sample-local logger after the calculator step finishes."""
        self.logger = None

    def _resolve_runtime_tokens(self, text: str, *, stage: str, field: str) -> str:
        if text is None:
            return ""
        raw = str(text)
        if self._command_parser is None:
            if any(token in raw for token in ("@SampleID", "@Sdir", "@PackID")):
                raise RuntimeError(
                    f"CommandParser is required for sample token resolution "
                    f"(stage '{stage}', field '{field}')"
                )
            return raw
        return self._command_parser.resolve_sample(
            text,
            sample_info=self.sample_info,
            pack_id=self.PackID,
            stage=stage,
            field=field,
        )

    def _next_command_index(self) -> int:
        self._command_counter += 1
        return self._command_counter

    def _require_scheduler(self) -> AsyncSubprocessScheduler:
        if self._scheduler is None:
            raise RuntimeError("calculator subprocess scheduler is not attached")
        return self._scheduler

    def _run_command(
        self,
        command: Mapping[str, Any],
        *,
        stage: str = "execution",
        command_index: int = 0,
        timeout_sec: float | None = None,
    ) -> None:
        scheduler = self._require_scheduler()
        cmd_text = self._resolve_runtime_tokens(str(command.get("cmd", "")), stage=stage, field="cmd")
        cwd_text = self._resolve_runtime_tokens(str(command.get("cwd", ".")), stage=stage, field="cwd")
        cwd = os.path.abspath(cwd_text or ".")
        stream_logger = self._stream_logger_for_stage(stage)
        # Pure ``cd <path>`` only exists to update inherited cwd for later entries
        # (baked into each command's ``cwd`` at normalize time). Spawning a shell
        # ``cd`` is a no-op across process boundaries — log the switch and return.
        if is_cwd_switch_command(cmd_text):
            target = next_cwd_from_command(cmd_text, cwd)
            if not os.path.isabs(target):
                target = os.path.abspath(target)
            if stream_logger is not None:
                stream_logger.info(
                    f" Cwd switch [{stage}#{int(command_index):05}] -> \n"
                    f"\t{cmd_text} \n following commands cwd -> \n\t{target} \n"
                )
            return
        os.makedirs(cwd, exist_ok=True)
        # Install → pack Installation_*.log; init/exec → Sample_running.log (V1).
        use_stream_log = stream_logger is not None
        result = scheduler.run(
            SubprocessJob(
                cmd=cmd_text,
                cwd=cwd,
                env=self._subprocess_env,
                timeout_sec=timeout_sec,
                log_policy="logger" if use_stream_log else "quiet",
                stream_logger=stream_logger if use_stream_log else None,
                meta={
                    "module": self.name,
                    "pack_id": self.PackID,
                    "stage": stage,
                    "command_index": int(command_index),
                    "sample_uuid": (
                        self.sample_info.get("uuid")
                        if isinstance(self.sample_info, dict)
                        else None
                    ),
                    "cwd": cwd,
                    "command_log_to_stream": use_stream_log,
                    "emit_command_summary": use_stream_log,
                },
            ),
            timeout=(float(timeout_sec) + 5.0) if timeout_sec is not None else None,
        )
        if not result.ok:
            raise RuntimeError(
                self._format_command_failure(
                    stage=stage,
                    command_index=command_index,
                    cmd_text=cmd_text,
                    cwd=cwd,
                    result=result,
                )
            )

    def _format_command_failure(
        self,
        *,
        stage: str,
        command_index: int,
        cmd_text: str,
        cwd: str,
        result: Any,
    ) -> str:
        """Build a failure string that stands alone on the worker log.

        Sample_running.log already has the multi-line "Run … command" header;
        worker ERROR lines only carry ``str(exc)``, so cwd / module / stderr
        must be embedded here (relative cmds like ``./main`` are useless alone).
        """
        pack = self.PackID
        pack_s = "" if pack in (None, "", "None") else str(pack)
        lines = [
            f"Command failed [{stage}#{int(command_index):05}]",
            (
                f"  rc={int(getattr(result, 'returncode', -1))} "
                f"timeout={bool(getattr(result, 'timed_out', False))} "
                f"dur={float(getattr(result, 'duration_sec', 0.0)):.3f}s"
            ),
            f"  module={self.name}"
            + (f" pack_id={pack_s}" if pack_s else ""),
            f"  cwd={cwd}",
            f"  cmd={cmd_text}",
        ]
        stderr_tail = str(getattr(result, "stderr_tail", "") or "").strip()
        if stderr_tail:
            # Keep one logical line so worker two-column logs stay greppable.
            stderr_one = " ".join(stderr_tail.split())
            if len(stderr_one) > 400:
                stderr_one = stderr_one[:397] + "..."
            lines.append(f"  stderr={stderr_one}")
        stdout_tail = str(getattr(result, "stdout_tail", "") or "").strip()
        if stdout_tail and not stderr_tail:
            stdout_one = " ".join(stdout_tail.split())
            if len(stdout_one) > 400:
                stdout_one = stdout_one[:397] + "..."
            lines.append(f"  stdout={stdout_one}")
        return "\n".join(lines)

    def _run_stage_commands(self, commands: list[Mapping[str, Any]], stage: str) -> None:
        for command in commands:
            command_index = self._next_command_index()
            self._run_command(
                self._preparer.stage_command(command, stage=stage),
                stage=stage,
                command_index=command_index,
            )

    def _io_context(self, input_data: Mapping[str, Any] | None = None):
        return build_io_context(
            sample_info=self.sample_info if isinstance(self.sample_info, dict) else {},
            pack_id=self.PackID,
            module=self.name,
            runtime_values=input_data or {},
            logger=self._logger(),
            evaluate_expression=self._evaluate_io_expression,
        )

    def _portal_type_label(self, spec: Mapping[str, Any]) -> str:
        raw = str(spec.get("type", "")).strip() or "UNKNOWN"
        return f"Portal:{raw}"

    def _ensure_sample_dir_for_io(self) -> None:
        """Materialize SAMPLE/<uuid> when FileOperation needs a save/temp copy."""
        info = self.sample_info if isinstance(self.sample_info, dict) else {}
        if info.get("save_dir"):
            return
        needs = False
        for spec in (*self.input_specs, *self.output_specs):
            if bool(spec.get("save", False)):
                needs = True
                break
        # V1 always lands outputs under SAMPLE/.temp even when save is false.
        if self.output_specs:
            needs = True
        if not needs:
            return
        from jarvishep2.sample import ensure_sample_materialized

        ensure_sample_materialized(info)
        self.sample_info = info

    def load_input(self, input_data: Mapping[str, Any]) -> dict[str, Any]:
        """Write calculator inputs via Portal; SAMPLE save via HEP FileOperation.

        Merge policy (HEP-owned):
        1. Start with the sample input params/observables (authoritative for their names).
        2. Overlay Portal write observables only for **new** keys (expression dumps).
        3. Apply ``save: true`` copy through FileOperation (not Portal).
        """
        merged: dict[str, Any] = {key: input_data.get(key) for key in input_data}
        if not self.input_specs:
            return merged
        self._ensure_sample_dir_for_io()
        protected = set(merged)
        context = self._io_context(input_data)
        logger = self._logger()
        for spec in self.input_specs:
            if logger is not None:
                name = str(spec.get("name") or "").strip() or "input"
                logger.debug(f"Adding the file {name} as '{self._portal_type_label(spec)}' type")
            path = self._resolve_runtime_tokens(str(spec.get("path", "")), stage="execution", field="path")
            portal_obs = write_io_input_sync(spec, input_data, context=context, path=path)
            if isinstance(portal_obs, dict):
                for key, value in portal_obs.items():
                    if key in protected:
                        continue
                    merged[key] = value
            # HEP FileOperation: SAMPLE copy when save:true
            saved = apply_hep_io_save(
                source_path=path,
                sample_info=self.sample_info,
                module=self.name,
                spec=spec,
                direction="input",
                file_ops=self._file_ops,
            )
            if saved and "name" in spec and bool(spec.get("save", False)):
                name = str(spec["name"])
                if name not in protected:
                    merged[name] = saved
        return merged

    def read_output(self) -> dict[str, Any]:
        merged: dict[str, Any] = {}
        if not self.output_specs:
            return merged
        self._ensure_sample_dir_for_io()
        context = self._io_context()
        logger = self._logger()
        for spec in self.output_specs:
            if logger is not None:
                name = str(spec.get("name") or "").strip() or "output"
                logger.debug(f"Loading the file {name} as '{self._portal_type_label(spec)}' type")
            path = self._resolve_runtime_tokens(str(spec.get("path", "")), stage="execution", field="path")
            portal_obs = read_io_output_sync(spec, context=context, path=path)
            if isinstance(portal_obs, dict):
                merged.update(portal_obs)
            # HEP FileOperation: save:true → SAMPLE/; else → SAMPLE/.temp/ (V1)
            saved = apply_hep_io_save(
                source_path=path,
                sample_info=self.sample_info,
                module=self.name,
                spec=spec,
                direction="output",
                file_ops=self._file_ops,
            )
            if saved and "name" in spec:
                merged[str(spec["name"])] = saved
        return merged
    def execute_commands(self, *, deadline: float | None = None) -> None:
        for command in self.commands_template:
            command_index = self._next_command_index()
            timeout_sec = None
            if deadline is not None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise RuntimeError(
                        f"Calculator execution timed out after {self.timeout:g}s for {self.name}"
                    )
                timeout_sec = remaining
            staged = self._preparer.stage_command(command, stage="execution")
            self._run_command(
                staged,
                stage="execution",
                command_index=command_index,
                timeout_sec=timeout_sec,
            )

    def _selection_values(self, sample_info: Mapping[str, Any]) -> dict[str, Any]:
        values: dict[str, Any] = {}
        params = sample_info.get("params")
        if isinstance(params, Mapping):
            values.update({str(key): value for key, value in params.items()})
        observables = sample_info.get("observables")
        if isinstance(observables, Mapping):
            values.update({str(key): value for key, value in observables.items()})
        return values

    def _should_run(self, sample_info: Mapping[str, Any]) -> bool:
        """V1 module-level selection: skip this calculator when the cut is false."""
        if not self.selection:
            return True
        try:
            return bool(
                evaluate_selection(
                    self.selection,
                    self._selection_values(sample_info),
                    context=self._expression_context,
                )
            )
        except BoolConversionError as exc:
            raise BoolConversionError(
                f"Calculator selection for '{self.name}' failed: {exc}"
            ) from exc

    def should_run(self, sample_info: Mapping[str, Any]) -> bool:
        """Return whether this sample passes the module selection gate."""
        return self._should_run(sample_info)

    def execute(self, sample_info: Mapping[str, Any], *, runtime_prepared: bool = False) -> dict[str, Any]:
        """load_input → run commands → read_output (scheduler required when selected)."""
        self.sample_info = dict(sample_info)
        if not self.should_run(sample_info):
            return {}
        # Prefer a single bind for prepare+execute (Worker may prepare first).
        if self.logger is None:
            self.update_sample_logger(sample_info)
        self._require_scheduler()
        try:
            if not runtime_prepared:
                self.prepare_runtime(sample_info)
            input_data = dict(sample_info.get("observables") or sample_info.get("params") or {})
            result: dict[str, Any] = {}
            deadline = None
            if self.timeout is not None:
                deadline = time.monotonic() + float(self.timeout)
            result.update(self.load_input(input_data))
            self.execute_commands(deadline=deadline)
            result.update(self.read_output())
            return result
        except Exception as exc:
            logger = self._logger()
            if logger is not None:
                logger.error(f"Error during execution: {exc}")
            raise
        finally:
            self.close_sample_logger()
            self.sample_info = {}

    @classmethod
    def from_config_list(cls, modules: list[Mapping[str, Any]]) -> dict[str, CalculatorModule]:
        loaded: dict[str, CalculatorModule] = {}
        for item in modules:
            name = str(item.get("name", "")).strip()
            if not name:
                continue
            module = cls(name, item)
            module.preload_templates()
            loaded[name] = module
        return loaded


def mint_pack_id() -> str:
    """Deprecated local helper — distributed runs must use Redis ``acquire_calc``.

    Kept for unit tests that do not stand up a free-list. Prefer stable ids
    like ``001`` when exercising clone_shadow paths.
    """
    return "001"


__all__ = ["CalculatorModule", "CalculatorSpec", "RuntimePreparer", "mint_pack_id"]
