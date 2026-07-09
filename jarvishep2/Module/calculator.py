#!/usr/bin/env python3
"""Calculator backend held long-term by Workers (sync scheduler path only)."""

from __future__ import annotations

import os
import time
from collections.abc import Mapping
from typing import Any
from uuid import uuid4

from jarvishep2.async_subprocess import AsyncSubprocessScheduler, SubprocessJob
from jarvishep2.command_parser import CommandParser
from jarvishep2.io_portal import (
    build_io_context,
    read_io_output_sync,
    write_io_input_sync,
)
from jarvishep2.Module.calculator_spec import CalculatorSpec
from jarvishep2.Module.runtime_preparer import RuntimePreparer


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
        self._subprocess_env: dict[str, str] | None = None
        self.sample_info: dict[str, Any] = {}
        self.PackID: str | None = None
        self._templates_loaded = False
        self._template_parse_count = 0
        self._command_counter = 0
        self._scheduler: AsyncSubprocessScheduler | None = None
        self._command_parser: CommandParser | None = None
        self._preparer = RuntimePreparer(self.spec)

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
        self.PackID = str(pack_id)
        self._preparer.acquire_pack_id(pack_id)

    def decode_shadow_path(self, path: str) -> str:
        return self._preparer.decode_shadow_path(path)

    def decode_shadow_commands(self, command: Mapping[str, Any]) -> dict[str, str]:
        return self._preparer.decode_shadow_commands(command)

    def shadow_runtime_path(self) -> str:
        return self._preparer.shadow_runtime_path()

    def ensure_shadow_installed(self) -> None:
        self._preparer.ensure_shadow_installed(run_stage=self._run_stage_commands)

    def ensure_symlink_runtime(self, sample_info: Mapping[str, Any]) -> str | None:
        return self._preparer.ensure_symlink_runtime(sample_info)

    def prepare_runtime(self, sample_info: Mapping[str, Any]) -> None:
        self.sample_info = dict(sample_info)
        self._preparer.prepare(sample_info, run_stage=self._run_stage_commands)

    def _logger(self):
        if isinstance(self.sample_info, dict):
            return self.sample_info.get("logger")
        return None

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
        os.makedirs(cwd, exist_ok=True)
        result = scheduler.run(
            SubprocessJob(
                cmd=cmd_text,
                cwd=cwd,
                env=self._subprocess_env,
                timeout_sec=timeout_sec,
                log_policy="quiet",
                meta={
                    "module": self.name,
                    "pack_id": self.PackID,
                    "stage": stage,
                    "command_index": command_index,
                },
            ),
            timeout=(float(timeout_sec) + 5.0) if timeout_sec is not None else None,
        )
        if not result.ok:
            raise RuntimeError(
                f"Command failed [{stage}#{command_index:05}] rc={result.returncode} "
                f"timeout={result.timed_out} cmd={cmd_text}"
            )

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
        )

    def load_input(self, input_data: Mapping[str, Any]) -> dict[str, Any]:
        merged: dict[str, Any] = {key: input_data.get(key) for key in input_data}
        if not self.input_specs:
            return merged
        context = self._io_context(input_data)
        for spec in self.input_specs:
            path = self._resolve_runtime_tokens(str(spec.get("path", "")), stage="execution", field="path")
            portal_obs = write_io_input_sync(spec, input_data, context=context, path=path)
            if isinstance(portal_obs, dict):
                merged.update(portal_obs)
        return merged

    def read_output(self) -> dict[str, Any]:
        merged: dict[str, Any] = {}
        if not self.output_specs:
            return merged
        context = self._io_context()
        for spec in self.output_specs:
            path = self._resolve_runtime_tokens(str(spec.get("path", "")), stage="execution", field="path")
            portal_obs = read_io_output_sync(spec, context=context, path=path)
            if isinstance(portal_obs, dict):
                merged.update(portal_obs)
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

    def execute(self, sample_info: Mapping[str, Any], *, runtime_prepared: bool = False) -> dict[str, Any]:
        """load_input → run commands → read_output (scheduler required)."""
        self._require_scheduler()
        self.sample_info = dict(sample_info)
        self._command_counter = 0
        if not runtime_prepared:
            self.prepare_runtime(sample_info)
        input_data = dict(sample_info.get("observables") or sample_info.get("params") or {})
        try:
            result: dict[str, Any] = {}
            deadline = None
            if self.timeout is not None:
                deadline = time.monotonic() + float(self.timeout)
            result.update(self.load_input(input_data))
            self.execute_commands(deadline=deadline)
            result.update(self.read_output())
            return result
        finally:
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
    """Generate a fresh pack_id for a calculator run (D1.2 local; D2 uses Redis pool)."""
    return str(uuid4())


__all__ = ["CalculatorModule", "CalculatorSpec", "RuntimePreparer", "mint_pack_id"]
