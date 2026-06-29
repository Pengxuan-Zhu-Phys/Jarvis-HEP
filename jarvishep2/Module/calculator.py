#!/usr/bin/env python3
"""Calculator backend held long-term by Workers (V2 delta over V1)."""

from __future__ import annotations

import asyncio
import os
import shutil
import signal
import time
from collections.abc import Mapping
from typing import Any
from uuid import uuid4

from jarvishep2.async_subprocess import AsyncSubprocessScheduler, SubprocessJob
from jarvishep2.command_parser import CommandParser
from jarvishep2.io_json import read_json_output, write_json_input
from jarvishep2.library import LibraryManager
from jarvishep2.sample import ensure_sample_materialized


class CalculatorModule:
    """External calculator runner with template preload and sync ``execute``."""

    def __init__(self, name: str, config: Mapping[str, Any]) -> None:
        self.name = str(name)
        self.config = dict(config)
        self.clone_shadow = bool(config.get("clone_shadow", False))
        self.installation = list(config.get("installation") or [])
        self.initialization = list(config.get("initialization") or [])
        execution = dict(config.get("execution") or {})
        self.execution = execution
        self.commands_template = list(execution.get("commands") or [])
        self.input_specs = list(execution.get("input") or [])
        self.output_specs = list(execution.get("output") or [])
        self.timeout = self._normalize_timeout(config.get("timeout"))
        self.basepath = str(config.get("path", execution.get("path", ".")))
        self.source = str(config.get("source", "")).strip()
        self.symlink_name = str(config.get("symlink_name", self.name)).strip() or self.name
        self.env_setup = list(config.get("env_setup") or [])
        self._subprocess_env: dict[str, str] | None = None
        self.sample_info: dict[str, Any] = {}
        self.PackID: str | None = None
        self._templates_loaded = False
        self._template_parse_count = 0
        self._command_counter = 0
        self._scheduler: AsyncSubprocessScheduler | None = None
        self._command_parser: CommandParser | None = None
        self._installed_shadows: set[str] = set()
        self._library = LibraryManager()

    @staticmethod
    def _normalize_timeout(value: Any) -> float | None:
        if value is None:
            return None
        timeout = float(value)
        return timeout if timeout > 0 else None

    @staticmethod
    def custom_format(record: Mapping[str, Any]) -> str:
        """V1-compatible log format hook (invariant #14 fix: accepts record)."""
        module = record.get("extra", {}).get("module", "No module")
        if "raw" in record.get("extra", {}):
            return "{message}\n"
        return f"{module} | {record.get('message', '')}"

    def preload_templates(self) -> None:
        """Parse and cache command/input templates once per Worker."""
        if self._templates_loaded:
            return
        # One-time structural read of cached template specs (no per-Sample work).
        _ = list(self.commands_template)
        _ = list(self.input_specs)
        _ = list(self.output_specs)
        self._template_parse_count += 1
        self._templates_loaded = True

    def attach_scheduler(self, scheduler: AsyncSubprocessScheduler | None) -> None:
        """Bind the per-Worker subprocess scheduler used for command execution."""
        self._scheduler = scheduler

    def attach_command_parser(self, parser: CommandParser | None) -> None:
        """Bind the per-Worker Phase-2 command resolver (WP-D3.1)."""
        self._command_parser = parser

    def bind_env(self, env: Mapping[str, str]) -> None:
        """Attach cached ``env_setup`` variables for subprocess execution."""
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
        """Tag the current run for traceability."""
        self.PackID = str(pack_id)

    def decode_shadow_path(self, path: str) -> str:
        """Resolve ``@PackID`` in calculator runtime paths."""
        if path is None:
            return ""
        resolved = str(path)
        if "@PackID" in resolved:
            if not self.PackID:
                raise RuntimeError("@PackID path requires pack_id before shadow decode")
            resolved = resolved.replace("@PackID", str(self.PackID))
        return os.path.abspath(resolved)

    def decode_shadow_commands(self, command: Mapping[str, Any]) -> dict[str, str]:
        """Resolve ``@PackID`` in install/init/execution command templates."""
        pack_id = str(self.PackID or "")
        if not pack_id:
            raise RuntimeError("@PackID command decode requires pack_id")
        return {
            "cmd": str(command.get("cmd", "")).replace("@PackID", pack_id),
            "cwd": str(command.get("cwd", ".")).replace("@PackID", pack_id),
        }

    def shadow_runtime_path(self) -> str:
        return self.decode_shadow_path(self.basepath)

    def _expand_install_tokens(self, text: str) -> str:
        runtime = self.shadow_runtime_path()
        source = os.path.abspath(self.source) if self.source else ""
        return str(text).replace("${source}", source).replace("${path}", runtime)

    def _stage_command(self, command: Mapping[str, Any], *, stage: str) -> dict[str, str]:
        raw = dict(command)
        if self.clone_shadow:
            raw = self.decode_shadow_commands(raw)
        if stage in {"install", "initialize"}:
            raw["cmd"] = self._expand_install_tokens(str(raw.get("cmd", "")))
            raw["cwd"] = self._expand_install_tokens(str(raw.get("cwd", ".")))
        return {"cmd": str(raw.get("cmd", "")), "cwd": str(raw.get("cwd", "."))}

    def _run_stage_commands_sync(self, commands: list[Mapping[str, Any]], *, stage: str) -> None:
        for command in commands:
            command_index = self._next_command_index()
            self._run_command_sync(
                self._stage_command(command, stage=stage),
                stage=stage,
                command_index=command_index,
            )

    def ensure_shadow_installed(self) -> None:
        """Physical per-pack install for ``clone_shadow`` calculators."""
        if not self.clone_shadow:
            return
        pack_id = str(self.PackID or "")
        if not pack_id:
            raise RuntimeError("clone_shadow install requires pack_id")
        if pack_id in self._installed_shadows:
            return
        runtime = self.shadow_runtime_path()
        os.makedirs(runtime, exist_ok=True)
        if self.installation:
            self._run_stage_commands_sync(self.installation, stage="install")
        elif self.source:
            shutil.copytree(os.path.abspath(self.source), runtime, dirs_exist_ok=True)
        else:
            raise RuntimeError(
                f"clone_shadow calculator '{self.name}' requires a source path or installation commands"
            )
        if self.initialization:
            self._run_stage_commands_sync(self.initialization, stage="initialize")
        self._installed_shadows.add(pack_id)

    def ensure_symlink_runtime(self, sample_info: Mapping[str, Any]) -> str | None:
        """Symlink a safe tool into the Sample dir when ``clone_shadow`` is false."""
        if self.clone_shadow or not self.source:
            return None
        save_dir = ensure_sample_materialized(dict(sample_info))
        if save_dir is None:
            raise RuntimeError(f"symlink runtime requires materialized save_dir for '{self.name}'")
        return self._library.link_into_sample(self.source, str(save_dir), self.symlink_name)

    def prepare_runtime(self, sample_info: Mapping[str, Any]) -> None:
        """Prepare per-run isolation before ``execute``."""
        self.sample_info = dict(sample_info)
        if self.clone_shadow:
            self.ensure_shadow_installed()
        else:
            self.ensure_symlink_runtime(sample_info)

    def _logger(self):
        if isinstance(self.sample_info, dict):
            return self.sample_info.get("logger")
        return None

    def _resolve_runtime_tokens(self, text: str, *, stage: str, field: str) -> str:
        if self._command_parser is not None:
            return self._command_parser.resolve_sample(
                text,
                sample_info=self.sample_info,
                pack_id=self.PackID,
                stage=stage,
                field=field,
            )
        if text is None:
            return ""
        raw = str(text)
        if stage == "install" or ("@SampleID" not in raw and "@Sdir" not in raw and "@PackID" not in raw):
            return raw
        if not isinstance(self.sample_info, dict):
            raise RuntimeError(
                f"Runtime token requires sample_info during stage '{stage}' for field '{field}'"
            )
        resolved = raw
        if "@PackID" in resolved:
            if not self.PackID:
                raise RuntimeError(f"@PackID requires pack_id during stage '{stage}' for field '{field}'")
            resolved = resolved.replace("@PackID", str(self.PackID))
        if "@SampleID" in resolved:
            sample_uuid = self.sample_info.get("uuid")
            if not sample_uuid:
                raise RuntimeError(f"@SampleID requires uuid during stage '{stage}' for field '{field}'")
            resolved = resolved.replace("@SampleID", str(sample_uuid))
        if "@Sdir" in resolved:
            save_dir = ensure_sample_materialized(self.sample_info)
            if save_dir is None:
                raise RuntimeError(f"@Sdir requires save_dir during stage '{stage}' for field '{field}'")
            resolved = resolved.replace("@Sdir", str(save_dir))
        return resolved

    def _next_command_index(self) -> int:
        self._command_counter += 1
        return self._command_counter

    @staticmethod
    async def _terminate_process(process: asyncio.subprocess.Process, grace_sec: float = 5.0) -> None:
        pid = int(getattr(process, "pid", 0) or 0)
        try:
            if pid > 0:
                os.killpg(pid, signal.SIGTERM)
            else:
                process.terminate()
        except Exception:
            process.terminate()
        try:
            await asyncio.wait_for(process.wait(), timeout=max(0.1, grace_sec))
        except Exception:
            try:
                if pid > 0:
                    os.killpg(pid, signal.SIGKILL)
                else:
                    process.kill()
            except Exception:
                process.kill()

    def _run_command_sync(
        self,
        command: Mapping[str, Any],
        *,
        stage: str = "execution",
        command_index: int = 0,
        timeout_sec: float | None = None,
    ) -> None:
        cmd_text = self._resolve_runtime_tokens(str(command.get("cmd", "")), stage=stage, field="cmd")
        cwd_text = self._resolve_runtime_tokens(str(command.get("cwd", ".")), stage=stage, field="cwd")
        cwd = os.path.abspath(cwd_text or ".")
        os.makedirs(cwd, exist_ok=True)
        if self._scheduler is None:
            raise RuntimeError("calculator subprocess scheduler is not attached")

        meta = {
            "module": self.name,
            "pack_id": self.PackID,
            "stage": stage,
            "command_index": command_index,
        }
        result = self._scheduler.run(
            SubprocessJob(
                cmd=cmd_text,
                cwd=cwd,
                env=self._subprocess_env,
                timeout_sec=timeout_sec,
                log_policy="quiet",
                meta=meta,
            ),
            timeout=(float(timeout_sec) + 5.0) if timeout_sec is not None else None,
        )
        if not result.ok:
            raise RuntimeError(
                f"Command failed [{stage}#{command_index:05}] rc={result.returncode} "
                f"timeout={result.timed_out} cmd={cmd_text}"
            )

    async def run_command(
        self,
        command: Mapping[str, Any],
        *,
        stage: str = "execution",
        command_index: int = 0,
        timeout_sec: float | None = None,
    ) -> None:
        if self._scheduler is not None:
            self._run_command_sync(
                command,
                stage=stage,
                command_index=command_index,
                timeout_sec=timeout_sec,
            )
            return

        cmd_text = self._resolve_runtime_tokens(str(command.get("cmd", "")), stage=stage, field="cmd")
        cwd_text = self._resolve_runtime_tokens(str(command.get("cwd", ".")), stage=stage, field="cwd")
        cwd = os.path.abspath(cwd_text or ".")
        os.makedirs(cwd, exist_ok=True)

        process = await asyncio.create_subprocess_shell(
            cmd_text,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=cwd,
            start_new_session=True,
        )
        timed_out = False
        try:
            if timeout_sec is None:
                rc = await process.wait()
            else:
                rc = await asyncio.wait_for(process.wait(), timeout=float(timeout_sec))
        except asyncio.TimeoutError:
            timed_out = True
            await self._terminate_process(process)
            rc = await process.wait()
        if int(rc) != 0 or timed_out:
            raise RuntimeError(
                f"Command failed [{stage}#{command_index:05}] rc={rc} timeout={timed_out} cmd={cmd_text}"
            )

    async def load_input(self, input_data: Mapping[str, Any]) -> dict[str, Any]:
        merged: dict[str, Any] = {}
        for spec in self.input_specs:
            if str(spec.get("type", "")).strip() != "JSON":
                continue
            path = self._resolve_runtime_tokens(str(spec.get("path", "")), stage="execution", field="path")
            actions = list(spec.get("actions") or [])
            write_json_input(path, actions=actions, param_values=input_data)
            merged.update({key: input_data.get(key) for key in input_data})
        return merged

    async def read_output(self) -> dict[str, Any]:
        merged: dict[str, Any] = {}
        for spec in self.output_specs:
            if str(spec.get("type", "")).strip() != "JSON":
                continue
            path = self._resolve_runtime_tokens(str(spec.get("path", "")), stage="execution", field="path")
            variables = list(spec.get("variables") or [])
            merged.update(read_json_output(path, variables=variables))
        return merged

    def _execute_commands_sync(self, *, deadline: float | None = None) -> None:
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
            staged = self._stage_command(command, stage="execution")
            self._run_command_sync(
                staged,
                stage="execution",
                command_index=command_index,
                timeout_sec=timeout_sec,
            )

    async def execute_commands(self, *, deadline: float | None = None) -> None:
        if self._scheduler is not None:
            self._execute_commands_sync(deadline=deadline)
            return
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
            staged = self._stage_command(command, stage="execution")
            await self.run_command(
                staged,
                stage="execution",
                command_index=command_index,
                timeout_sec=timeout_sec,
            )

    def _load_input_sync(self, input_data: Mapping[str, Any]) -> dict[str, Any]:
        merged: dict[str, Any] = {}
        for spec in self.input_specs:
            if str(spec.get("type", "")).strip() != "JSON":
                continue
            path = self._resolve_runtime_tokens(str(spec.get("path", "")), stage="execution", field="path")
            actions = list(spec.get("actions") or [])
            write_json_input(path, actions=actions, param_values=input_data)
            merged.update({key: input_data.get(key) for key in input_data})
        return merged

    def _read_output_sync(self) -> dict[str, Any]:
        merged: dict[str, Any] = {}
        for spec in self.output_specs:
            if str(spec.get("type", "")).strip() != "JSON":
                continue
            path = self._resolve_runtime_tokens(str(spec.get("path", "")), stage="execution", field="path")
            variables = list(spec.get("variables") or [])
            merged.update(read_json_output(path, variables=variables))
        return merged

    def _execute_sync(self, input_data: Mapping[str, Any]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        deadline = None
        if self.timeout is not None:
            deadline = time.monotonic() + float(self.timeout)

        input_obs = self._load_input_sync(input_data)
        if isinstance(input_obs, dict):
            result.update(input_obs)

        self._execute_commands_sync(deadline=deadline)

        output_obs = self._read_output_sync()
        if isinstance(output_obs, dict):
            result.update(output_obs)
        return result

    async def _execute_async(self, input_data: Mapping[str, Any]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        deadline = None
        if self.timeout is not None:
            deadline = time.monotonic() + float(self.timeout)

        input_obs = await self.load_input(input_data)
        if isinstance(input_obs, dict):
            result.update(input_obs)

        await self.execute_commands(deadline=deadline)

        output_obs = await self.read_output()
        if isinstance(output_obs, dict):
            result.update(output_obs)
        return result

    def execute(self, sample_info: Mapping[str, Any], *, runtime_prepared: bool = False) -> dict[str, Any]:
        """Sync convenience entry: load_input → run commands → read_output."""
        self.sample_info = dict(sample_info)
        self._command_counter = 0
        if not runtime_prepared:
            self.prepare_runtime(sample_info)
        input_data = dict(sample_info.get("observables") or sample_info.get("params") or {})
        try:
            if self._scheduler is not None:
                return self._execute_sync(input_data)
            return asyncio.run(self._execute_async(input_data))
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


__all__ = ["CalculatorModule", "mint_pack_id"]