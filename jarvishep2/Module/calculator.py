#!/usr/bin/env python3
"""Calculator backend held long-term by Workers (V2 delta over V1)."""

from __future__ import annotations

import asyncio
import os
import signal
import time
from collections.abc import Mapping
from typing import Any
from uuid import uuid4

from jarvishep2.io_json import read_json_output, write_json_input
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
        self.sample_info: dict[str, Any] = {}
        self.PackID: str | None = None
        self._templates_loaded = False
        self._command_counter = 0

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
        self._templates_loaded = True

    def acquire_pack_id(self, pack_id: str) -> None:
        """Tag the current run for traceability."""
        self.PackID = str(pack_id)

    def _logger(self):
        if isinstance(self.sample_info, dict):
            return self.sample_info.get("logger")
        return None

    def _resolve_runtime_tokens(self, text: str, *, stage: str, field: str) -> str:
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

    async def run_command(
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

    async def execute_commands(self, *, deadline: float | None = None) -> None:
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
            await self.run_command(
                command,
                stage="execution",
                command_index=command_index,
                timeout_sec=timeout_sec,
            )

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

    def execute(self, sample_info: Mapping[str, Any]) -> dict[str, Any]:
        """Sync convenience entry: load_input → run commands → read_output."""
        self.sample_info = dict(sample_info)
        self._command_counter = 0
        input_data = dict(sample_info.get("observables") or sample_info.get("params") or {})
        try:
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