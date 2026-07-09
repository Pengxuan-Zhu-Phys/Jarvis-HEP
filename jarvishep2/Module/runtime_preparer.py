#!/usr/bin/env python3
"""Per-pack shadow install / symlink preparation (D9.3 RuntimePreparer)."""

from __future__ import annotations

import os
import shutil
from collections.abc import Callable, Mapping
from typing import Any

from jarvishep2.library import LibraryManager
from jarvishep2.Module.calculator_spec import CalculatorSpec
from jarvishep2.sample import ensure_sample_materialized


class RuntimePreparer:
    """Owns clone_shadow physical install and safe-tool symlink isolation."""

    def __init__(
        self,
        spec: CalculatorSpec,
        *,
        library: LibraryManager | None = None,
    ) -> None:
        self.spec = spec
        self._library = library or LibraryManager()
        self._installed_shadows: set[str] = set()
        self.pack_id: str | None = None

    def acquire_pack_id(self, pack_id: str) -> None:
        self.pack_id = str(pack_id)

    def decode_shadow_path(self, path: str) -> str:
        if path is None:
            return ""
        resolved = str(path)
        if "@PackID" in resolved:
            if not self.pack_id:
                raise RuntimeError("@PackID path requires pack_id before shadow decode")
            resolved = resolved.replace("@PackID", str(self.pack_id))
        return os.path.abspath(resolved)

    def decode_shadow_commands(self, command: Mapping[str, Any]) -> dict[str, str]:
        pack_id = str(self.pack_id or "")
        if not pack_id:
            raise RuntimeError("@PackID command decode requires pack_id")
        return {
            "cmd": str(command.get("cmd", "")).replace("@PackID", pack_id),
            "cwd": str(command.get("cwd", ".")).replace("@PackID", pack_id),
        }

    def shadow_runtime_path(self) -> str:
        return self.decode_shadow_path(self.spec.basepath)

    def _expand_install_tokens(self, text: str) -> str:
        runtime = self.shadow_runtime_path()
        source = os.path.abspath(self.spec.source) if self.spec.source else ""
        return str(text).replace("${source}", source).replace("${path}", runtime)

    def stage_command(self, command: Mapping[str, Any], *, stage: str) -> dict[str, str]:
        raw = dict(command)
        if self.spec.clone_shadow:
            raw = self.decode_shadow_commands(raw)
        if stage in {"install", "initialize"}:
            raw["cmd"] = self._expand_install_tokens(str(raw.get("cmd", "")))
            raw["cwd"] = self._expand_install_tokens(str(raw.get("cwd", ".")))
        return {"cmd": str(raw.get("cmd", "")), "cwd": str(raw.get("cwd", "."))}

    def ensure_shadow_installed(
        self,
        *,
        run_stage: Callable[[list[Mapping[str, Any]], str], None],
    ) -> None:
        if not self.spec.clone_shadow:
            return
        pack_id = str(self.pack_id or "")
        if not pack_id:
            raise RuntimeError("clone_shadow install requires pack_id")
        if pack_id in self._installed_shadows:
            return
        runtime = self.shadow_runtime_path()
        os.makedirs(runtime, exist_ok=True)
        if self.spec.installation:
            run_stage(list(self.spec.installation), "install")
        elif self.spec.source:
            shutil.copytree(os.path.abspath(self.spec.source), runtime, dirs_exist_ok=True)
        else:
            raise RuntimeError(
                f"clone_shadow calculator '{self.spec.name}' requires a source path or installation commands"
            )
        if self.spec.initialization:
            run_stage(list(self.spec.initialization), "initialize")
        self._installed_shadows.add(pack_id)

    def ensure_symlink_runtime(self, sample_info: Mapping[str, Any]) -> str | None:
        if self.spec.clone_shadow or not self.spec.source:
            return None
        save_dir = ensure_sample_materialized(dict(sample_info))
        if save_dir is None:
            raise RuntimeError(
                f"symlink runtime requires materialized save_dir for '{self.spec.name}'"
            )
        return self._library.link_into_sample(
            self.spec.source, str(save_dir), self.spec.symlink_name
        )

    def prepare(
        self,
        sample_info: Mapping[str, Any],
        *,
        run_stage: Callable[[list[Mapping[str, Any]], str], None],
    ) -> None:
        if self.spec.clone_shadow:
            self.ensure_shadow_installed(run_stage=run_stage)
        else:
            self.ensure_symlink_runtime(sample_info)


__all__ = ["RuntimePreparer"]
