#!/usr/bin/env python3
"""Per-pack shadow install / symlink preparation (D9.3 RuntimePreparer)."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from typing import Any

from jarvishep2.library import LibraryManager
from jarvishep2.calculators.calculator_spec import CalculatorSpec
from jarvishep2.sample import ensure_sample_materialized

# Written under the pack shadow directory after a successful install.
INSTALL_STAMP_BASENAME = ".jarvis_install_stamp.json"
INSTALL_STAMP_SCHEMA = "jarvishep2.calc_install/v1"
INSTALL_CONTROL_BASENAME = "jarvis_install.json"
INSTALL_CONTROL_SCHEMA = "jarvishep2.calc_install/v2"


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _source_stat_payload(source: str) -> dict[str, Any] | None:
    """Cheap source identity: path + root stat (not a full tree hash)."""
    text = str(source or "").strip()
    if not text:
        return None
    path = os.path.abspath(text)
    if not os.path.exists(path):
        return {"path": path, "exists": False}
    try:
        st = os.stat(path)
        return {
            "path": path,
            "exists": True,
            "mtime_ns": int(getattr(st, "st_mtime_ns", int(st.st_mtime * 1e9))),
            "size": int(st.st_size),
            "is_dir": os.path.isdir(path),
        }
    except OSError:
        return {"path": path, "exists": False}


def build_install_fingerprint(
    spec: CalculatorSpec,
    *,
    pack_id: str | None = None,
) -> str:
    """Stable hash of calculator *install settings* (not per-sample init).

    Includes: module name, basepath template, source root identity, and the
    installation command list. PackID is intentionally **not** mixed into the
    hash so 001/002 share the same config identity; each pack still has its own
    stamp file under its shadow directory.
    """
    del pack_id  # reserved for future per-pack overrides
    commands: list[dict[str, str]] = []
    for item in spec.installation:
        if not isinstance(item, Mapping):
            continue
        commands.append(
            {
                "cmd": str(item.get("cmd", "")),
                "cwd": str(item.get("cwd", ".")),
            }
        )
    payload = {
        "schema": INSTALL_STAMP_SCHEMA,
        "module": str(spec.name),
        "clone_shadow": bool(spec.clone_shadow),
        "basepath": str(spec.basepath),
        "source": str(spec.source or ""),
        "source_stat": _source_stat_payload(spec.source),
        "installation": commands,
        "mode_parent": str(spec.mode_parent or ""),
        "mode_name": str(spec.mode_name or ""),
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def build_parent_install_fingerprint(spec: CalculatorSpec) -> str:
    """Hash the physical shared base, excluding mode-specific rebuild commands."""
    commands: list[dict[str, str]] = []
    for item in spec.parent_installation:
        if not isinstance(item, Mapping):
            continue
        commands.append({
            "cmd": str(item.get("cmd", "")),
            "cwd": str(item.get("cwd", ".")),
        })
    payload = {
        "schema": "jarvishep2.calc_shared_base/v1",
        "module": str(spec.mode_parent or spec.name),
        "clone_shadow": bool(spec.clone_shadow),
        "basepath": str(spec.basepath),
        "source": str(spec.source or ""),
        "source_stat": _source_stat_payload(spec.source),
        "parent_installation": commands,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def install_stamp_path(runtime_dir: str) -> str:
    return os.path.join(os.path.abspath(str(runtime_dir)), INSTALL_STAMP_BASENAME)


def install_control_path(spec: CalculatorSpec) -> str | None:
    """Return the control-file path beside the ``@PackID`` directories."""
    if not spec.clone_shadow or "@PackID" not in str(spec.basepath):
        return None
    pack_path = os.path.abspath(str(spec.basepath).replace("@PackID", "__jarvis_pack__"))
    return os.path.join(os.path.dirname(pack_path), INSTALL_CONTROL_BASENAME)


def read_install_stamp(runtime_dir: str) -> dict[str, Any] | None:
    path = install_stamp_path(runtime_dir)
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def write_install_stamp(
    runtime_dir: str,
    *,
    fingerprint: str,
    module: str,
    pack_id: str,
    epoch: int = 0,
    extra: Mapping[str, Any] | None = None,
) -> str:
    path = install_stamp_path(runtime_dir)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    payload: dict[str, Any] = {
        "schema": INSTALL_STAMP_SCHEMA,
        "module": str(module),
        "pack_id": str(pack_id),
        "epoch": max(0, int(epoch)),
        "fingerprint": str(fingerprint),
        "installed_at_utc": _utc_now_iso(),
    }
    if extra:
        payload["extra"] = dict(extra)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(tmp, path)
    return path


def read_install_control(spec: CalculatorSpec) -> dict[str, Any] | None:
    path = install_control_path(spec)
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def write_install_control(spec: CalculatorSpec, control: Mapping[str, Any]) -> str | None:
    """Atomically write the control-process-owned calculator state file."""
    path = install_control_path(spec)
    if not path:
        return None
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(dict(control), handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(tmp, path)
    return path


def force_calc_install_requested(env: Mapping[str, str] | None = None) -> bool:
    """True when operators force a full reinstall (env or future CLI)."""
    environ = env if env is not None else os.environ
    raw = str(environ.get("JARVIS_FORCE_CALC_INSTALL", "") or "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


def prepare_install_controls(
    modules: Sequence[Mapping[str, Any]] | None,
    *,
    logger: Any | None = None,
    env: Mapping[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Bump reinstall epochs once in the control process and plumb them to Workers."""
    force = force_calc_install_requested(env)
    prepared: list[dict[str, Any]] = []
    # A multi-mode parent has several logical CalculatorSpecs but exactly one
    # physical ``jarvis_install.json`` beside its shared @PackID directories.
    # Remember the epoch after the first child so ``reinstall: true`` advances
    # it once, not once per mode.
    shared_controls: dict[str, tuple[int, bool]] = {}
    for raw in modules or []:
        config = dict(raw)
        spec = CalculatorSpec.from_config(str(config.get("name", "Calculator")), config)
        if not spec.clone_shadow:
            prepared.append(config)
            continue
        control_path = install_control_path(spec)
        if spec.shared_mode_pack and control_path in shared_controls:
            next_epoch, requested = shared_controls[control_path]
            config["_install_epoch"] = next_epoch
            prepared.append(config)
            continue
        state = dict(read_install_control(spec) or {})
        try:
            epoch = max(0, int(state.get("reinstall_epoch", 0) or 0))
        except (TypeError, ValueError):
            epoch = 0
        requested = bool(state.get("reinstall")) or force
        next_epoch = epoch + 1 if requested else epoch
        control = {
            "schema": INSTALL_CONTROL_SCHEMA,
            "module": spec.mode_parent if spec.shared_mode_pack else spec.name,
            "source": os.path.abspath(spec.source) if spec.source else "",
            "reinstall": False,
            "reinstall_epoch": next_epoch,
            "packs": dict(state.get("packs") or {}),
            "_help": "Edited your calculator source? Set \"reinstall\": true and run again.",
        }
        write_install_control(spec, control)
        config["_install_epoch"] = next_epoch
        prepared.append(config)
        if spec.shared_mode_pack and control_path:
            shared_controls[control_path] = (next_epoch, requested)
        if logger is not None:
            try:
                logger.info(
                    "%s: %s calculator install epoch %d%s",
                    spec.name,
                    "reinstalling" if requested else "reusing",
                    next_epoch,
                    " (operator requested)" if requested else "",
                )
            except Exception:
                pass
    return prepared


def refresh_install_control_summaries(
    modules: Sequence[Mapping[str, Any]] | None,
) -> None:
    """Refresh best-effort per-pack summaries (control process only)."""
    seen_controls: set[str] = set()
    for raw in modules or []:
        config = dict(raw)
        spec = CalculatorSpec.from_config(str(config.get("name", "Calculator")), config)
        control_path = install_control_path(spec)
        if control_path and control_path in seen_controls:
            continue
        if control_path:
            seen_controls.add(control_path)
        control = read_install_control(spec)
        if not control or not spec.clone_shadow or "@PackID" not in str(spec.basepath):
            continue
        template = os.path.abspath(str(spec.basepath).replace("@PackID", "__jarvis_pack__"))
        parent = os.path.dirname(template)
        packs: dict[str, Any] = {}
        try:
            entries = os.listdir(parent)
        except OSError:
            entries = []
        for entry in entries:
            if not entry.isdigit():
                continue
            stamp = read_install_stamp(os.path.join(parent, entry))
            if stamp:
                extra = stamp.get("extra")
                packs[str(entry)] = {
                    "fingerprint": stamp.get("fingerprint", ""),
                    "epoch": int(stamp.get("epoch", 0) or 0),
                    "installed_at_utc": stamp.get("installed_at_utc", ""),
                    "mode": str(extra.get("mode") or "") if isinstance(extra, Mapping) else "",
                }
        control["packs"] = packs
        write_install_control(spec, control)


class RuntimePreparer:
    """Owns clone_shadow physical install and safe-tool symlink isolation."""

    def __init__(
        self,
        spec: CalculatorSpec,
        *,
        library: LibraryManager | None = None,
        force_reinstall: bool = False,
        install_epoch: int = 0,
    ) -> None:
        self.spec = spec
        self._library = library or LibraryManager()
        self._installed_shadows: set[str] = set()
        self.pack_id: str | None = None
        # Module/task override; also ORed with JARVIS_FORCE_CALC_INSTALL.
        self.force_reinstall = bool(force_reinstall)
        self.install_epoch = max(0, int(install_epoch))
        self.last_install_action: str | None = None  # "installed" | "reused" | None

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

    def install_fingerprint(self) -> str:
        return build_install_fingerprint(self.spec, pack_id=self.pack_id)

    def parent_install_fingerprint(self) -> str:
        return build_parent_install_fingerprint(self.spec)

    def _shared_parent_is_ready(
        self,
        runtime: str,
        *,
        stamp: Mapping[str, Any] | None = None,
    ) -> bool:
        """Whether a shared pack's common base can be retained for a mode switch."""
        if self.force_reinstall:
            return False
        current = stamp if stamp is not None else read_install_stamp(runtime)
        if not isinstance(current, Mapping):
            return False
        extra = current.get("extra")
        if not isinstance(extra, Mapping):
            return False
        try:
            parent_epoch = max(0, int(extra.get("parent_install_epoch", 0) or 0))
        except (TypeError, ValueError):
            parent_epoch = 0
        return (
            str(extra.get("parent_install_fingerprint") or "")
            == self.parent_install_fingerprint()
            and parent_epoch >= self.install_epoch
        )

    def can_reuse_install(self, runtime: str | None = None) -> bool:
        """True when pack dir has a matching install stamp and force is off."""
        if self.force_reinstall:
            return False
        root = runtime or self.shadow_runtime_path()
        if not os.path.isdir(root):
            return False
        stamp = read_install_stamp(root)
        if not stamp:
            return False
        expected = self.install_fingerprint()
        try:
            stamp_epoch = max(0, int(stamp.get("epoch", 0) or 0))
        except (TypeError, ValueError):
            stamp_epoch = 0
        stamp_extra = stamp.get("extra")
        stamp_mode = (
            str(stamp_extra.get("mode") or "")
            if isinstance(stamp_extra, Mapping)
            else ""
        )
        return (
            str(stamp.get("fingerprint") or "") == expected
            and stamp_epoch >= self.install_epoch
            and (
                not self.spec.shared_mode_pack
                or stamp_mode == str(self.spec.mode_name or "")
            )
        )

    def ensure_shadow_installed(
        self,
        *,
        run_stage: Callable[[list[Mapping[str, Any]], str], None],
        logger: Any | None = None,
    ) -> None:
        """Physical install of a shadow pack (once per matching config).

        Skip policy
        -----------
        If ``{runtime}/.jarvis_install_stamp.json`` exists and its fingerprint
        matches the current calculator install settings, reuse the pack without
        re-running ``installation`` commands. Force with
        ``force_reinstall=True``. The control-process environment flag is
        converted into a monotone epoch before Workers spawn.

        Process-local ``_installed_shadows`` still short-circuits within one
        Worker after the first successful install/reuse of a PackID.
        """
        if not self.spec.clone_shadow:
            return
        pack_id = str(self.pack_id or "")
        if not pack_id:
            raise RuntimeError("clone_shadow install requires pack_id")
        # Ordinary modules own their pool, so a process-local reuse shortcut is
        # safe. Shared multi-mode packs can have been rebuilt by another Worker
        # since this module last touched the numeric PackID; always re-read its
        # on-disk mode stamp in that case.
        if pack_id in self._installed_shadows and not self.spec.shared_mode_pack:
            self.last_install_action = "reused"
            return

        runtime = self.shadow_runtime_path()
        os.makedirs(runtime, exist_ok=True)
        fingerprint = self.install_fingerprint()

        if self.can_reuse_install(runtime):
            if not self.spec.shared_mode_pack:
                self._installed_shadows.add(pack_id)
            self.last_install_action = "reused"
            if logger is not None:
                try:
                    logger.info(
                        f"Skip calculator install (settings unchanged) -> "
                        f"module={self.spec.name} pack_id={pack_id}\n"
                        f"\t runtime -> {runtime}\n"
                        f"\t fingerprint -> {fingerprint[:16]}…"
                    )
                except Exception:
                    pass
            return

        old_stamp = read_install_stamp(runtime)
        # Shared packs may currently hold a different mode. Clear the old
        # stamp before touching the directory: a failed rebuild must never be
        # mistaken for a healthy copy of the prior mode on the next acquire.
        try:
            os.unlink(install_stamp_path(runtime))
        except FileNotFoundError:
            pass
        if self.spec.shared_mode_pack:
            # The parent block builds the common physical base. A normal mode
            # switch retains that base and runs only this mode's installation.
            parent_needed = not self._shared_parent_is_ready(runtime, stamp=old_stamp)
            if parent_needed:
                if self.spec.parent_installation:
                    run_stage(list(self.spec.parent_installation), "install")
                elif self.spec.source:
                    shutil.copytree(
                        os.path.abspath(self.spec.source), runtime, dirs_exist_ok=True
                    )
                elif not self.spec.mode_installation:
                    raise RuntimeError(
                        f"shared calculator '{self.spec.mode_parent or self.spec.name}' "
                        "requires a source path, parent installation, or mode installation"
                    )
            if self.spec.mode_installation:
                run_stage(list(self.spec.mode_installation), "install")
        # Normal calculators retain the original one-stage contract.
        elif self.spec.installation:
            run_stage(list(self.spec.installation), "install")
        elif self.spec.source:
            shutil.copytree(
                os.path.abspath(self.spec.source), runtime, dirs_exist_ok=True
            )
        else:
            raise RuntimeError(
                f"clone_shadow calculator '{self.spec.name}' requires a source path "
                "or installation commands"
            )
        write_install_stamp(
            runtime,
            fingerprint=fingerprint,
            module=self.spec.name,
            pack_id=pack_id,
            epoch=self.install_epoch,
            extra=(
                {
                    "mode_parent": self.spec.mode_parent,
                    "mode": self.spec.mode_name,
                    "parent_install_fingerprint": self.parent_install_fingerprint(),
                    "parent_install_epoch": self.install_epoch,
                }
                if self.spec.shared_mode_pack
                else None
            ),
        )
        if not self.spec.shared_mode_pack:
            self._installed_shadows.add(pack_id)
        self.last_install_action = "installed"
        if logger is not None:
            try:
                logger.info(
                    f"Calculator install finished -> module={self.spec.name} "
                    f"pack_id={pack_id}\n"
                    f"\t runtime -> {runtime}\n"
                    f"\t stamp -> {install_stamp_path(runtime)}"
                )
            except Exception:
                pass

    def run_initialization(
        self,
        *,
        run_stage: Callable[[list[Mapping[str, Any]], str], None],
    ) -> None:
        """V1 contract: re-run initialization for every sample on a pack."""
        if not self.spec.initialization:
            return
        if self.spec.clone_shadow:
            pack_id = str(self.pack_id or "")
            if not pack_id:
                raise RuntimeError("clone_shadow initialization requires pack_id")
        run_stage(list(self.spec.initialization), "initialize")

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
        logger: Any | None = None,
    ) -> None:
        if self.spec.clone_shadow:
            self.ensure_shadow_installed(run_stage=run_stage, logger=logger)
            self.run_initialization(run_stage=run_stage)
        else:
            self.ensure_symlink_runtime(sample_info)
            self.run_initialization(run_stage=run_stage)


__all__ = [
    "INSTALL_STAMP_BASENAME",
    "INSTALL_STAMP_SCHEMA",
    "INSTALL_CONTROL_BASENAME",
    "INSTALL_CONTROL_SCHEMA",
    "RuntimePreparer",
    "build_install_fingerprint",
    "build_parent_install_fingerprint",
    "force_calc_install_requested",
    "install_control_path",
    "read_install_control",
    "write_install_control",
    "prepare_install_controls",
    "refresh_install_control_summaries",
    "install_stamp_path",
    "read_install_stamp",
    "write_install_stamp",
]
