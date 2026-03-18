#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict

import dill
import numpy as np
from loguru import logger as loguru_logger


def make_json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): make_json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [make_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


@dataclass
class RLSamplerPaths:
    task_root: str
    scan_name: str
    outputs_root: str
    checkpoints_root: str
    logs_root: str
    sampler_log: str

    def as_dict(self) -> Dict[str, str]:
        return asdict(self)


class RLStateSaver:
    def __init__(self, root_dir: str, logger=None) -> None:
        self.root_dir = os.path.abspath(root_dir)
        self.logger = logger
        os.makedirs(self.root_dir, exist_ok=True)

    def save_checkpoint(self, name: str, payload: Dict[str, Any]) -> str:
        path = os.path.join(self.root_dir, f"{name}.pt")
        with open(path, "wb") as handle:
            dill.dump(payload, handle)
        if self.logger is not None:
            self.logger.info(f"RL checkpoint saved -> {path}")
        return path

    def load_checkpoint(self, path_or_name: str, default: Any = None) -> Any:
        path = path_or_name
        if not os.path.isabs(path):
            path = os.path.join(self.root_dir, f"{path_or_name}.pt")
        if not os.path.exists(path):
            return default
        with open(path, "rb") as handle:
            payload = dill.load(handle)
        if self.logger is not None:
            self.logger.info(f"RL checkpoint loaded -> {path}")
        return payload

    def save_json(self, name: str, payload: Dict[str, Any]) -> str:
        path = os.path.join(self.root_dir, f"{name}.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(make_json_safe(payload), handle, indent=2, sort_keys=True)
        if self.logger is not None:
            self.logger.info(f"RL metadata saved -> {path}")
        return path


class RLDataWriter:
    DEFAULT_STREAMS = {
        "decision_metrics": "decision_metrics.jsonl",
        "training_metrics": "training_metrics.jsonl",
        "state_action_reward": "state_action_reward.jsonl",
    }

    def __init__(self, root_dir: str, logger=None) -> None:
        self.root_dir = os.path.abspath(root_dir)
        self.logger = logger
        os.makedirs(self.root_dir, exist_ok=True)
        self._stream_paths = {
            name: os.path.join(self.root_dir, filename)
            for name, filename in self.DEFAULT_STREAMS.items()
        }
        self._manifest: Dict[str, Any] = {
            "streams": dict(self._stream_paths),
            "writes": {name: 0 for name in self._stream_paths},
        }

    def append(self, stream: str, payload: Dict[str, Any]) -> str:
        if stream not in self._stream_paths:
            raise KeyError(f"Unknown RL data stream: {stream}")
        path = self._stream_paths[stream]
        with open(path, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(make_json_safe(payload), sort_keys=True) + "\n")
        self._manifest["writes"][stream] = int(self._manifest["writes"].get(stream, 0)) + 1
        return path

    def write_manifest(self, extra: Dict[str, Any] | None = None) -> str:
        manifest = dict(self._manifest)
        if extra:
            manifest.update(make_json_safe(extra))
        path = os.path.join(self.root_dir, "manifest.json")
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
        return path


class RLSamplerFileLogger:
    def __init__(self, module_name: str, file_path: str) -> None:
        self.module_name = str(module_name)
        self.file_path = os.path.abspath(file_path)
        self._handler_id: int | None = None
        os.makedirs(os.path.dirname(self.file_path), exist_ok=True)

    def bind(self):
        def _filter(record):
            return record["extra"].get("module") == self.module_name

        if self._handler_id is None:
            self._handler_id = loguru_logger.add(
                self.file_path,
                format="{time:MM-DD HH:mm:ss.SSS} | {level} | {message}",
                level="DEBUG",
                rotation=None,
                retention=None,
                filter=_filter,
            )
        return loguru_logger.bind(
            module=self.module_name,
            to_console=False,
            Jarvis=True,
            _log_domain="jarvis_hep",
        )

    def close(self) -> None:
        if self._handler_id is not None:
            try:
                loguru_logger.remove(self._handler_id)
            except Exception:
                pass
            self._handler_id = None


class RLSamplerBase:
    def _resolve_rl_task_root(self) -> str:
        task_root = (getattr(self, "path", {}) or {}).get("task_root")
        if task_root:
            return os.path.abspath(task_root)

        task_result_dir = (
            (getattr(self, "info", {}) or {})
            .get("sample", {})
            .get("task_result_dir", os.getcwd())
        )
        task_result_dir = os.path.abspath(task_result_dir)
        parent = os.path.dirname(task_result_dir)
        if os.path.basename(parent) == "outputs":
            return os.path.dirname(parent)
        return parent

    def _resolve_rl_scan_name(self) -> str:
        task_result_dir = (
            (getattr(self, "info", {}) or {})
            .get("sample", {})
            .get("task_result_dir", "")
        )
        task_result_dir = os.path.abspath(task_result_dir) if task_result_dir else ""
        if task_result_dir:
            name = os.path.basename(task_result_dir.rstrip(os.sep))
            if name:
                return name
        config = getattr(self, "config", {}) or {}
        scan_cfg = config.get("Scan", {}) if isinstance(config, dict) else {}
        return str(scan_cfg.get("name", getattr(self, "method", "rlsampler")))

    def init_rl_runtime(self, sampler_slug: str) -> None:
        task_root = self._resolve_rl_task_root()
        scan_name = self._resolve_rl_scan_name()
        outputs_root = os.path.join(task_root, "outputs", scan_name, sampler_slug)
        checkpoints_root = os.path.join(task_root, "checkpoints", scan_name, sampler_slug)
        logs_root = os.path.join(task_root, "logs", scan_name)
        sampler_log = os.path.join(logs_root, f"{sampler_slug}.log")

        os.makedirs(outputs_root, exist_ok=True)
        os.makedirs(checkpoints_root, exist_ok=True)
        os.makedirs(logs_root, exist_ok=True)

        self.rl_paths = RLSamplerPaths(
            task_root=task_root,
            scan_name=scan_name,
            outputs_root=outputs_root,
            checkpoints_root=checkpoints_root,
            logs_root=logs_root,
            sampler_log=sampler_log,
        )
        self.rl_data_writer = RLDataWriter(outputs_root)
        self.rl_file_logger = RLSamplerFileLogger(
            module_name=f"Jarvis-HEP.{getattr(self, 'method', sampler_slug)}.RL",
            file_path=sampler_log,
        )
        self.rl_logger = self.rl_file_logger.bind()
        self.rl_state_saver = RLStateSaver(checkpoints_root, logger=self.rl_logger)

        self.info.setdefault("rl", {})
        self.info["rl"].update(self.rl_paths.as_dict())

        self.rl_logger.warning(
            "RL runtime initialized -> "
            f"outputs={outputs_root} checkpoints={checkpoints_root} log={sampler_log}"
        )

    def rl_log(self, level: str, message: str) -> None:
        main_logger = getattr(self, "logger", None)
        rl_logger = getattr(self, "rl_logger", None)
        log_method = str(level).strip().lower()
        if main_logger is not None and hasattr(main_logger, log_method):
            getattr(main_logger, log_method)(message)
        if rl_logger is not None and hasattr(rl_logger, log_method):
            getattr(rl_logger, log_method)(message)

    def write_rl_record(self, stream: str, payload: Dict[str, Any]) -> str:
        return self.rl_data_writer.append(stream, payload)

    def save_rl_checkpoint(self, name: str, payload: Dict[str, Any]) -> str:
        return self.rl_state_saver.save_checkpoint(name, payload)

    def save_rl_metadata(self, name: str, payload: Dict[str, Any]) -> str:
        return self.rl_state_saver.save_json(name, payload)

    def close_rl_runtime(self, extra_manifest: Dict[str, Any] | None = None) -> None:
        if getattr(self, "rl_data_writer", None) is not None:
            manifest_path = self.rl_data_writer.write_manifest(extra=extra_manifest or {})
            if getattr(self, "rl_logger", None) is not None:
                self.rl_logger.info(f"RL manifest saved -> {manifest_path}")
        if getattr(self, "rl_file_logger", None) is not None:
            self.rl_file_logger.close()
