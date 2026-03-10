#!/usr/bin/env python3
from __future__ import annotations

import os
import json
from pathlib import Path
from typing import Dict, Optional

_TASK_ROOT_ENV_VARS = ("JARVIS_HEP_TASK_ROOT", "JHEP_TASK_ROOT")
_PROJECT_MARKERS = (".jarvis-project.json", "jarvis.project.yaml")


class Base():
    def __init__(self) -> None:
        package_root = os.path.abspath(os.path.dirname(__file__))
        task_root = self._resolve_initial_task_root(default=os.getcwd())
        self.path: Dict[str, str]   = {
            "src_root": package_root,
            "package_root": package_root,
            "task_root": task_root,
            # Backward-compatible alias for legacy code/config.
            "jpath": task_root,
        }
        self.path['args_info'] = self.decode_path("&SRC/card/argparser.json")
        self.path['logger_config_path'] = self.decode_path("&SRC/card/jarvis_logging_config.yaml")
        self.schemablock = {}
        self.load_schema_paths()

    def _resolve_initial_task_root(self, default: str) -> str:
        env_root = self._env_task_root()
        if env_root:
            return env_root
        return os.path.abspath(os.path.expanduser(default))

    def _env_task_root(self) -> Optional[str]:
        for env_name in _TASK_ROOT_ENV_VARS:
            value = os.getenv(env_name, "").strip()
            if value:
                return os.path.abspath(os.path.expanduser(value))
        return None

    def _infer_task_root_from_config_dir(self, config_dir: str) -> str:
        path = Path(config_dir).expanduser().resolve()
        # Preferred project-root detection: explicit marker file in current or ancestor dirs.
        for candidate in [path, *path.parents]:
            for marker in _PROJECT_MARKERS:
                if (candidate / marker).exists():
                    return str(candidate)
        # Keep the historical `project/bin/*.yaml` style by mapping any ancestor
        # named `bin` to its parent as task root.
        if path.name.lower() == "bin":
            return str(path.parent)
        for parent in path.parents:
            if parent.name.lower() == "bin":
                return str(parent.parent)
        return str(path)

    def configure_runtime_context(self, config_path: Optional[str] = None, cwd: Optional[str] = None) -> None:
        config_file = None
        config_dir = None
        if config_path:
            config_file = os.path.abspath(os.path.expanduser(config_path))
            config_dir = os.path.dirname(config_file)

        env_root = self._env_task_root()
        if config_dir:
            # An explicit config path should re-anchor runtime paths to that
            # project instead of inheriting stale task-root env from a prior run.
            task_root = self._infer_task_root_from_config_dir(config_dir)
        elif env_root:
            task_root = env_root
        else:
            task_root = os.path.abspath(os.path.expanduser(cwd or os.getcwd()))

        self.path["task_root"] = task_root
        self.path["jpath"] = task_root
        # Keep subprocess/module-level helpers consistent even if they don't share this Base instance.
        os.environ["JARVIS_HEP_TASK_ROOT"] = task_root
        os.environ["JHEP_TASK_ROOT"] = task_root
        if config_file:
            self.path["config_file"] = config_file
        if config_dir:
            self.path["config_dir"] = config_dir

    def _is_uri(self, path: str) -> bool:
        return "://" in path

    def decode_path(self, path, *, base_dir: Optional[str] = None): 
        """
        Resolves special markers in the provided path.

        Parameters:
        - path: The path to be resolved.

        Returns:
        - The resolved full path.
        """
        if path is None or not isinstance(path, str):
            return path

        normalized = path.replace("\\", "/")
        if (normalized.startswith("&SRC/") or normalized.startswith("&J/")) and "/src/card/" in normalized:
            raise ValueError(
                "Legacy card path prefix is no longer supported: "
                f"{path}. Use '&SRC/card/...' in packaged runtime."
            )

        source_root = self.path.get("package_root") or self.path.get("src_root", "")
        task_root = self.path.get("task_root") or self.path.get("jpath") or os.getcwd()

        # Internal project marker: resolve against installed jarvishep package root.
        if "&SRC" in path:
            path = path.replace("&SRC", source_root)

        # Runtime project marker: resolve against external task/workspace root.
        if "&J" in path:
            path = path.replace("&J", task_root)

        # Replace the user home directory marker ~
        if "~" in path:
            path = os.path.expanduser(path)

        if not self._is_uri(path) and not os.path.isabs(path):
            if base_dir:
                anchor = base_dir
            elif path.startswith("./") or path.startswith("../"):
                anchor = self.path.get("config_dir") or task_root
            else:
                anchor = task_root or self.path.get("config_dir")
            path = os.path.abspath(os.path.join(anchor, path))

        return path  

    def load_schema_paths(self) -> None:
        self.path['preference'] = self.decode_path("&SRC/card/preference.json")
        with open(self.path['preference'], 'r') as f1:
            js = json.loads(f1.read())
            for kk, vv in js["Schema"].items():
                js["Schema"][kk] = self.decode_path(vv)
            for kk, vv in js['SchemaBlock'].items():
                js["SchemaBlock"][kk] = f"file://{self.decode_path(vv)}"
            
            self.path.update(js['Schema'])
            self.schemablock.update(js['SchemaBlock'])
            self.path['logo'] = self.decode_path(js['logo'])


    def manage_directories(self, base_path):
        """
        Manages numbered directories within a given base directory. If the highest-numbered directory contains
        more than 200 files, a new directory with the next number is created.

        Args:
        base_path (str): The path to the base directory to manage.
        """
        # Step 1: Identify all numbered directories
        numbered_dirs = [d for d in os.listdir(base_path) if d.isdigit() and os.path.isdir(os.path.join(base_path, d))]

        if not numbered_dirs:
            # If no numbered directories exist, create the first one and exit
            os.makedirs(os.path.join(base_path, '1'), exist_ok=True)
            # print("Created directory '1' as no numbered directories existed.")
            return os.path.join(base_path, '1')

        # Step 2: Find the highest-numbered directory
        max_number = max([int(num) for num in numbered_dirs])
        max_dir = os.path.join(base_path, str(max_number))

        # Step 3: Check the number of files in the highest-numbered directory
        file_count = len([name for name in os.listdir(max_dir) if os.path.isdir(os.path.join(max_dir, name))])

        # Step 4: Create a new directory if necessary
        if file_count >= 200:
            new_dir_number = max_number + 1
            new_dir_path = os.path.join(base_path, str(new_dir_number))
            os.makedirs(new_dir_path, exist_ok=True)
            # print(f"Created new directory '{new_dir_number}' due to file count greater than 200 in directory '{max_number}'.")
            return new_dir_path
        else:
            # print(f"\nDirectory '{max_dir}' contains {file_count} files, no new directory created.\n")
            return max_dir
