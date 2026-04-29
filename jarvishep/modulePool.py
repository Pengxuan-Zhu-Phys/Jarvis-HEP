#!/usr/bin/env python3
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from jarvishep.Module.calculator import CalculatorModule
from copy import deepcopy
import json
import os
import threading
from typing import Dict, Any

from loguru import logger
import asyncio


class ModulePool:
    def __init__(self, module, max_workers=2):
        self.name = module.name 
        self.config = module.config 
        # Keep a lightweight executor for background installation tasks only.
        # Runtime module execution stays on the caller thread (Factory worker)
        # to avoid nested executor scheduling in the hot path.
        self.executor = ThreadPoolExecutor(max_workers=max_workers)
        self.module_name = module.name 
        self.type = module.type 
        self.instances = []
        self.id_counter = 0 
        self.instances_info_file = self.get_instances_info_file_path()
        self._info_loaded = False 
        self._funcs = {}
        self._subprocess_scheduler = None
        self._io_manager = None
        self._run_summary_collector = None

        # Single-process, multi-thread safety
        self._inst_lock = threading.Lock()   # protects instances list + is_busy transitions
        self._io_lock = threading.Lock()     # protects info-file writes

    def _maybe_call(self, instance, method_name, *args, **kwargs):
        method = getattr(instance, method_name, None)
        if callable(method):
            method(*args, **kwargs)

    def set_logger(self):
        self.logger = logger.bind(
            module=f"Jarvis-HEP.Workflow.{self.name}",
            to_console=True,
            Jarvis=True,
            _log_domain="jarvis_hep",
        )
        self.load_installed_instances()

    def create_instance(self):
        self.id_counter += 1
        id = f"{self.id_counter:03}"
        config = deepcopy(self.config)
        self.logger.info(f"Create the {id} Instance for Module {self.name}")
        instance = CalculatorModule(self.name, config=config)
        instance.assign_ID(id)
        instance.is_installed = False
        instance.is_busy = False
        instance.installation_event = threading.Event()
        instance._funcs = self.funcs
        self._maybe_call(instance, "set_subprocess_scheduler", self._subprocess_scheduler)
        self._maybe_call(instance, "set_io_manager", self._io_manager)
        self._maybe_call(instance, "set_run_summary_collector", self._run_summary_collector)
        self.instances.append(instance)
        return instance

    def reload_instance(self, pack_id, pack_info):
        # Keep id_counter in sync with existing numeric PackIDs
        try:
            self.id_counter = max(self.id_counter, int(str(pack_id)))
        except Exception:
            pass

        config = deepcopy(self.config)
        self.logger.info(f"Re-loading the {pack_id} Instance for Module {self.name}")
        instance = CalculatorModule(self.name, config=config)
        instance.assign_ID(pack_id)

        # Persisted state: only installation matters for reload.
        instance.is_installed = bool((pack_info or {}).get("is_installed", True))

        # Runtime-only state (never persisted)
        instance.is_busy = False
        instance.installation_event = threading.Event()
        if instance.is_installed:
            instance.installation_event.set()

        instance._funcs = self.funcs
        self._maybe_call(instance, "set_subprocess_scheduler", self._subprocess_scheduler)
        self._maybe_call(instance, "set_io_manager", self._io_manager)
        self._maybe_call(instance, "set_run_summary_collector", self._run_summary_collector)
        self.instances.append(instance)
        return instance

    def export_blueprint(self) -> Dict[str, Any]:
        instances = {}
        for instance in self.instances:
            pack_id = getattr(instance, "PackID", None)
            if pack_id is None:
                continue
            instances[str(pack_id)] = {
                "is_installed": bool(getattr(instance, "is_installed", False)),
            }
        return {
            "name": self.name,
            "type": self.type,
            "config": deepcopy(self.config),
            "instances": instances,
        }

    def restore_blueprint(self, blueprint: Dict[str, Any]) -> None:
        if not isinstance(blueprint, dict):
            return

        instances = blueprint.get("instances", {})
        if not isinstance(instances, dict):
            instances = {}

        with self._inst_lock:
            self.instances = []
            self.id_counter = 0
            self._info_loaded = False

        for pack_id in sorted(instances.keys(), key=lambda value: str(value)):
            pack_info = instances.get(pack_id, {})
            self.reload_instance(pack_id=pack_id, pack_info=pack_info)

        self.update_instances_info_file()

    def set_funcs(self, funcs):
        self._funcs = funcs
        for instance in self.instances:
            instance._funcs = funcs
        # self.logger.warning(f"ModulePool {self.name} funcs -> {self.funcs}")

    def set_subprocess_scheduler(self, scheduler):
        self._subprocess_scheduler = scheduler
        for instance in self.instances:
            self._maybe_call(instance, "set_subprocess_scheduler", scheduler)

    def set_io_manager(self, io_manager):
        self._io_manager = io_manager
        for instance in self.instances:
            self._maybe_call(instance, "set_io_manager", io_manager)

    def set_run_summary_collector(self, collector):
        self._run_summary_collector = collector
        for instance in self.instances:
            self._maybe_call(instance, "set_run_summary_collector", collector)

    @property
    def funcs(self):
        return self._funcs

    def execute(self, params, sample_info):
        """Submit parameters to a module instance for calculation and return the calculation results."""

        with self._inst_lock:
            instance = self.get_available_instance()

            if instance is None:
                # If no available instance found, create a new one (installed asynchronously)
                self.logger.warning(
                    f"No availiable instance for Module {self.name} found, trying to install a new one"
                )
                instance = self.create_instance()
                self.executor.submit(self.install_instance, instance)
            # Reserve the instance immediately to avoid double-rent races.
            instance.is_busy = True

        # Wait for installation to complete (install_instance always sets the event)
        if not instance.is_installed:
            instance.installation_event.wait()
            if not instance.is_installed:
                with self._inst_lock:
                    instance.is_busy = False
                raise RuntimeError(
                    f"Installation failed for Module {self.name} instance {getattr(instance, 'PackID', 'UNKNOWN')}"
                )
            # Only now persist installation state (to avoid reloading half-installed instances)
            self.update_instances_info_file()

        try:
            self._maybe_call(instance, "set_subprocess_scheduler", self._subprocess_scheduler)
            self._maybe_call(instance, "set_io_manager", self._io_manager)
            self._maybe_call(instance, "set_run_summary_collector", self._run_summary_collector)
            # Execute directly on current worker thread to avoid nested
            # ThreadPoolExecutor submit/result overhead.
            return instance.execute(params, sample_info)
        finally:
            with self._inst_lock:
                instance.is_busy = False
            # Do NOT write info file here; busy is runtime-only.

    @staticmethod
    def install_instance(instance):
        # Always release waiters, even if installation fails.
        try:
            if not instance.is_installed:
                asyncio.run(instance.install())
                instance.is_installed = True
        except Exception:
            instance.is_installed = False
        finally:
            # Ensure waiters are released
            if hasattr(instance, "installation_event") and instance.installation_event is not None:
                instance.installation_event.set()
        return instance

    def get_available_instance(self):
        # print(self.instances)
        for instance in self.instances:
            if not instance.is_busy and instance.is_installed:
                return instance
        return None

    def _atomic_write_json(self, path: str, data: Dict[str, Any]) -> None:
        """Write JSON atomically to avoid partial/corrupted files."""
        tmp_path = path + ".tmp"
        with open(tmp_path, "w") as f:
            json.dump(data, f, indent=4)
            f.flush()
        os.replace(tmp_path, path)

    def update_instances_info_file(self):
        """Persist only installed instances to support reload and avoid redundant installs."""
        installed_instances_info = {
            instance.PackID: {"is_installed": bool(instance.is_installed)}
            for instance in self.instances
            if getattr(instance, "PackID", None) is not None
        }
        data_to_save = {"installed_instances": installed_instances_info}

        # Atomic write to avoid partial/corrupted JSON under concurrent threads
        with self._io_lock:
            self._atomic_write_json(self.instances_info_file, data_to_save)

    def load_installed_instances(self):
        if self._info_loaded:
            self.logger.warning("Instance information has already been loaded.")
            return
        
        if not os.path.exists(self.instances_info_file):
            self.logger.warning(f"{self.instances_info_file} not found. Initializing a new instances info file.")
            self.init_instances_info_file()
        else: 
            self.logger.warning(f"Trying to loac the installed instance information for mudule -> {self.name}")
            try:
                with open(self.instances_info_file, 'r') as file:
                    installed_instances_info = json.load(file)
                    installed_ids = installed_instances_info.get("installed_instances", {})
                    self.restore_blueprint({"instances": installed_ids})
                    self.logger.warning("Instance reloaded!")
            except FileNotFoundError:
                self.logger.info(f"{self.instances_info_file} not found. Assuming no instances are installed.")
        self._info_loaded = True

    def get_instances_info_file_path(self):
        return os.path.join(self.config['path'].replace("/@PackID", ""), f"{self.name}_instance_info.json")

    def init_instances_info_file(self):
        if not os.path.exists(os.path.dirname(self.instances_info_file)):
            self.logger.warning(f"Init instance info file -> {os.path.dirname(self.instances_info_file)}")
            os.makedirs(os.path.dirname(self.instances_info_file))
        data_to_save = {"installed_instances": {}}
        with self._io_lock:
            self._atomic_write_json(self.instances_info_file, data_to_save)
