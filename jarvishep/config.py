#!/usr/bin/env python3
from __future__ import annotations

import yaml 
import platform
import os, sys 
import re 
import shlex
import subprocess
from pathlib import Path
from urllib.parse import unquote, urlparse
from importlib import metadata as importlib_metadata
from jarvishep.base import Base
import json
from jsonschema.exceptions import ValidationError
from jsonschema.validators import validator_for
from referencing import Registry, Resource
from referencing.jsonschema import DRAFT7
from jarvishep.Module.module import Module

class ConfigLoader(Base):
    def __init__(self) -> None:
        super().__init__()
        self.filepath   = None
        self.config     = None
        self.schema     = None
        self.logger     = None 
        self.legal      = False 
        self.summary    = {}

    def get_cmd_output(self, cmd):
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.stderr:
            self.logger.warning("An error happend when excute: \n -> {} \n\t{}".format(cmd, result.stderr))    
        return result.stdout.strip()

    def set_schema(self, schema) -> None:
        self.schema = schema

    def load_config(self, filepath):
        try:
            self.configure_runtime_context(config_path=filepath)
            self.logger.info(f"Start configuring the input file -> {filepath}")
            with open(filepath, 'r') as file:
                self.config = yaml.safe_load(file)
                self._normalize_optional_sections()
                self.filepath = filepath
                env_reqs = self.config.get("EnvReqs", {}) if isinstance(self.config, dict) else {}
                check_default = env_reqs.get("Check_default_dependences", {}) if isinstance(env_reqs, dict) else {}
                if isinstance(check_default, dict) and bool(check_default.get("required", False)):
                    self.update_dependences()
        except Exception as e: 
            self.logger.error(f"Error: loading config file -> {e}")
            sys.exit(2)

    def _normalize_optional_sections(self):
        """Normalize optional top-level sections before schema validation.

        Current schemas require a Calculators section. We inject an empty default
        so Operas-only tasks remain schema-compatible.

        Directory settings defaults:
        - New preferred location: Scan.sample_directory
        - Legacy compatible location: Directory_Setting
        - Keys:
          - archive_samples: bool, default True
          - archive_prune: bool, default False
        """
        if not isinstance(self.config, dict):
            return
        if "Calculators" not in self.config:
            self.config["Calculators"] = {
                "make_paraller": 4,
                "Modules": [],
            }
        if "Operas" not in self.config:
            self.config["Operas"] = {}

        defaults = {
            "limit": 200,
            "width": 6,
            "archive_samples": True,
            "archive_prune": False,
        }

        legacy_cfg = self.config.get("Directory_Setting")
        if not isinstance(legacy_cfg, dict):
            legacy_cfg = {}

        scan_cfg = self.config.get("Scan")
        if not isinstance(scan_cfg, dict):
            scan_cfg = {}
        scan_dir_cfg = scan_cfg.get("sample_directory")
        if not isinstance(scan_dir_cfg, dict):
            scan_dir_cfg = {}

        merged = dict(defaults)
        merged.update(legacy_cfg)
        merged.update(scan_dir_cfg)

        scan_cfg["sample_directory"] = merged
        self.config["Scan"] = scan_cfg
        # Keep legacy key populated for backward compatibility with older paths.
        self.config["Directory_Setting"] = dict(merged)

    def check_dependency_installed(self) -> None:
        dependencies = self.config["EnvReqs"]
        # Update the environment requirement by default setting 
        if "Check_default_dependences" in dependencies:
            self.update_dependences()
        # Check the system version
        if "OS" in dependencies:
            self.check_OS_requirement(dependencies["OS"])
        # Check the CERN-ROOT setting
        if "CERN_ROOT" in dependencies:
            self.check_ROOT()
        # Check the Python Dependences 
        if "Python" in dependencies:
            self.check_PYTHON_env()
        
    
    def get_sampling_method(self) -> str:
        try: 
            return self.config['Sampling']['Method']
        except:
            self.logger.error(f"Iegal configparser file, No Sampling method founded -> {self.filepath}")
            sys.exit(2)

    def validate_config(self) -> None: 
        validator = ConfigValidator()
        validator.logger = self.logger
        validator.set_config(self.config)
        validator.set_schema(self.schema)
        validator.validate_yaml()

    def check_PYTHON_env(self) -> None: 
        py_env = self.config['EnvReqs']['Python']

        def to_version_tuple(version_text):
            numbers = tuple(int(x) for x in re.findall(r"\d+", str(version_text)))
            return numbers if numbers else (0,)

        def compare_versions(version1, version2):
            v1 = to_version_tuple(version1)
            v2 = to_version_tuple(version2)
            width = max(len(v1), len(v2))
            v1 += (0,) * (width - len(v1))
            v2 += (0,) * (width - len(v2))
            return v1 >= v2

        missing_required = []
        incompatible_required = []

        def check_package_requirement(name, required, min_version):
            try:
                installed_version = importlib_metadata.version(name)
            except importlib_metadata.PackageNotFoundError:
                if required:
                    missing_required.append((name, min_version))
                    self.summary[name] = None
                else:
                    self.logger.warning(f"{name} is optional and not installed")
                    self.summary[name] = None
                return

            if compare_versions(installed_version, min_version):
                self.logger.info(f"{name} is found, version: {installed_version} meets the requirement")
                self.summary[name] = installed_version
                return

            self.summary[name] = installed_version
            if required:
                incompatible_required.append((name, installed_version, min_version))
            else:
                self.logger.warning(
                    f"{name} version: {installed_version} does not meet optional requirement >= {min_version}"
                )
       
        vf = sys.version_info
        py_version = f"{vf.major}.{vf.minor}.{vf.micro}"
        version_matched = re.match(r">=(\d+(?:\.\d+)?)", py_env['version'])
        if version_matched:
            min_version = version_matched.group(1)
            if compare_versions(py_version, min_version):
                self.summary["Python version"] = py_version
                self.logger.info(f"Python-{py_version} found.")
            else:
                self.logger.error(f"Python-{py_version} version does not meet the requirement >= {min_version}.")
                sys.exit(2)
        
        if "Dependencies" in py_env:
            for lib in py_env["Dependencies"]:
                name = lib.get("name")
                required = bool(lib.get("required", False))
                version_spec = str(lib.get("version", "")).strip()
                if not name:
                    self.logger.warning("Encountered Python dependency entry without name; skipping.")
                    continue
                version_match = re.match(r">=\s*(\d+(?:\.\d+)*)", version_spec)
                if not version_match:
                    self.logger.warning(f"Python library {name} has unsupported version spec {version_spec!r}; skipping.")
                    continue
                min_version = version_match.group(1)
                check_package_requirement(name, required, min_version)

        if missing_required or incompatible_required:
            if missing_required:
                for name, min_version in missing_required:
                    self.logger.error(f"Missing required package: {name}>={min_version}")
            if incompatible_required:
                for name, installed_version, min_version in incompatible_required:
                    self.logger.error(
                        f"Incompatible package version: {name}=={installed_version} (requires >= {min_version})"
                    )
            install_targets = [f"{name}>={min_version}" for name, min_version in missing_required]
            install_targets.extend(
                [f"{name}>={min_version}" for name, _, min_version in incompatible_required]
            )
            if install_targets:
                install_cmd = "python3 -m pip install -U " + " ".join(f"'{item}'" for item in install_targets)
                self.logger.error(f"Install/upgrade required dependencies via: {install_cmd}")
            sys.exit(2)

    def check_ROOT(self) -> None:
        def compare_versions(version1, version2):
            version1 = version1.replace("/", '.')
            v1 = tuple(map(int, version1.split('.')))
            v2 = tuple(map(int, version2.split('.')))
            return v1 >= v2
        
        root_req = self.config['EnvReqs']['CERN_ROOT']
        if root_req['required']:
            try:
                result = subprocess.run("root-config --version", shell=True, capture_output=True, text=True)
                self.summary['ROOT'] = True
            except subprocess.CalledProcessError:
                self.logger.error("CERN ROOT is not installed or root-config is not in the PATH")
                self.summary['ROOT'] = False
            except FileNotFoundError:
                self.logger.error("CERN ROOT is not installed or root-config is not in the PATH")
                self.summary['ROOT'] = False

            if self.summary['ROOT']:
                if "get_path_command" in root_req:
                    root_prefix = self.get_cmd_output(root_req['get_path_command'])
                elif "path" in root_req:
                    root_prefix = self.decode_path(root_req['path'])
                self.summary['ROOT path'] = root_prefix

                if "version" in root_req:
                    root_version = self.get_cmd_output("root-config --version")
                    version_matched = re.match(r">=(\d+(?:\.\d+)?)", root_req['version'])
                    if version_matched:
                        min_version = version_matched.group(1)
                        if compare_versions(root_version, min_version):
                            self.summary["ROOT version"] = root_version
                            self.logger.info(f"CERN ROOT-{root_version} found.")
                        else:
                            self.logger.error(f"ROOT Version-{root_version} version does not meet the requirement >= {min_version}.")
                            sys.exit(2)
                
                if "Dependencies" in root_req:
                    for dep in root_req['Dependencies']:
                        if dep['required']:
                            result = subprocess.run(dep['check_command'], shell=True, capture_output=True, text=True)
                            if not result.stdout.strip() == dep['expected_output']:
                                self.logger.error("ROOT-{} not founded, please check the CERN ROOT installation".format(dep['name']))
                                self.summary['ROOT-{}'.format(dep['name'])] = False
                                self.summary['ROOT'] = False
                            else:
                                self.summary['ROOT-{}'.format(dep['name'])] = True 

            if not self.summary['ROOT']:
                self.logger.error("ROOT installation is not meeting the requirements!")
                sys.exit(2)
            else:
                self.logger.info("CERN ROOT meets the requirements!")

    def check_OS_requirement(self, os_requirement) -> None: 
        # def compare_versions(version1, version2):
        #     v1 = tuple(map(int, version1.split('.')))
        #     v2 = tuple(map(int, version2.split('.')))
        #     return v1 >= v2

        def compare_versions(version1, version2):
            # Extract numeric parts of the version strings
            numeric_version1 = re.findall(r'\d+', version1)
            numeric_version2 = re.findall(r'\d+', version2)

            # Convert numeric strings to integers
            v1 = tuple(map(int, numeric_version1))
            v2 = tuple(map(int, numeric_version2))

            # Compare version tuples
            return v1 >= v2



        current_os = platform.system()
        current_version = platform.release()
        self.summary['OS'] = "{}-{}".format(current_os, current_version)

        os_matched = False
        for os_req in os_requirement:
            if os_req['name'].lower() == current_os.lower(): 
                version_matched = re.match(r">=(\d+(?:\.\d+)?)", os_req['version'])
                if version_matched:
                    min_version = version_matched.group(1)
                    if compare_versions(current_version, min_version):
                        os_matched = True
                        self.logger.info(f"{current_os} version {current_version} meets the requirement >= {min_version}.")
                    else:
                        os_matched = False
                        self.logger.critical(f"{current_os} version {current_version} does not meet the requirement >= {min_version}.")
                break
        if not os_matched:
            sys.exit(2)

    def update_dependences(self) -> None:
        if "required" in self.config['EnvReqs']['Check_default_dependences'] and "default_yaml_path" in self.config['EnvReqs']['Check_default_dependences']: 
            if self.config['EnvReqs']['Check_default_dependences']['required']:
                try:
                    self.config['EnvReqs']['Check_default_dependences']['default_yaml_path'] = self.decode_path(self.config['EnvReqs']['Check_default_dependences']['default_yaml_path'])
                    with open(self.decode_path(self.config['EnvReqs']['Check_default_dependences']['default_yaml_path']), 'r') as file:
                        default_env = yaml.safe_load(file)
                        self.config['EnvReqs'].update(default_env['EnvReqs'])
                except FileExistsError:
                    self.logger.error("Jarvis-HEP load file error: {} not found".format(self.config['EnvReqs']['default_yaml_path']))
        self.logger.info("Updating the Environment Requirements from default setting file \n\t{}".format(self.config['EnvReqs']['Check_default_dependences']['default_yaml_path']))

    def get_modules(self):
        sampling_cfg = self.config.get('Sampling', {}) if isinstance(self.config, dict) else {}
        parameters = sampling_cfg.get('Variables', [])
        if parameters is None:
            parameters = []
        if not isinstance(parameters, list):
            parameters = []
        modules = {
            "Parameter": parameters
        }
        
        if sampling_cfg.get("Nuisance", False): 
            modules['Nuisance'] = sampling_cfg['Nuisance']['Variables']
        
        if "LibDeps" in self.config:
            if self.config["LibDeps"]["Modules"]:
                self.analysis_Library()
                modules['Library'] = self.config["LibDeps"]["Modules"]

        calculators = (self.config.get("Calculators", {}) or {}).get("Modules", []) or []
        operas = (self.config.get("Operas", {}) or {}).get("Modules", []) or []
        if calculators:
            self.analysis_calculator()
            modules['Calculator'] = self.config["Calculators"]["Modules"]
        if operas:
            self.analysis_operas()
            modules["Operas"] = self.config["Operas"]["Modules"]

        if not calculators and not operas:
            self.logger.error("Invalid config: at least one of Calculators.Modules or Operas.Modules must be non-empty.")
            sys.exit(2)
        return modules
        
    def decode_lib_command_via_config(self, command, pwd, pack):
        cmd = self._resolve_command_text(command=command, scopes=(self.config, pack))
        current_cwd = self.decode_path(pwd)
        cmdd = {"cmd": cmd, "cwd": current_cwd}
        next_cwd = self._next_cwd_from_command(command=cmd, current_cwd=current_cwd)
        return cmdd, next_cwd

    def analysis_Library(self):
        def analysis_path(lib):
            lib['installation']['path'] = self.decode_path(lib['installation']['path'])
            lib['installation']['source'] = self.decode_path(lib['installation']['source'])
            return lib 

        def analysis_commands(lib):
            cmds = lib['installation']['commands']
            cwd  = self.config['LibDeps']['path']
            commands = []
            for command in cmds:
                cmd, cwd = self.decode_lib_command_via_config(command, cwd, lib['installation'])
                commands.append(cmd)
            lib['installation']['commands'] = commands
            return lib 

        self.config["LibDeps"]['path'] = self.decode_path(self.config["LibDeps"]['path'])
        libs_config = self.config["LibDeps"]["Modules"]
        for lib in libs_config:
            lib.update(analysis_path(lib))
            lib.update(analysis_commands(lib))

    def decode_calc_command_via_config(self, command, pwd, pack):
        cmd = self._resolve_command_text(command=command, scopes=(self.config, pack))
        current_cwd = self.decode_path(pwd)
        cmdd = {"cmd": cmd, "cwd": current_cwd}
        next_cwd = self._next_cwd_from_command(command=cmd, current_cwd=current_cwd)
        return cmdd, next_cwd

    def analysis_calculator(self):
        def analysis_path(calc):
            if 'path' in calc:
                calc['path'] = self.decode_path(calc['path'])
            if "source" in calc:
                calc['source'] = self.decode_path(calc['source'])
            return calc 

        def analysis_install_commands(calc):
            cmds = calc['installation']
            cwd  = self.config['Calculators']['path']
            commands = []
            for command in cmds:
                cmd, cwd = self.decode_calc_command_via_config(command, cwd, calc)
                commands.append(cmd)
            calc['installation'] = commands
            return calc 
        
        def analysis_initialization_commands_single(calc):
            cmds = calc['initialization']
            cwd  = calc['path']
            commands = []
            for command in cmds:
                cmd, cwd = self.decode_calc_command_via_config(command, cwd, calc)
                commands.append(cmd)
            calc['initialization'] = commands
            return calc
        
        def analysis_execution_single(calc):
            calc_setting = calc['execution']

            calc_setting['path'] = self.decode_path(calc_setting['path'])

            commands = []
            cwd = calc_setting['path']
            for command in calc_setting['commands']:
                cmd, cwd = self.decode_calc_command_via_config(command, cwd, calc_setting)
                commands.append(cmd)
            calc_setting['commands'] = commands

            for ipf in calc_setting['input']:
                ipf['path'] = self.decode_path(ipf['path'])

            for opf in calc_setting['output']:
                opf['path'] = self.decode_path(opf['path'])

            return calc

        calc_root = self.config.get("Calculators", {}) or {}
        calc_root.setdefault("make_paraller", 4)
        if "path" in calc_root:
            calc_root['path'] = self.decode_path(calc_root['path'])
        else:
            calc_root["path"] = self.decode_path("&J/calculators/runtime/program")
        self.config["Calculators"] = calc_root
        calc_config = calc_root.get('Modules', []) or []
        for calc in calc_config:
            calc.update(analysis_path(calc))
            calc.update(analysis_install_commands(calc))
            if "modes" in calc:
                self.logger.info(f"Loading configuration of Module {calc['name']} in multiple modes")
            else:
                self.logger.info(f"Loading configuration of Module {calc['name']} in single mode")
                calc['modes'] = False
                calc.update(analysis_initialization_commands_single(calc))
                calc.update(analysis_execution_single(calc))

    def analysis_operas(self):
        """Normalize Operas section into a workflow-ready structure."""
        def normalize_call_mode(value, *, context: str) -> str:
            if value is None:
                return "call"
            if not isinstance(value, str):
                self.logger.error(
                    f"Invalid Operas call_mode in {context}: expected string 'call' or 'acall', got {value!r}"
                )
                sys.exit(2)
            mode = value.strip().lower()
            if mode not in {"call", "acall"}:
                self.logger.error(
                    f"Invalid Operas call_mode in {context}: '{value}'. Expected 'call' or 'acall'."
                )
                sys.exit(2)
            return mode

        operas_root = self.config.get("Operas", {}) or {}
        operas_root.setdefault("make_paraller", 4)
        modules = operas_root.get("Modules", []) or []
        for module in modules:
            if not isinstance(module, dict):
                self.logger.error(f"Invalid Operas module entry: {module}")
                sys.exit(2)
            if not isinstance(module.get("name"), str) or not module.get("name"):
                self.logger.error(f"Invalid Operas module name: {module}")
                sys.exit(2)
            if not isinstance(module.get("operator"), str) or not module.get("operator"):
                self.logger.error(f"Invalid Operas operator for module '{module.get('name', '<unknown>')}': {module}")
                sys.exit(2)
            module["call_mode"] = normalize_call_mode(
                module.get("call_mode", "call"),
                context=f"Operas.Modules[{module.get('name', '<unknown>')}]",
            )
            if "required_modules" not in module:
                module["required_modules"] = []
            if isinstance(module["required_modules"], str):
                module["required_modules"] = [module["required_modules"]]
            module.setdefault("input", [])
            module.setdefault("kwargs", {})
            if not isinstance(module["input"], list):
                self.logger.error(
                    f"Invalid Operas module input for '{module.get('name', '<unknown>')}': expected list."
                )
                sys.exit(2)
            if not isinstance(module["kwargs"], dict):
                self.logger.error(
                    f"Invalid Operas module kwargs for '{module.get('name', '<unknown>')}': expected object."
                )
                sys.exit(2)
            normalized_input = []
            for item in module["input"]:
                if isinstance(item, str):
                    normalized_input.append({"name": item})
                elif isinstance(item, dict):
                    if not isinstance(item.get("name"), str) or not item.get("name"):
                        self.logger.error(
                            f"Invalid Operas input mapping in module '{module.get('name', '<unknown>')}': {item}"
                        )
                        sys.exit(2)
                    if "expression" in item and not isinstance(item["expression"], str):
                        self.logger.error(
                            f"Invalid Operas input expression in module '{module.get('name', '<unknown>')}': {item}"
                        )
                        sys.exit(2)
                    normalized_input.append(item)
                else:
                    self.logger.error(
                        f"Invalid Operas input item in module '{module.get('name', '<unknown>')}': {item}"
                    )
                    sys.exit(2)
            module["input"] = normalized_input

            normalized_output = []
            for item in module.get("output", []) or []:
                if isinstance(item, str):
                    normalized_output.append({"name": item, "entry": item})
                elif isinstance(item, dict):
                    name = item.get("name")
                    if not name:
                        self.logger.error(f"Invalid Operas output mapping in module '{module.get('name', '<unknown>')}': {item}")
                        sys.exit(2)
                    normalized_output.append(
                        {
                            "name": name,
                            "entry": item.get("entry", name),
                        }
                    )
                else:
                    self.logger.error(f"Invalid Operas output item in module '{module.get('name', '<unknown>')}': {item}")
                    sys.exit(2)
            if not normalized_output:
                self.logger.error(
                    f"Invalid Operas module '{module.get('name', '<unknown>')}': output mapping cannot be empty."
                )
                sys.exit(2)
            module["output"] = normalized_output
        operas_root["Modules"] = modules
        self.config["Operas"] = operas_root

    def get_worker_parallel(self) -> int:
        """Return unified worker parallelism from Calculators/Operas sections."""
        candidates = []
        c_para = (self.config.get("Calculators", {}) or {}).get("make_paraller")
        o_para = (self.config.get("Operas", {}) or {}).get("make_paraller")
        for val in (c_para, o_para):
            if isinstance(val, int) and val > 0:
                candidates.append(val)
        return max(candidates) if candidates else 4

    def get_subprocess_runtime_options(self, worker_parallel: int | None = None) -> dict:
        """Return normalized runtime knobs for async subprocess scheduler.

        Preferred user config path (schema-friendly, optional):
        - Runtime.Subprocess
        Legacy-compatible fallback:
        - Runtime
        """
        cpu_count = int(os.cpu_count() or 4)
        io_default = min(128, max(4, cpu_count * 2))
        worker_parallel = int(worker_parallel or self.get_worker_parallel())
        default_max_concurrency = max(1, min(worker_parallel, io_default))

        hard_max_concurrency = 256

        defaults = {
            "max_concurrency": default_max_concurrency,
            "max_pending": max(128, default_max_concurrency * 16),
            "per_task_timeout_sec": None,
            "progress_interval_sec": 5.0,
            "log_policy": "logger",
            "diagnostics_enabled": False,
            "diagnostics_interval_sec": 10.0,
            "terminate_grace_sec": 5.0,
        }

        runtime_cfg = self.config.get("Runtime", {}) if isinstance(self.config, dict) else {}
        if not isinstance(runtime_cfg, dict):
            runtime_cfg = {}
        sub_cfg = runtime_cfg.get("Subprocess", runtime_cfg)
        if not isinstance(sub_cfg, dict):
            sub_cfg = {}

        opts = dict(defaults)
        for kk in defaults.keys():
            if kk in sub_cfg and sub_cfg.get(kk) is not None:
                opts[kk] = sub_cfg[kk]

        try:
            opts["max_concurrency"] = max(1, int(opts["max_concurrency"]))
        except Exception:
            opts["max_concurrency"] = default_max_concurrency
        if opts["max_concurrency"] > hard_max_concurrency:
            if self.logger is not None:
                self.logger.warning(
                    f"Runtime.Subprocess.max_concurrency={opts['max_concurrency']} exceeds hard cap {hard_max_concurrency}; clamp to cap."
                )
            opts["max_concurrency"] = hard_max_concurrency
        try:
            opts["max_pending"] = max(1, int(opts["max_pending"]))
        except Exception:
            opts["max_pending"] = defaults["max_pending"]

        if opts.get("per_task_timeout_sec") is not None:
            try:
                opts["per_task_timeout_sec"] = max(0.1, float(opts["per_task_timeout_sec"]))
            except Exception:
                opts["per_task_timeout_sec"] = None

        try:
            opts["progress_interval_sec"] = max(0.2, float(opts["progress_interval_sec"]))
        except Exception:
            opts["progress_interval_sec"] = defaults["progress_interval_sec"]

        mode = str(opts.get("log_policy", "logger")).strip().lower()
        if mode not in {"logger", "file", "quiet", "tee-limited"}:
            mode = "logger"
        opts["log_policy"] = mode

        opts["diagnostics_enabled"] = bool(opts.get("diagnostics_enabled", False))
        try:
            opts["diagnostics_interval_sec"] = max(0.2, float(opts["diagnostics_interval_sec"]))
        except Exception:
            opts["diagnostics_interval_sec"] = defaults["diagnostics_interval_sec"]
        try:
            opts["terminate_grace_sec"] = max(0.1, float(opts["terminate_grace_sec"]))
        except Exception:
            opts["terminate_grace_sec"] = defaults["terminate_grace_sec"]

        return opts

    def get_operas_function_whitelist(self) -> list[dict]:
        funcs = (self.config.get("Operas", {}) or {}).get("Functions", []) or []
        if not isinstance(funcs, list):
            self.logger.error("Invalid config: Operas.Functions must be a list.")
            sys.exit(2)
        normalized = []
        for item in funcs:
            if not isinstance(item, dict):
                self.logger.error(f"Invalid Operas.Functions item: {item}")
                sys.exit(2)
            name = item.get("name")
            if not isinstance(name, str) or not name:
                self.logger.error(f"Invalid Operas function name: {item}")
                sys.exit(2)
            alias = item.get("alias", name.split(":")[-1])
            if not isinstance(alias, str) or not alias:
                self.logger.error(f"Invalid Operas alias: {item}")
                sys.exit(2)
            inputs = item.get("inputs", [])
            if not isinstance(inputs, list) or any(not isinstance(x, str) for x in inputs):
                self.logger.error(f"Invalid Operas function inputs: {item}")
                sys.exit(2)
            normalized.append(
                {
                    "name": name,
                    "alias": alias,
                    "inputs": inputs,
                }
            )
        return normalized

    def decode_path(self, path, *, base_dir=None) -> None:
        return super().decode_path(path, base_dir=base_dir)

    def _resolve_command_text(self, command: str, scopes: tuple) -> str:
        text = str(command)
        for scope in scopes:
            if isinstance(scope, dict):
                text = self.resolve_placeholders_config(scope, text)
        text = self.resolve_placeholders_summary(text)
        unresolved = sorted(set(re.findall(r"\$\{[^}]+\}", text)))
        if unresolved:
            missing = ", ".join(unresolved)
            self.logger.error(
                f"Unresolved placeholder(s) in command: {missing}. Command: {command}"
            )
            sys.exit(2)
        return text

    def _next_cwd_from_command(self, command: str, current_cwd: str) -> str:
        line = str(command).strip()
        if not line:
            return current_cwd
        try:
            parts = shlex.split(line, posix=True)
        except ValueError:
            return current_cwd
        if not parts or parts[0] != "cd":
            return current_cwd
        # Keep a strict contract: only standalone `cd <path>` updates inherited cwd.
        if len(parts) != 2:
            return current_cwd
        target = parts[1]
        return self.decode_path(target, base_dir=current_cwd)

    def resolve_placeholders_config(self, data, text, separator=':'):
        pattern = re.compile(r'\$\{([^}]+)\}')
        def find_value(d, keys):
            if keys and keys[0] in d:
                return find_value(d[keys[0]], keys[1:]) if len(keys) > 1 else d[keys[0]]
            return None
        for match in pattern.finditer(text):
            key_path = match.group(1)
            keys = key_path.split(separator)
            value = find_value(data, keys)
            if value is not None:
                text = text.replace(match.group(0), str(value))
        return text
    
    def resolve_placeholders_summary(self, text):
        pattern = re.compile(r'\@\{([^}]+)\}')
        def find_value(d, keys):
            if keys and keys[0] in d:
                return find_value(d[keys[0]], keys[1:]) if len(keys) > 1 else d[keys[0]]
            return None

        for match in pattern.finditer(text):
            key_path = match.group(1)
            keys = key_path.split(":")
            value = find_value(self.summary, keys)
            if value is not None:
                text = text.replace(match.group(0), str(value).strip())
        return text





class ConfigValidator(Base):
    def __init__(self) -> None:
        super().__init__()
        self.logger = None
        self.config = None
        self.schema = None
        self.passcheck = False

    def set_config(self, config):
        self.config = config

    def set_schema(self, schema):
        with open(schema, 'r') as file:
            self.schema = json.load(file)

    def _build_local_schema_registry(self, schema_block):
        registry = Registry()
        for block_name, block_schema in schema_block.items():
            if not isinstance(block_schema, dict):
                continue
            ref_uri = block_schema.get("$ref")
            if not isinstance(ref_uri, str) or not ref_uri:
                continue

            parsed = urlparse(ref_uri)
            if parsed.scheme and parsed.scheme != "file":
                continue

            if parsed.scheme == "file":
                ref_path = Path(unquote(parsed.path))
            else:
                ref_path = Path(ref_uri)

            if not ref_path.is_file():
                raise FileNotFoundError(
                    f"schemaBlock.{block_name} points to missing file: {ref_path}"
                )

            with open(ref_path, "r", encoding="utf-8") as file:
                ref_schema = json.load(file)
            registry = registry.with_resource(
                ref_uri,
                Resource.from_contents(ref_schema, default_specification=DRAFT7),
            )
        return registry

    def validate_yaml(self):
        if self.config is None or self.schema is None:
            self.logger.error("Error: Config or schema not set.")
            sys.exit(2)
        try:
            schema_block = self.schema.get('schemaBlock', {})
            if not isinstance(schema_block, dict):
                schema_block = {}
            # Only patch refs that are declared by the active sampler schema.
            for kk, vv in self.schemablock.items():
                if kk in schema_block and isinstance(schema_block.get(kk), dict):
                    schema_block[kk]['$ref'] = vv

            validator_cls = validator_for(self.schema)
            validator_cls.check_schema(self.schema)
            registry = self._build_local_schema_registry(schema_block)
            validator = validator_cls(self.schema, registry=registry)
            validator.validate(self.config)
            self.logger.warning("Validation successful. The input YAML file meets the schema requirement.")
            self.passcheck = True
        except ValidationError as e:
            self.logger.error(f"Validation error: {e.message},\n the problematic module data: {e.instance}")
            self.passcheck = False
            sys.exit(2)
        except Exception as e:
            self.logger.error(f"Validation failed unexpectedly: {e}")
            self.passcheck = False
            sys.exit(2)


class ConfigTemplateGenerator():
    def __init__(self) -> None:
        self.schema = None
        self.output_path = None

    def set_schema(self, schema):
        self.schema = schema

    def set_output_path(self, output_path):
        self.output_path = output_path

    def generate_template(self):
        if self.schema is None or self.output_path is None:
            print("Error: Schema or output path not set.")
            sys.exit(2)
        
        template = {}
        for section, params in self.schema.items():
            template[section] = {param: '' for param in params}
        
        try:
            with open(self.output_path, 'w') as file:
                yaml.dump(template, file, default_flow_style=False)
            print(f"Template file is generated: {self.output_path}")
        except Exception as e:
            print(f"An error happened in generating template file: {e}")
            sys.exit(2)
