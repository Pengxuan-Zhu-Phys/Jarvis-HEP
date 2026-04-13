#!/usr/bin/env python3
from __future__ import annotations

from copy import deepcopy
import json
import os
import sys
import threading
import time 
import argparse
from jarvishep.config import ConfigLoader
import logging
from jarvishep.base import Base
from jarvishep.distributor import Distributor
from jarvishep.library import Library
from jarvishep.workflow import Workflow
from jarvishep.factory import WorkerFactory 
from jarvishep.sample import Sample
from jarvishep.moduleManager import ModuleManager
import asyncio 
from jarvishep.log_kv import format_two_column_log
from jarvishep.io_manager import IOManager
from loguru import logger
import setproctitle 
logger.remove()
from jarvishep.plot import JarvisPLOT as PlotterClass
from jarvishep.versioning import render_logo_with_version
from jarvishep.Sampling.Source.MCMC.runtime_checkpoint import StateSaver, stable_json_hash, validate_state_payload
# from monitor import Monitor

JARVIS_HEP_LOG_DOMAIN = "jarvis_hep"


class Core(Base):
    def __init__(self) -> None:
        # print("Testing the first logging line")
        super().__init__()
        self.argparser: Any             = None
        self.info: Dict[str, Any]       = {}
        self.yaml: ConfigLoader         = ConfigLoader()
        self.sampler: Any               = None
        self.libraries: Any             = None
        self.factory: Any               = None 
        self.likelihood: Any            = None 
        self.logger                     = None
        self.__state_saver              = None
        self.module_manager             = None
        self._funcs                     = {}
        self.tasks                      = []
        self.scan_mode                  = True
        self.plotter                    = PlotterClass()
        self.mode                       = None
        self.subprocess_scheduler       = None
        self.io_manager                 = None
        self.run_summary_collector      = None
        self.run_summary_renderer       = None
        self._resume_checkpoint_payload = None
        self._resume_run_spec           = None
        self._resume_factory_blueprint  = None
        self._resume_checkpoint_policy  = "auto"
        # self.monitor                    = Monitor()
        # self.monitor.start()

    def _render_logo_banner(self) -> str:
        return render_logo_with_version(self.path["logo"])

    def init_argparser(self) -> None:
        self.argparser = argparse.ArgumentParser(description="Jarvis Program Help Center", formatter_class=argparse.RawTextHelpFormatter)
        self.info['args'] = load_args_config(self.path['args_info'])
        
        for pos_arg in self.info['args'].get('positionals', []):
            heltp = pos_arg['help'].replace("$n", "\n")
            kwargs = {'help': heltp}
            # Keep `file` optional at argparse level; enforce requirements after parsing by mode.
            if pos_arg['name'] == "file":
                kwargs['nargs'] = '?'
            self.argparser.add_argument(pos_arg['name'], **kwargs)
        for opt in self.info['args'].get('options', []):
            kwargs = {
                'help': opt['help'],
                'action': opt.get('action', 'store'),
                'dest': opt.get('dest')
            }
            if 'default' in opt:
                kwargs['default'] = opt['default']
            if 'nargs' in opt:
                kwargs['nargs'] = opt['nargs']
            if 'const' in opt:
                kwargs['const'] = opt['const']
            if 'choices' in opt:
                kwargs['choices'] = opt['choices']
            if "type" in opt:
                if opt['type'] == 'int':
                    kwargs['type'] = int
                elif opt['type'] == 'float':
                    kwargs['type'] = float
                else:
                    kwargs['type'] = str
            if 'short' in opt and 'long' in opt:
                self.argparser.add_argument(
                    opt['short'],
                    opt['long'],
                    **kwargs
                )
            elif 'long' in opt:
                self.argparser.add_argument(
                    opt['long'], **kwargs
                )
        self.check_init_args()
        if getattr(self.args, "version", False):
            self.scan_mode = False
            self.mode = "VERSION"
            return
        if self.args.cvtDB:
            self.scan_mode = False
            self.mode      = "CDB"      # CDB means Convert DataBase from hdf5 into csv
        if self.args.plot:
            self.scan_mode = False
            self.mode      = "PLOT"
        if self.args.OPC: 
            self.scan_mode = True
            self.args.debug = True
            self.mode      = "1PC"      # 1PC means one point check mode 
        if self.args.monitor: 
            self.scan_mode = False
            self.mode       = "Monitor"
        
    def init_project(self) -> None: 
        import yaml 
        with open(os.path.abspath(self.args.file), 'r') as file:
            config = yaml.safe_load(file)
        self.info["project_name"] = os.path.splitext(os.path.basename(self.args.file))[0]
        self.info['config_file'] = os.path.abspath(self.args.file)

        # If YAML has a Scan section, use it; otherwise build a lightweight project layout for plotting-only YAML
        if 'Scan' in config and not self.args.plot:
            # Standard scan project layout
            self.info['scan_name'] = config['Scan']['name']
            self.info['proctitle'] = "Jarvis-HEP@{}".format(self.info['scan_name'])
            task_result_dir = os.path.join(config['Scan']['save_dir'], self.info['scan_name'])
            task_result_dir = self.decode_path(task_result_dir)
        elif 'Scan' in config and self.args.plot:
            # Plotting with a full scan YAML: reuse scan directories
            self.info['scan_name'] = config['Scan']['name']
            self.info['proctitle'] = "Jarvis-HEP@{}".format(self.info['scan_name'])
            task_result_dir = os.path.join(config['Scan']['save_dir'], self.info['scan_name'])
            task_result_dir = self.decode_path(task_result_dir)
        else:
            # Plotting-only YAML (no 'Scan' section): create a minimal project under the YAML's directory
            self.info['scan_name'] = self.info["project_name"]
            self.info['proctitle'] = "Jarvis-HEP@{}".format(self.info['scan_name'])
            yaml_dir = os.path.dirname(self.info['config_file'])
            # Use a 'PLOT' folder next to the YAML to hold outputs
            task_result_dir = os.path.join(yaml_dir, "PLOT", self.info['scan_name'])

        logs_dir = self.decode_path(os.path.join("&J", "logs", self.info["scan_name"]))
        images_dir = self.decode_path(os.path.join("&J", "images", self.info["scan_name"]))

        # Common paths (exist regardless of mode)
        self.info['logs_dir'] = logs_dir
        self.info['images_dir'] = images_dir
        self.info['jarvis_log'] = os.path.join(logs_dir, f"{self.info['scan_name']}.log")
        self.info['pickle_file'] = os.path.join(task_result_dir, f"{self.info['project_name']}.pkl")
        self.info['flowchart_path'] = os.path.join(images_dir, "flowchart.png")
        self.info['flowchart_semantic_path'] = os.path.join(images_dir, "flowchart.json")
        # Sampling method may be absent for plotting-only YAML
        sampling_method = None
        try:
            sampling_method = config['Sampling']['Method']
        except Exception:
            sampling_method = "Sampler"
        self.info['sampler_log'] = os.path.join(logs_dir, f"{sampling_method}.log")
        self.info['factory_log'] = os.path.join(logs_dir, "Factory.log")

        directory_cfg = {}
        scan_cfg = config.get("Scan", {})
        if isinstance(scan_cfg, dict):
            scan_dir_cfg = scan_cfg.get("sample_directory", {})
            if isinstance(scan_dir_cfg, dict):
                directory_cfg.update(scan_dir_cfg)
        if not directory_cfg:
            legacy_cfg = config.get("Directory_Setting", {})
            if isinstance(legacy_cfg, dict):
                directory_cfg.update(legacy_cfg)
        archive_samples = bool(directory_cfg.get("archive_samples", True))

        self.info['sample'] = {
            "task_result_dir": self.decode_path(task_result_dir),
            "sample_dirs": os.path.join(task_result_dir, "SAMPLE"),
            "jarvis_log":  self.info['jarvis_log'],
            "archive_samples": archive_samples,
        }
        self.info["db"] = {
            "path":  os.path.join(task_result_dir, "DATABASE", "samples.hdf5"),
            "info":  os.path.join(task_result_dir, "DATABASE", "running.json")
        }
        self.info['proc'] = {
            "path": os.path.join(task_result_dir, "DATABASE", ".pid.txt")
        }

        # Ensure directories exist (create minimal layout for plotting-only YAML)
        os.makedirs(task_result_dir, exist_ok=True)
        os.makedirs(os.path.join(task_result_dir, "SAMPLE"), exist_ok=True)
        os.makedirs(os.path.join(task_result_dir, "DATABASE"), exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)
        os.makedirs(images_dir, exist_ok=True)

        # Plot folder and config reference
        if self.args.plot:
            # Keep plot configs/assets at project scope instead of run-output scope.
            plot_dir = images_dir
            os.makedirs(plot_dir, exist_ok=True)
            if 'Scan' in config:
                # With full scan YAML, emit plotting config under project images root.
                self.info['plot'] = {
                    "save_path": plot_dir,
                    "config":    os.path.join(plot_dir, f"{self.info['scan_name']}.yaml")
                }
            else:
                # Plotting-only YAML: treat the provided YAML as plotting config
                self.info['plot'] = {
                    "save_path": plot_dir,
                    "config":    self.info['config_file']
                }

    def init_logger(self) -> None:
        logger.configure(extra={"_log_domain": JARVIS_HEP_LOG_DOMAIN})

        def global_log_filter(record):
            extra = record.get("extra", {})
            return extra.get("Jarvis", False) and (
                extra.get("_log_domain", "") == JARVIS_HEP_LOG_DOMAIN
            )

        def stream_filter(record):
            extra = record.get("extra", {})
            return extra.get("to_console", False) and (
                extra.get("_log_domain", "") == JARVIS_HEP_LOG_DOMAIN
            )

        def custom_format(record):
            module = record["extra"].get("module", "No module")
            if "raw" in record["extra"]:
                return "{message}\n"
            elif module == "Jarvis-HEP.hdf5-Writter":
                return f"\nϠ <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{{message}}</level> "
            elif "Sample@" in module:
                return f"\n·•· <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{{message}}</level>"
            else:
                return f"\n·•· <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{{message}}</level> "

        jarvislog = self.info["jarvis_log"]

        logger.add(jarvislog, 
                   rotation="5 MB", 
                   encoding='utf-8', 
                   format=custom_format, 
                   enqueue=True, 
                   level="WARNING",
                   filter=global_log_filter
                   )
        
        logger.add(
                sys.stdout, 
                filter=stream_filter, 
                format=custom_format,
                colorize=True,
                enqueue=True,
                level="DEBUG" if self.args.debug else "WARNING",
        )

        self.logger = logger.bind(
            module="Jarvis-HEP",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        self.logger.warning(f"\n{self._render_logo_banner()}")        
        self.logger.warning("Jarvis-HEP logging system initialized successful!")
        if self.args.debug:
            self.logger.info(f"Jarvis-HEP write into main log file -> {jarvislog}")
            self.logger.info("Jarvis-HEP in debug mode currently!")
        if self.args.plot:
            plogger = logger.bind(
                module="Jarvis-PLOT",
                to_console=True,
                Jarvis=False,
                _log_domain=JARVIS_HEP_LOG_DOMAIN,
            )
            self.plotter.logger = plogger

    def init_configparser(self) -> None: 
        self.yaml.logger = logger.bind(
            module="Jarvis-HEP.ConfigParser",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        from copy import deepcopy
        self.yaml.path = deepcopy(self.path)
        self.yaml.load_config(os.path.abspath(self.args.file))
        self.yaml.check_dependency_installed()
        self.sampler = Distributor.set_method(self.yaml.get_sampling_method()) 
        self.yaml.set_schema(self.sampler.schema)
        self.yaml.validate_config()

    def init_configparser_light(self) -> None:
        """Lightweight config load for plot/convert/monitor modes: load YAML only, skip scan-only setup."""
        self.yaml.logger = logger.bind(
            module="Jarvis-HEP.ConfigParser",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        from copy import deepcopy
        self.yaml.path = deepcopy(self.path)
        # Only load YAML; do not install dependencies, set sampler/schema, or validate scan config
        self.yaml.load_config(os.path.abspath(self.args.file))

    def init_utils(self) -> None: 
        if "Utils" in self.yaml.config:
            if "interpolations_1D" in self.yaml.config['Utils']:
                from jarvishep.utils import get_interpolate_1D_function_from_config
                for item in self.yaml.config['Utils']['interpolations_1D']:
                    if "file" in item:
                        item['file'] = self.decode_path(item['file'])
                    func = get_interpolate_1D_function_from_config(item)
                    if func is not None: 
                        self._funcs[item['name']] = func
                        self.logger.info(f"Succefully resolve the interpolate function -> {item['name']}")
                    else: 
                        self.logger.info(f"Illegal setting for interpolate function -> {item['name']}")
        self.init_operas_functions()

    def _build_operas_func_wrapper(self, registry, operator_name, input_names):
        input_names = input_names if isinstance(input_names, list) else []

        def _wrapped(*args, **kwargs):
            call_kwargs = dict(kwargs)
            if "observables" not in call_kwargs:
                if len(args) == 1 and isinstance(args[0], dict):
                    call_kwargs["observables"] = args[0]
                elif args:
                    if input_names and len(input_names) == len(args):
                        call_kwargs["observables"] = {k: v for k, v in zip(input_names, args)}
                    else:
                        raise ValueError(
                            f"Operas function '{operator_name}' received positional args but no valid input mapping."
                        )
            observables_payload = call_kwargs.get("observables")
            if isinstance(observables_payload, dict):
                for key, value in observables_payload.items():
                    call_kwargs.setdefault(str(key), value)
            return registry.call(operator_name, logger=self.logger, **call_kwargs)

        return _wrapped

    def init_operas_functions(self) -> None:
        operas_cfg = (self.yaml.config.get("Operas", {}) or {})
        if not operas_cfg:
            return

        has_operas = bool(operas_cfg.get("Modules")) or bool(operas_cfg.get("Functions"))
        if not has_operas:
            return

        try:
            from jarvis_operas import get_global_operas_registry
            from importlib.metadata import PackageNotFoundError, version as dist_version
        except Exception as exc:
            self.logger.error(f"Jarvis-Operas is required by Operas config but unavailable: {exc}")
            sys.exit(2)

        try:
            operas_version = dist_version("Jarvis-Operas")
        except PackageNotFoundError:
            operas_version = "unknown(local-source)"
        self.logger.warning(f"Loading Jarvis-Operas version -> {operas_version}")

        registry = get_global_operas_registry()
        for item in self.yaml.get_operas_function_whitelist():
            alias = item["alias"]
            if not alias.isidentifier():
                self.logger.error(f"Illegal Operas alias '{alias}'. Alias must be a valid identifier.")
                sys.exit(2)
            try:
                registry.resolve_name(item["name"])
            except Exception as exc:
                self.logger.error(f"Operas function '{item['name']}' is not registered: {exc}")
                sys.exit(2)
            if alias in self._funcs:
                self.logger.warning(f"Operas alias '{alias}' overrides existing function with the same name.")
            wrapper = self._build_operas_func_wrapper(registry, item["name"], item.get("inputs", []))
            self._funcs[alias] = wrapper
            self.logger.info(f"Register Operas function '{item['name']}' as alias '{alias}'")

    def init_StateSaver(self) -> None:
        slogger = logger.bind(
            module="Jarvis-HEP.Checkpoint",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        checkpoint_root = self._checkpoint_root_for_sampler()
        self.info["checkpoint_root"] = checkpoint_root
        self.info["checkpoint_file"] = os.path.join(checkpoint_root, "state.pkl")

        if not getattr(self.sampler, "supports_runtime_checkpointing", lambda: False)():
            slogger.info(
                f"Runtime checkpointing skipped for sampler -> {getattr(self.sampler, 'method', 'unknown')}"
            )
            return

        os.makedirs(checkpoint_root, exist_ok=True)
        slogger.warning(f"Configuring runtime checkpointing -> {checkpoint_root}")
        self.sampler.configure_runtime_checkpointing(
            checkpoint_root,
            interval_seconds=30.0,
            auto_resume=True,
            logger=slogger,
        )
        if self._resume_run_spec is not None:
            self.sampler.set_runtime_checkpoint_context(run_spec=self._resume_run_spec)
        if self._resume_factory_blueprint is not None:
            self.sampler.set_runtime_checkpoint_context(factory_blueprint=self._resume_factory_blueprint)
        resumed = self.sampler.restore_runtime_checkpoint_if_available(self._resume_checkpoint_payload)
        if resumed:
            slogger.warning("Breakpoint resume state restored into sampler runtime")

    def _checkpoint_root_for_sampler(self) -> str:
        scan_name = self.info.get("scan_name", self.info.get("project_name", "scan"))
        sampler_name = getattr(self.sampler, "method", None)
        if not sampler_name:
            try:
                sampler_name = self.yaml.get_sampling_method()
            except Exception:
                sampler_name = "sampler"
        return os.path.join(self.path.get("task_root", os.getcwd()), "checkpoints", str(scan_name), str(sampler_name))

    def _checkpoint_file_for_sampler(self) -> str:
        return os.path.join(self._checkpoint_root_for_sampler(), "state.pkl")

    def _prompt_resume_from_checkpoint(self, checkpoint_file: str, timeout_seconds: float = 30.0) -> bool:
        prompt = "Detected checkpoint file. Re-run from scratch? [y/N] (default: resume in 30s): "
        clogger = getattr(self, "logger", None)
        if not getattr(sys.stdin, "isatty", lambda: False)():
            if clogger is not None:
                clogger.warning(
                    f"Detected checkpoint file -> {checkpoint_file}, but no interactive terminal is available; resuming from checkpoint."
                )
            return True

        response = {"text": None, "error": None}

        def _read_response():
            try:
                response["text"] = sys.stdin.readline()
            except Exception as exc:
                response["error"] = exc

        print(prompt, end="", flush=True)
        reader = threading.Thread(target=_read_response, daemon=True)
        reader.start()
        reader.join(timeout_seconds)
        if reader.is_alive():
            if clogger is not None:
                clogger.warning(
                    f"Checkpoint prompt timed out after {int(timeout_seconds)}s; resuming from checkpoint."
                )
            return True

        if response["error"] is not None:
            if clogger is not None:
                clogger.warning(
                    f"Checkpoint prompt failed ({response['error']}); resuming from checkpoint."
                )
            return True

        answer = str(response["text"] or "").strip().lower()
        if answer in {"y", "yes"}:
            return False
        return True

    def _discard_existing_checkpoint(self, checkpoint_file: str) -> None:
        clogger = getattr(self, "logger", None)
        try:
            os.remove(checkpoint_file)
            if clogger is not None:
                clogger.warning(f"Discarded existing checkpoint file -> {checkpoint_file}")
        except FileNotFoundError:
            return
        except Exception as exc:
            if clogger is not None:
                clogger.warning(
                    f"Failed to remove existing checkpoint file -> {checkpoint_file}: {exc}"
                )

    def _build_run_spec(self) -> dict:
        raw_yaml_text = None
        config_file = self.info.get("config_file")
        if isinstance(config_file, str) and os.path.exists(config_file):
            try:
                with open(config_file, "r", encoding="utf-8") as handle:
                    raw_yaml_text = handle.read()
            except Exception:
                raw_yaml_text = None

        return {
            "raw_yaml_text": raw_yaml_text,
            "normalized_config": deepcopy(self.yaml.config),
            "scan_name": self.info.get("scan_name"),
            "task_root": self.path.get("task_root"),
            "task_result_dir": self.info.get("sample", {}).get("task_result_dir"),
            "logs_dir": self.info.get("logs_dir"),
            "images_dir": self.info.get("images_dir"),
            "worker_parallel": int(getattr(self.factory, "_max_workers", self.yaml.get_worker_parallel())),
            "sampler_method": getattr(self.sampler, "method", None),
            "workflow": deepcopy(getattr(self.workflow, "workflow", {})),
            "workflow_layers": deepcopy(getattr(self.workflow, "calc_layer", {})),
        }

    def _build_factory_blueprint(self) -> dict:
        module_pools = {}
        if getattr(self.module_manager, "module_pools", None):
            for name, pool in self.module_manager.module_pools.items():
                if hasattr(pool, "export_blueprint"):
                    module_pools[str(name)] = pool.export_blueprint()
        return {
            "max_workers": int(getattr(self.factory, "_max_workers", 0) or 0),
            "workflow": deepcopy(getattr(self.workflow, "workflow", {})),
            "calc_layer": deepcopy(getattr(self.workflow, "calc_layer", {})),
            "module_pools": module_pools,
        }

    def _apply_resume_runtime_spec(self, run_spec: dict) -> None:
        if not isinstance(run_spec, dict):
            return
        normalized_config = run_spec.get("normalized_config")
        if isinstance(normalized_config, dict):
            self.yaml.config = deepcopy(normalized_config)

        scan_name = run_spec.get("scan_name")
        if scan_name:
            self.info["scan_name"] = str(scan_name)

        task_root = run_spec.get("task_root")
        if task_root:
            self.path["task_root"] = os.path.abspath(str(task_root))
            self.path["jpath"] = self.path["task_root"]

        task_result_dir = run_spec.get("task_result_dir")
        if isinstance(task_result_dir, str) and task_result_dir:
            self.info.setdefault("sample", {})
            self.info["sample"]["task_result_dir"] = task_result_dir
            self.info["sample"]["sample_dirs"] = os.path.join(task_result_dir, "SAMPLE")
            self.info.setdefault("db", {})
            self.info["db"]["path"] = os.path.join(task_result_dir, "DATABASE", "samples.hdf5")
            self.info["db"]["info"] = os.path.join(task_result_dir, "DATABASE", "running.json")
            self.info.setdefault("proc", {})
            self.info["proc"]["path"] = os.path.join(task_result_dir, "DATABASE", ".pid.txt")

        logs_dir = run_spec.get("logs_dir")
        if isinstance(logs_dir, str) and logs_dir:
            self.info["logs_dir"] = logs_dir
            self.info["jarvis_log"] = os.path.join(logs_dir, f"{self.info.get('scan_name', 'scan')}.log")
            self.info["factory_log"] = os.path.join(logs_dir, "Factory.log")

        images_dir = run_spec.get("images_dir")
        if isinstance(images_dir, str) and images_dir:
            self.info["images_dir"] = images_dir

        self._resume_run_spec = deepcopy(run_spec)

    def _preload_resume_checkpoint(self) -> None:
        slogger = logger.bind(
            module="Jarvis-HEP.Checkpoint",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        supported_methods = {
            "MCMC",
            "ToyMCMC",
            "PTMCMC",
            "AMMCMC",
            "RobustAM",
            "DRAM",
            "DEMCMC",
            "DREAM",
            "DREAMLite",
            "EnsembleMCMC",
            "PTEnsemble",
            "SliceMCMC",
            "ESS",
            "MALA",
            "HMC",
            "NUTS",
            "Grid",
            "Random",
            "Bridson",
            "CSV",
            "DNN",
            "Diver",
            "Dynesty",
            "MultiNest",
            "RLTPMCMC",
        }
        try:
            method_name = str(self.yaml.get_sampling_method())
        except Exception:
            method_name = ""
        if method_name not in supported_methods:
            self._resume_checkpoint_payload = None
            self._resume_run_spec = None
            self._resume_factory_blueprint = None
            self._resume_checkpoint_policy = "auto"
            return

        checkpoint_root = self._checkpoint_root_for_sampler()
        current_checkpoint = os.path.join(checkpoint_root, "state.pkl")
        if self._resume_checkpoint_policy == "fresh":
            self._resume_checkpoint_payload = None
            self._resume_run_spec = None
            self._resume_factory_blueprint = None
            return

        if not os.path.exists(current_checkpoint):
            if getattr(self.args, "resume", False):
                slogger.warning(
                    f"Checkpoint resume requested but no checkpoint was found -> {current_checkpoint}. Starting a fresh run."
                )
            self._resume_checkpoint_payload = None
            self._resume_run_spec = None
            self._resume_factory_blueprint = None
            return

        if not getattr(self.args, "resume", False):
            resume_from_checkpoint = self._prompt_resume_from_checkpoint(current_checkpoint, timeout_seconds=30.0)
            if not resume_from_checkpoint:
                self._resume_checkpoint_policy = "fresh"
                self._resume_checkpoint_payload = None
                self._resume_run_spec = None
                self._resume_factory_blueprint = None
                self._discard_existing_checkpoint(current_checkpoint)
                slogger.warning("Starting a fresh run from user confirmation; existing checkpoint was discarded.")
                return
            self._resume_checkpoint_policy = "resume"
        else:
            self._resume_checkpoint_policy = "resume"

        payload = None
        saver = StateSaver(current_checkpoint, logger=slogger)
        payload = saver.load(default=None)
        if not isinstance(payload, dict):
            self._resume_checkpoint_payload = None
            self._resume_run_spec = None
            self._resume_factory_blueprint = None
            return

        ok, reason = validate_state_payload(payload)
        if not ok:
            slogger.warning(f"StateSaver checkpoint rejected -> {reason}")
            self._resume_checkpoint_payload = None
            self._resume_run_spec = None
            self._resume_factory_blueprint = None
            return

        integrity = dict(payload.get("integrity", {}))
        saved_hash = integrity.get("config_hash")
        current_hash = stable_json_hash(self.yaml.config)
        if saved_hash and saved_hash != current_hash:
            slogger.warning(
                "Resume checkpoint config drift detected; using frozen checkpoint run_spec instead of current YAML."
            )

        self._resume_checkpoint_payload = payload
        run_spec = payload.get("run_spec", {})
        if isinstance(run_spec, dict):
            self._apply_resume_runtime_spec(run_spec)
            self._resume_factory_blueprint = deepcopy(payload.get("factory_blueprint", {}))

    def init_sampler(self) -> None:
        logger_name = f"Jarvis-HEP.{self.sampler.method}"
        def filte_func(record):
            return logger_name in record['extra']['module']
        def custom_format(record):
            module = record["extra"].get("module", "No module")
            if "raw" in record["extra"]:
                return "{message}"
            return f"\n·•· <red>{module}</red> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{{message}}</level>"
    
        # logger = self.logger.create_dynamic_logger(self.sampler.method)
        slogger = logger.bind(
            module=logger_name,
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        slogger.add(self.info['sampler_log'], format=custom_format, level="DEBUG", rotation=None, retention=None, filter=filte_func )
        self.sampler.info['logfile'] = self.info['sampler_log']
        self.sampler.info['sample'] = deepcopy(self.info['sample'])
        self.sampler.set_logger(slogger)
        self.sampler.set_config(self.yaml.config)
        if self._resume_checkpoint_payload is not None and hasattr(self.sampler, "set_runtime_checkpoint_resume_hint"):
            self.sampler.set_runtime_checkpoint_resume_hint(True)
        self.sampler.initialize()
        self.yaml.vars = self.sampler.vars 

    def init_librarys(self) -> None:
        if hasattr(self.yaml, "SupportLibrary"):
            self.libraries = Library()
            self.libraries._skip_library = self.args.skiplibrary
            # logger = self.logger.create_dynamic_logger("Library", logging.INFO)
            slogger = logger.bind(
                module="Jarvis-HEP.Library",
                to_console=True,
                Jarvis=True,
                _log_domain=JARVIS_HEP_LOG_DOMAIN,
            )
            self.libraries.set_logger(slogger)
            self.libraries.set_config(self.yaml.config)
    
            self.libraries.display_installation_summary()
            for module in self.libraries.modules.values():
                module.install()
            
    def init_workflow(self) -> None: 
        self.workflow = Workflow()
        modules = self.yaml.get_modules()
        self.workflow.set_modules(modules)
        self.workflow.resolve_dependencies()
        try:
            self.logger.warning(
                f"Export workflow semantic graph into {self.info['flowchart_semantic_path']}"
            )
            asyncio.run(
                self.workflow.export_flowchart_semantics(
                    save_path=self.info['flowchart_semantic_path'],
                    workflow_name=self.info.get("project_name", "Workflow"),
                )
            )
        except Exception as exc:
            self.logger.warning(
                f"Skipping semantic flowchart export after failure: {exc}"
            )
        if not self.args.skipFC:
            self.logger.warning(f"Draw workflow chart into {self.info['flowchart_path']}")
            asyncio.run(
                self.workflow.draw_flowchart(save_path=self.info['flowchart_path'])
            )
        self.workflow.get_workflow_dict()

    def init_WorkerFactory(self) -> None: 
        self.factory = WorkerFactory()
        self.module_manager = ModuleManager()
        max_workers = self.yaml.get_worker_parallel()
        self.io_manager = IOManager(max_workers=max_workers)
        self.factory.configure(module_manager=self.module_manager,
            max_workers=max_workers
            )
        self.sampler.set_max_workers(max_workers)
        self.module_manager.set_io_manager(self.io_manager)

        if hasattr(self.yaml, "get_subprocess_runtime_options"):
            runtime_opts = self.yaml.get_subprocess_runtime_options(worker_parallel=max_workers)
        else:
            fallback_concurrency = max(1, int(max_workers))
            runtime_opts = {
                "max_concurrency": fallback_concurrency,
                "max_pending": max(128, fallback_concurrency * 16),
                "per_task_timeout_sec": None,
                "progress_interval_sec": 5.0,
                "log_policy": "logger",
                "diagnostics_enabled": False,
                "diagnostics_interval_sec": 10.0,
                "terminate_grace_sec": 5.0,
            }

        from jarvishep.async_subprocess import AsyncSubprocessScheduler, SubprocessRuntimeConfig

        task_result_dir = self.info.get("sample", {}).get("task_result_dir", "/tmp")
        runtime_dir = os.path.join(task_result_dir, "RUNTIME")
        status_path = os.path.join(runtime_dir, "status.jsonl")

        splogger = logger.bind(
            module="Jarvis-HEP.Subprocess",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        self.subprocess_scheduler = AsyncSubprocessScheduler(
            config=SubprocessRuntimeConfig(**runtime_opts),
            logger=splogger,
            status_path=status_path,
        )
        self.factory.set_subprocess_scheduler(self.subprocess_scheduler)
        self.module_manager.set_subprocess_scheduler(self.subprocess_scheduler)
        splogger.warning(
            format_two_column_log(
                "Subprocess scheduler configured",
                [
                    ("max_concurrency", runtime_opts.get("max_concurrency")),
                    ("max_pending", runtime_opts.get("max_pending")),
                    ("timeout", runtime_opts.get("per_task_timeout_sec")),
                    ("log_policy", runtime_opts.get("log_policy")),
                    ("diagnostics", runtime_opts.get("diagnostics_enabled")),
                ],
            )
        )

        factory_logger_name = "Jarvis-HEP.Factory"

        def _factory_filter(record):
            module = record["extra"].get("module", "")
            return module.startswith(factory_logger_name)

        def _factory_format(record):
            module = record["extra"].get("module", "No module")
            if "raw" in record["extra"]:
                return "{message}"
            return (
                f"\n·•· <magenta>{module}</magenta> \n\t-> "
                f"<green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - "
                "[<level>{level}</level>] >>> \n<level>{message}</level>"
            )

        flogger = logger.bind(
            module=factory_logger_name,
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        factory_log_path = self.info.get("factory_log")
        if not factory_log_path:
            logs_dir = self.info.get("logs_dir")
            if logs_dir:
                factory_log_path = os.path.join(logs_dir, "Factory.log")
            else:
                sample_cfg = self.info.get("sample", {}) if isinstance(self.info, dict) else {}
                task_result_dir = sample_cfg.get("task_result_dir", "/tmp")
                factory_log_path = os.path.join(task_result_dir, "Factory.log")
            self.info["factory_log"] = factory_log_path
        flogger.add(
            factory_log_path,
            format=_factory_format,
            level="DEBUG",
            rotation="10 MB",
            retention=None,
            filter=_factory_filter,
        )
        self.factory.set_logger(flogger)
        mlogger = logger.bind(
            module="Jarvis-HEP.Factory.Manager",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        self.module_manager.set_logger(mlogger)
        self.module_manager.set_max_worker(max_workers)
        self.module_manager.set_config(self.yaml.config)
        self.module_manager.set_funcs(self._funcs)
        self.module_manager.workflow = deepcopy(self.workflow.workflow)
        for kk, layer in self.workflow.calc_layer.items():
            if kk > 1: 
                for module in layer['module']:
                    self.module_manager.add_module_pool(self.workflow.modules[module])
        self.factory.info['sample'] = deepcopy(self.info['sample'])
        if self.sampler._with_nuisance: 
            self.module_manager.nuisance_loglikelihoods = self.sampler.nuisance_sampler.loglikelihoods
            self.module_manager.nuisance_passconditions = self.sampler.nuisance_sampler.passconditions
            self.module_manager._with_nuisance = True
        if self._resume_factory_blueprint is not None:
            self.module_manager.restore_factory_blueprint(self._resume_factory_blueprint)

    def shutdown_io_manager(self) -> None:
        manager = getattr(self, "io_manager", None)
        if manager is None:
            return
        manager.shutdown(wait=True, cancel_futures=True)
        self.io_manager = None

    def init_likelihood(self) -> None: 
        if self.yaml.config['Sampling'].get("LogLikelihood", False):
            self.module_manager.set_likelihood()
            self.module_manager.set_funcs(self._funcs)
            from copy import deepcopy
            self.sampler.set_likelihood(deepcopy(self.module_manager.loglikelihood))

    def init_database(self) -> None: 
        from jarvishep.hdf5writer import GlobalHDF5Writer
        self.module_manager._database = GlobalHDF5Writer(self.info['db'])
        self.module_manager.database.start()

    def init_run_summary(self) -> None:
        from jarvishep.monitoring.run_summary import (
            RunSummaryCollector,
            RunSummaryRenderer,
        )

        sample_cfg = self.info.get("sample", {}) if isinstance(self.info, dict) else {}
        task_result_dir = sample_cfg.get("task_result_dir")
        if not task_result_dir:
            return

        sampler_name = getattr(getattr(self, "sampler", None), "method", None)
        if not sampler_name:
            try:
                sampler_name = self.yaml.get_sampling_method()
            except Exception:
                sampler_name = None

        configured_workers = None
        try:
            configured_workers = self.yaml.get_worker_parallel()
        except Exception:
            configured_workers = None

        self.run_summary_collector = RunSummaryCollector(
            output_dir=task_result_dir,
            project_name=self.info.get("project_name"),
            sampler_name=sampler_name,
            configured_workers=configured_workers,
            run_label=self.info.get("scan_name") or self.info.get("project_name"),
        )
        self.run_summary_renderer = RunSummaryRenderer()

        if getattr(self, "factory", None) is not None:
            self.factory.set_run_summary_collector(self.run_summary_collector)
            self.run_summary_collector.attach_factory(self.factory)
        if getattr(self, "subprocess_scheduler", None) is not None:
            self.run_summary_collector.attach_scheduler(self.subprocess_scheduler)
        if getattr(self, "module_manager", None) is not None:
            for pool in self.module_manager.module_pools.values():
                if hasattr(pool, "set_run_summary_collector"):
                    pool.set_run_summary_collector(self.run_summary_collector)

    def _start_run_summary(self) -> None:
        if self.run_summary_collector is not None:
            self.run_summary_collector.start()

    def _emit_run_summary(self) -> None:
        if self.run_summary_collector is None or self.run_summary_renderer is None:
            return
        try:
            summary = self.run_summary_collector.finish()
            rendered = self.run_summary_renderer.render(summary)
            if getattr(self, "logger", None) is not None:
                self.logger.warning(rendered.rstrip("\n"))
            self.run_summary_renderer.write_outputs(
                summary,
                self.info["sample"]["task_result_dir"],
                rendered_text=rendered,
            )
        except Exception as exc:
            if getattr(self, "logger", None) is not None:
                self.logger.error(f"Run summary emission failed -> {exc}")

    def initialization(self) -> None:
        # Parse CLI and decide mode
        self.init_argparser()
        self.configure_runtime_context(
            config_path=(os.path.abspath(self.args.file) if getattr(self.args, "file", None) else None),
            cwd=os.getcwd(),
        )

        # Common initializations for all modes that need project context/logging
        def _init_common_project_and_logger():
            self.init_project()
            self.init_logger()


        # Mode switch: PLOT / CDB / Monitor / SCAN (default) / 1PC
        if self.mode == "VERSION" or getattr(self.args, "version", False):
            # Version mode: no project/config/runtime initialization required
            return

        elif self.mode == "PLOT" or getattr(self.args, 'plot', False):
            # Plot mode: static config emission only; keep optional runtime
            # integrations such as Jarvis-Operas out of this path.
            _init_common_project_and_logger()
            self.init_configparser_light()
            # main() will call self.plot() afterwards
            return
        
        elif self.mode == "CDB" or getattr(self.args, 'cvtDB', False):
            # Convert mode: set up paths/logging; optionally light-load config
            _init_common_project_and_logger()
            try:
                self.init_configparser_light()
            except Exception:
                pass
            # main() will call self.convert() afterwards
            return

        elif self.mode == "Monitor" or getattr(self.args, 'monitor', False):
            # Monitor mode: minimal init for logger + project paths
            _init_common_project_and_logger()
            # main() will call self.monitor() afterwards
            return

        else:
            # SCAN / 1PC modes (full pipeline)
            _init_common_project_and_logger()
            self.init_configparser()
            self.init_utils()
            self._preload_resume_checkpoint()
            
            setproctitle.setproctitle(self.info['proctitle'])
            self.logger.warning("Setting process title -> {}".format(self.info['proctitle']))
            self.save_pid()
            self.init_sampler()
            self.init_StateSaver()
            self.init_workflow()
            self.init_librarys()
            self.init_WorkerFactory()
            self.init_likelihood()
            self.init_database()
            self.init_run_summary()
            return

    def save_pid(self):
        pid = os.getpid()
        with open(self.info['proc']['path'], "w") as f1:
            f1.write("{}".format(pid))
        
    def run_sampling(self)->None:
        if self.mode == "1PC":
            self.test_assembly_line()
        else:
            self.run_until_finished()

    def test_assembly_line(self):
        self.logger.warning("Start testing assembly line")
        self._start_run_summary()
        try:
            for ii in range(10):
                param = next(self.sampler)
                sample = Sample(param)
                sample.set_config(deepcopy(self.info['sample']))
                self.logger.warning(f"Run test assembly line for sample -> {sample.info['uuid']}\n")
                try:
                    future = self.factory.submit_task(sample.info)
                    output = future.result()
                    self.module_manager.database.add_data(output)
                finally:
                    sample.close()
        except Exception as e:
            # 异常处理
            self.logger.error(f"An error occurred: {e}")

        finally:
            self.factory.shutdown()
            self.module_manager.database.stop()
            self.shutdown_io_manager()

            from time import time
            start = time()
            tot = 1000 * (time() - start)
            self.logger.info(f"{tot} millisecond -> All samples have been processed.")
            self._emit_run_summary()
            # self.monitor.stop()

    def run_until_finished(self):
        self.sampler.set_factory(factory = self.factory)
        if hasattr(self.sampler, "set_runtime_checkpoint_context"):
            self.sampler.set_runtime_checkpoint_context(
                run_spec_getter=self._build_run_spec,
                factory_blueprint_getter=self._build_factory_blueprint,
            )
        run_exc = None
        self._start_run_summary()
        try:
            self.sampler.run_nested()
        except Exception as exc:
            run_exc = exc
            raise
        finally:
            from time import time
            cleanup_errors = []

            def _cleanup_step(step_name, fn):
                try:
                    fn()
                except Exception as exc:
                    cleanup_errors.append((step_name, exc))
                    self.logger.error(f"Cleanup step failed -> {step_name} -> {exc}")

            start = time()
            _cleanup_step(
                "sampler.persist_runtime_checkpoint",
                lambda: self.sampler.persist_runtime_checkpoint(force=True, reason="core_cleanup"),
            )
            _cleanup_step("database.stop", self.factory.module_manager.database.stop)
            tot = 1000 * (time() - start)
            self.logger.info(f"{tot} millisecond -> All samples have been processed.")
            _cleanup_step("sampler.finalize", self.sampler.finalize)
            _cleanup_step(
                "sampler.shutdown_runtime_checkpointing",
                lambda: getattr(self.sampler, "shutdown_runtime_checkpointing", lambda: None)(),
            )
            _cleanup_step("sampler.combine_data", lambda: self.sampler.combine_data(self.info['db']['path']))
            _cleanup_step("factory.shutdown", self.factory.shutdown)
            _cleanup_step("sampler.finalize_sample_archive", self.sampler.finalize_sample_archive)
            _cleanup_step("io_manager.shutdown", self.shutdown_io_manager)
            self._emit_run_summary()

            if run_exc is None and cleanup_errors:
                first_step, first_exc = cleanup_errors[0]
                raise RuntimeError(f"Cleanup failed after sampling -> {first_step}") from first_exc
            # self.monitor.stop()

    def check_init_args(self) -> None:
        try:
            self.args = self.argparser.parse_args()
        except argparse.ArgumentError as e:
            print(str(e))
            self.argparser.print_help()
            sys.exit(2)

        if getattr(self.args, "version", False):
            return

        if not self.args.file:
            self.argparser.error("the following arguments are required: file")

        config_file = os.path.abspath(os.path.expanduser(str(self.args.file)))
        if not os.path.isfile(config_file):
            self.argparser.error(f"YAML file not found: {config_file}")
        self.args.file = config_file

    def convert(self) -> None: 
        if os.path.exists(self.info['db']['info']):
            with open(self.info['db']['info'], 'r') as f1:
                dbinfo = json.load(f1)
            out_csv = "".join([dbinfo['pathroot'], ".{}".format(dbinfo['activeNO']), ".csv"])

            if out_csv not in dbinfo['converted'] or not os.path.exists(out_csv):
                if os.path.exists(dbinfo['active path']):
                    self.logger.warning("Jarvis-HEP not find the -> {}.\nStarting converting from hdf5.".format(out_csv))
                    snapshot_path = dbinfo['active path'] + ".snap"
                    import shutil 
                    shutil.copy2(dbinfo['active path'], snapshot_path)
                    try: 
                        from jarvishep.utils import convert_hdf5_to_csv
                        convert_hdf5_to_csv(snapshot_path, out_csv)
                        self.logger.warning(f"Data has been successfully written to -> {out_csv}")
                    finally:
                        try:
                            os.remove(snapshot_path)
                            self.logger.warning(f"Snapshot -> {snapshot_path} has been successfully deleted.")
                        except Exception as e:
                            self.logger.error(f"Failed to delete snapshot {snapshot_path}: {str(e)}")

    def plot(self) -> None:
        self.plotter.logger.warning(f"Generate JarvisPLOT YAML for {self.info['project_name']}")
        self.plotter.get_plot_config_from_Jarvis(self.info, self.yaml)
        
    def monitor(self) -> None: 
        from jarvishep.monitor import JarvisMonitor
        import curses
        self.logger.warning("Start monitoring -> {}".format(self.info['proctitle']))
        if os.path.exists(self.info['proc']['path']):
            pid = None
            with open(self.info['proc']['path'], 'r') as f1:
                pid = f1.read().strip()
            self.monitor = JarvisMonitor(pid)
            asyncio.run(curses.wrapper(self.monitor.main))

def load_args_config(json_file):
    with open(json_file, 'r') as file:
        config = json.load(file)
    return config
