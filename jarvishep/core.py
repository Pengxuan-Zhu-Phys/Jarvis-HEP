#!/usr/bin/env python3
from __future__ import annotations

from copy import deepcopy
import json
import os
from subprocess import Popen, run
import sys
import time 
import argparse
from jarvishep.config import ConfigLoader
import logging
from jarvishep.base import Base
from jarvishep.distributor import Distributor
import dill 
from threading import Timer
from jarvishep.library import Library
from jarvishep.workflow import Workflow
from jarvishep.factory import WorkerFactory 
from jarvishep.sample import Sample
from jarvishep.moduleManager import ModuleManager
from pprint import pprint
import pandas as pd 
import concurrent.futures
import asyncio 
from jarvishep.project_scaffold import PROJECT_SUBDIRS, create_project_scaffold
from loguru import logger
import setproctitle 
logger.remove()
from jarvishep.plot import JarvisPLOT as PlotterClass
from jarvishep.versioning import render_logo_with_version
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
        self.async_loop                 = asyncio.get_event_loop()
        self.scan_mode                  = True
        self.plotter                    = PlotterClass()
        self.mode                       = None
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
        if self.args.mkproject:
            self.scan_mode = False
            self.mode      = "MKPROJECT"
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

        # Common paths (exist regardless of mode)
        self.info['jarvis_log'] = os.path.join(task_result_dir, "LOG", f"{self.info['scan_name']}.log")
        self.info['pickle_file'] = os.path.join(task_result_dir, f"{self.info['project_name']}.pkl")
        self.info['flowchart_path'] = os.path.join(task_result_dir, "flowchart.png")
        # Sampling method may be absent for plotting-only YAML
        sampling_method = None
        try:
            sampling_method = config['Sampling']['Method']
        except Exception:
            sampling_method = "Sampler"
        self.info['sampler_log'] = os.path.join(task_result_dir, "LOG", f"{sampling_method}.log")

        self.info['sample'] = {
            "task_result_dir": self.decode_path(task_result_dir),
            "sample_dirs": os.path.join(task_result_dir, "SAMPLE"),
            "jarvis_log":  self.info['jarvis_log']
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
        os.makedirs(os.path.join(task_result_dir, "LOG"), exist_ok=True)
        os.makedirs(os.path.join(task_result_dir, "DATABASE"), exist_ok=True)

        # Plot folder and config reference
        if self.args.plot:
            plot_dir = os.path.join(task_result_dir, "IMAGE")
            os.makedirs(plot_dir, exist_ok=True)
            if 'Scan' in config:
                # With full scan YAML, keep generated plotting config under IMAGE
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
        from datetime import datetime
        current_time = datetime.now().strftime("%Y-%m-%d[%H:%M:%S]")
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
                return f"\n <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{{message}}</level> "

        jarvislog = f"{self.info['jarvis_log']}@{current_time}"

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
        # logger = self.logger.create_dynamic_logger("StateSaver", logging.INFO)
        slogger = logger.bind(
            module="Jarvis-HEP.StateSaver",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
        )
        slogger.warning("Enabling breakpoint resume function ... ")
        self.__state_saver = self.__StateSaver(self, filename=self.info['pickle_file'] , logger=slogger, save_interval_seconds=60)

    def init_sampler(self) -> None:
        logger_name = f"Jarvis-HEP.{self.sampler.method}"
        def filte_func(record):
            return logger_name in record['extra']['module']
        def custom_format(record):
            module = record["extra"].get("module", "No module")
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
        self.factory.configure(module_manager=self.module_manager,
            max_workers=max_workers
            )
        self.sampler.set_max_workers(max_workers)
        flogger = logger.bind(
            module="Jarvis-HEP.Factory",
            to_console=True,
            Jarvis=True,
            _log_domain=JARVIS_HEP_LOG_DOMAIN,
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


        # Mode switch: MKPROJECT / PLOT / CDB / Monitor / SCAN (default) / 1PC
        if self.mode == "MKPROJECT" or getattr(self.args, 'mkproject', None):
            # Project scaffold mode: no YAML/project logger required
            return

        elif self.mode == "VERSION" or getattr(self.args, "version", False):
            # Version mode: no project/config/runtime initialization required
            return

        elif self.mode == "PLOT" or getattr(self.args, 'plot', False):
            # Plot mode: lightweight setup; no scan-only preprocessing
            _init_common_project_and_logger()
            self.init_configparser_light()
            self.init_utils()
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
            
            setproctitle.setproctitle(self.info['proctitle'])
            self.logger.warning("Setting process title -> {}".format(self.info['proctitle']))
            self.save_pid()
            self.init_StateSaver()

            self.init_sampler()
            self.init_workflow()
            self.init_librarys()
            self.init_WorkerFactory()
            self.init_likelihood()
            self.init_database()
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
        try:
            for ii in range(10):
                param = next(self.sampler)
                sample = Sample(param)
                sample.set_config(deepcopy(self.info['sample']))
                self.logger.warning(f"Run test assembly line for sample -> {sample.info['uuid']}\n")
                future = self.factory.submit_task(sample.params, sample.info)
                output = future.result()
                self.module_manager.database.add_data(output)
        except Exception as e:
            # 异常处理
            self.logger.error(f"An error occurred: {e}")

        finally:
            self.factory.executor.shutdown()
            self.module_manager.database.stop()

            from time import time
            start = time()
            self.module_manager.database.hdf5_to_csv()
            # self.module_manager.database.hdf5_to_csv(self.info['db']['out_csv'])
            tot = 1000 * (time() - start)
            self.logger.info(f"{tot} millisecond -> All samples have been processed.")
            # self.monitor.stop()

    def run_until_finished(self):
        self.sampler.set_factory(factory = self.factory)
        try:
            self.sampler.run_nested()
        finally:
            from time import time
            start = time()
            self.factory.module_manager.database.stop()
            self.factory.module_manager.database.hdf5_to_csv()
            # self.factory.module_manager.database.hdf5_to_csv(self.info['db']['out_csv'])
            tot = 1000 * (time() - start)
            self.logger.info(f"{tot} millisecond -> All samples have been processed.")
            self.sampler.finalize()
            self.sampler.combine_data(self.info['db']['path'])
            

            self.factory.executor.shutdown()
            # self.monitor.stop()

    def check_init_args(self) -> None:
        try:
            self.args = self.argparser.parse_args()
        except argparse.ArgumentError as e:
            print(str(e))
            self.argparser.print_help()
            sys.exit(2)

        if self.args.mkproject:
            if self.args.file:
                self.argparser.error("positional argument `file` cannot be used with --mkproject")
            if any([
                getattr(self.args, "plot", False),
                getattr(self.args, "cvtDB", False),
                getattr(self.args, "monitor", False),
                getattr(self.args, "OPC", False),
            ]):
                self.argparser.error("--mkproject cannot be combined with workflow mode options")
            return

        if getattr(self.args, "version", False):
            return

        if not self.args.file:
            self.argparser.error("the following arguments are required: file")

    def mkproject(self) -> None:
        try:
            project_root = create_project_scaffold(str(self.args.mkproject), cwd=os.getcwd())
        except ValueError as exc:
            print(f"[Jarvis-HEP] {exc}")
            sys.exit(2)
        except FileExistsError as exc:
            print(f"[Jarvis-HEP] Project directory already exists: {exc}")
            sys.exit(1)

        print(f"[Jarvis-HEP] Project scaffold created at: {project_root}")
        print(f"[Jarvis-HEP] Created folders: {', '.join(PROJECT_SUBDIRS)}")

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

    class __StateSaver:
        def __init__(self, 
                     obj, 
                     filename='my_object.pkl', 
                     save_interval_seconds=30, 
                     logger=logging.getLogger('MyClassStateSaver')
            ):
            self.obj = obj  # Save the outside object 
            self.filename = os.path.abspath(filename)
            self.save_interval_seconds = save_interval_seconds
            self.logger = logger
            self.timer = None
            self.start_auto_save()
            self.logger.warning(f"Started successfully, Jarvis-HEP create the storage station -> {self.filename}")

        def save_state(self):
            try:
                with open(self.filename, 'wb') as f:
                    dill.dump(self.obj, f)
                self.logger.info(f"Progress has been saved to hard disk space -> {self.filename}")
            except Exception as e:
                self.logger.error(f"Failed to save state: {e}")

        def start_auto_save(self):
            """Start auto saving mission"""
            self.timer = Timer(self.save_interval_seconds, self.auto_save)
            self.timer.daemon = True
            self.timer.start()

        def auto_save(self):
            """Auto saving the object, and restart the timer"""
            self.save_state()
            self.start_auto_save()  # save the state, restart the auto save 

        def stop_auto_save(self):
            """Stop auto saving method"""
            if self.timer is not None:
                self.timer.cancel()
                self.timer = None



def load_args_config(json_file):
    with open(json_file, 'r') as file:
        config = json.load(file)
    return config
