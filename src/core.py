#!/usr/bin/env python3

import configparser
from copy import deepcopy
import json
import os
from subprocess import Popen, run
import sys
import time 
import argparse
from config import ConfigLoader
import logging
from base import Base
from distributor import Distributor
import dill 
from threading import Timer
from library import Library
from workflow import Workflow
from factory import WorkerFactory 
from sample import Sample
from moduleManager import ModuleManager
from pprint import pprint
import pandas as pd 
import concurrent.futures
import asyncio 
from loguru import logger
logger.remove()
from plot import BudingPLOT
from monitor import Monitor


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
        self.plotter                    = BudingPLOT()
        self.mode                       = None
        # self.monitor                    = Monitor()
        # self.monitor.start()

    def init_argparser(self) -> None:
        self.argparser = argparse.ArgumentParser(description="Jarvis Program Help Center")
        self.info['args'] = load_args_config(self.path['args_info'])
        
        for pos_arg in self.info['args'].get('positionals', []):
            self.argparser.add_argument(pos_arg['name'], help=pos_arg['help'])
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

        
    def init_project(self) -> None: 
        import yaml 
        with open(os.path.abspath(self.args.file), 'r') as file:
            config = yaml.safe_load(file)
        self.info['scan_name'] = config['Scan']['name']
        self.info["project_name"] = os.path.splitext(os.path.basename(self.args.file))[0]
        task_result_dir = os.path.join(config['Scan']['save_dir'], self.info['scan_name'])
        task_result_dir = self.decode_path(task_result_dir)
        self.info['config_file'] = os.path.abspath(self.args.file)

        self.info['jarvis_log'] = os.path.join(task_result_dir, "LOG", f"{self.info['scan_name']}.log")
        self.info['pickle_file'] = os.path.join(task_result_dir, f"{self.info['project_name']}.pkl")
        self.info['flowchart_path'] = os.path.join(task_result_dir, "flowchart.png")
        self.info['sampler_log'] = os.path.join(task_result_dir, "LOG", f"{config['Sampling']['Method']}.log")

        self.info['sample'] = {
            "task_result_dir": self.decode_path(task_result_dir),
            "sample_dirs": os.path.join(task_result_dir, "SAMPLE"),
            "jarvis_log":  self.info['jarvis_log']
        }
        self.info["db"] = {
            "path":  os.path.join(task_result_dir, "DATABASE", "samples.hdf5"),
            "info":  os.path.join(task_result_dir, "DATABASE", "running.json")
            # "out_csv":  os.path.join(task_result_dir, "DATABASE", "samples.csv")
        }

        os.makedirs(task_result_dir, exist_ok=True)
        os.makedirs(os.path.join(task_result_dir, "SAMPLE"), exist_ok=True)
        os.makedirs(os.path.join(task_result_dir, "LOG"), exist_ok=True)
        with open(os.path.join(task_result_dir, "LOG", "monitor.log"), "w") as f1: 
            f1.write(os.getpid())
        
        os.makedirs(os.path.join(task_result_dir, "DATABASE"), exist_ok=True)
        if self.args.plot:
            self.info['plot'] = {
                "save_path":    os.path.join(task_result_dir, "IMAGE"),
                "config":       os.path.join(task_result_dir, "IMAGE", f"{self.info['scan_name']}.yaml")
            }
            os.makedirs(os.path.join(task_result_dir, "IMAGE"), exist_ok=True)
        
        # pprint(self.yaml.config)

    def init_logger(self) -> None:
        from datetime import datetime
        current_time = datetime.now().strftime("%Y-%m-%d[%H:%M:%S]")

        def global_log_filter(record):
            # 只有包含 'global_log' 标记的日志消息才被写入
            return record["extra"].get("Jarvis", False)

        def stream_filter(record):
            # 检查日志记录是否包含 'to_console' 标记，并且该标记为 True
            return record["extra"].get("to_console", False)

        def custom_format(record):
            module = record["extra"].get("module", "No module")
            if "raw" in record["extra"]:
                return "{message}\n"
            elif module == "Jarvis-HEP.hdf5-Writter":
                return f"\nϠ <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{record['message']}</level> "
            elif "Sample@" in module:
                return f"\n·•· <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{record['message']}</level>"
            else:
                return f"\n <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{record['message']}</level> "

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

        self.logger = logger.bind(module="Jarvis-HEP", to_console=True, Jarvis=True)
        with open(self.path['logo'], 'r') as f:
            self.logger.warning(f"\n{f.read()}")        
        self.logger.warning("Jarvis-HEP logging system initialized successful!")
        if self.args.debug:
            self.logger.info(f"Jarvis-HEP write into main log file -> {jarvislog}")
            self.logger.info("Jarvis-HEP in debug mode currently!")
        if self.args.plot:
            plogger = logger.bind(module="Jarvis-PLOT", to_console=True, Jarvis=False)
            self.plotter.logger = plogger

    def init_configparser(self) -> None: 
        self.yaml.logger = logger.bind(module="Jarvis-HEP.ConfigParser", to_console=True, Jarvis=True)
        from copy import deepcopy
        self.yaml.path = deepcopy(self.path)
        self.yaml.load_config(os.path.abspath(self.args.file))
        self.yaml.check_dependency_installed()
        self.sampler = Distributor.set_method(self.yaml.get_sampling_method()) 
        self.yaml.set_schema(self.sampler.schema)
        self.yaml.validate_config()

    def init_utils(self) -> None: 
        if "Utils" in self.yaml.config:
            if "interpolations_1D" in self.yaml.config['Utils']:
                from utils import get_interpolate_1D_function_from_config
                for item in self.yaml.config['Utils']['interpolations_1D']:
                    if "file" in item:
                        item['file'] = self.decode_path(item['file'])
                    func = get_interpolate_1D_function_from_config(item)
                    if func is not None: 
                        self._funcs[item['name']] = func
                        self.logger.info(f"Succefully resolve the interpolate function -> {item['name']}")
                    else: 
                        self.logger.info(f"Illegal setting for interpolate function -> {item['name']}")

    def init_StateSaver(self) -> None:
        # logger = self.logger.create_dynamic_logger("StateSaver", logging.INFO)
        slogger = logger.bind(module="Jarvis-HEP.StateSaver", to_console=True, Jarvis=True)
        slogger.warning("Enabling breakpoint resume function ... ")
        self.__state_saver = self.__StateSaver(self, filename=self.info['pickle_file'] , logger=slogger, save_interval_seconds=60)

    def init_sampler(self) -> None:
        logger_name = f"Jarvis-HEP.{self.sampler.method}"
        def filte_func(record):
            return record['extra']['module'] == logger_name
        def custom_format(record):
            module = record["extra"].get("module", "No module")
            return f"\n·•· <red>{module}</red> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{record['message']}</level>"
    
        self.sampler.set_config(self.yaml.config)
        # logger = self.logger.create_dynamic_logger(self.sampler.method)
        slogger = logger.bind(module=logger_name, to_console=True, Jarvis=True)
        slogger.add(self.info['sampler_log'], format=custom_format, level="DEBUG", rotation=None, retention=None, filter=filte_func )
        self.sampler.info['logfile'] = self.info['sampler_log']
        self.sampler.set_logger(slogger)
        self.sampler.initialize()
        self.yaml.vars = self.sampler.vars 
        self.sampler.info['sample'] = deepcopy(self.info['sample'])

    def init_librarys(self) -> None:
        self.libraries = Library()
        self.libraries._skip_library = self.args.skiplibrary
        # logger = self.logger.create_dynamic_logger("Library", logging.INFO)
        slogger = logger.bind(module="Jarvis-HEP.Library", to_console=True, Jarvis=True)
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
            from threading import Thread
            self.logger.warning(f"Draw workflow chart into {self.info['flowchart_path']}")
            asyncio.run(
                self.workflow.draw_flowchart(save_path=self.info['flowchart_path'])
            )
        self.workflow.get_workflow_dict()
        self.logger.warning(self.workflow.workflow.keys())

    def init_WorkerFactory(self) -> None: 
        self.factory = WorkerFactory()
        self.module_manager = ModuleManager()
        self.factory.configure(module_manager=self.module_manager,
            max_workers=self.yaml.config['Calculators']['make_paraller']
            )
        self.sampler.set_max_workers(self.yaml.config['Calculators']['make_paraller'])
        flogger = logger.bind(module="Jarvis-HEP.Factory", to_console=True, Jarvis=True)
        self.factory.set_logger(flogger)
        mlogger = logger.bind(module="Jarvis-HEP.Factory.Manager", to_console=True, Jarvis=True)
        self.module_manager.set_logger(mlogger)
        self.module_manager.set_max_worker(self.yaml.config['Calculators']['make_paraller'])
        self.module_manager.set_config(self.yaml.config)
        self.module_manager.set_funcs(self._funcs)
        self.module_manager.workflow = deepcopy(self.workflow.workflow)
        for kk, layer in self.workflow.calc_layer.items():
            if kk > 1: 
                for module in layer['module']:
                    self.module_manager.add_module_pool(self.workflow.modules[module])
        self.factory.info['sample'] = deepcopy(self.info['sample'])

    def init_likelihood(self) -> None: 
        self.module_manager.set_likelihood()
        self.module_manager.set_funcs(self._funcs)

    def init_database(self) -> None: 
        from hdf5writer import GlobalHDF5Writer
        self.module_manager._database = GlobalHDF5Writer(self.info['db'])
        self.module_manager.database.start()

    def initialization(self) -> None:
        self.init_argparser()
        self.init_project()
        self.init_logger()
        self.init_configparser()
        self.init_utils()
        
        if self.scan_mode:
            self.init_StateSaver()
            
            self.init_sampler()
            self.init_workflow()
            self.init_librarys()
            self.init_WorkerFactory()
            self.init_likelihood()
            self.init_database()
        # elif self.args.cvtDB:
        #     self.convert()

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
                        from utils import convert_hdf5_to_csv
                        convert_hdf5_to_csv(snapshot_path, out_csv)
                        self.logger.warning(f"Data has been successfully written to -> {out_csv}")
                    finally:
                        try:
                            os.remove(snapshot_path)
                            self.logger.warning(f"Snapshot -> {snapshot_path} has been successfully deleted.")
                        except Exception as e:
                            self.logger.error(f"Failed to delete snapshot {snapshot_path}: {str(e)}")

    def plot(self) -> None:
        self.plotter.logger.warning(f"Data visulization for {self.info['project_name']}")
        if "Plot" not in self.yaml.config:
            if not os.path.exists(self.info['plot']['config']):
                self.plotter.get_plot_config_from_Jarvis(self.info, self.yaml)
                self.plotter.plot()
            else:
                self.plotter.load_config(self.info['plot']['config'])
                self.plotter.plot()
        

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
