#!/usr/bin/env python3

from concurrent.futures import ThreadPoolExecutor
from Module.parameters import Parameters
from Module.library import LibraryModule
from Module.calculator import CalculatorModule
from copy import deepcopy
from concurrent.futures import as_completed
import json
import os, sys 
import threading 

class ModulePool:
    def __init__(self, module, max_workers=2):
        self.name = module.name 
        self.config = module.config 
        self.executor = ThreadPoolExecutor(max_workers=max_workers)
        self.module_name = module.name 
        self.type = module.type 
        self.instances = []
        self.id_counter = 0 
        self.instances_info_file = self.get_instances_info_file_path()
        self._info_loaded = False 
        self._funcs = {}

    def set_logger(self, logger):
        self.logger = logger
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
        self.instances.append(instance)
        self.update_instances_info_file()
        return instance

    def reload_instance(self, pack_id, pack_info):
        self.id_counter += 1
        config = deepcopy(self.config)
        self.logger.info(f"Re-loading the {pack_id} Instance for Module {self.name}")
        instance = CalculatorModule(self.name, config=config)
        instance.assign_ID(pack_id)
        instance.is_installed = True
        instance.is_busy = False
        instance._funcs = self.funcs
        self.instances.append(instance)
        self.update_instances_info_file()
        return instance

    def set_funcs(self, funcs):
        self._funcs = funcs
        # self.logger.warning(f"ModulePool {self.name} funcs -> {self.funcs}")

    @property
    def funcs(self):
        return self._funcs

    def execute(self, params, sample_info):
        """Submit parameters to a module instance for calculation and return the calculation results.

        Args:
            params (dict): The parameters required for the calculation.

        Returns:
            dict: The updated observables dictionary.
        """

        instance = self.get_available_instance()
        # self.logger.warning(f"Line 64 -> {instance}")
        if instance is None:
            # if no availiable instance found, create a new one!
            self.logger.warning(f"No availiable instance for Module {self.name} found, trying to install a new one")
            instance = self.create_instance()
            self.executor.submit(self.install_instance, instance)
            instance.installation_event.wait()

        instance.is_busy = True 
        # submit the task and tagging the instance into busy
        future = self.executor.submit(instance.execute, params, sample_info)

        # waiting for calculation finished
        result = future.result()        
        # Update instance status 
        instance.is_busy = False
        self.update_instances_info_file()

        return result


    @staticmethod
    def install_instance(instance):
        if not instance.is_installed:
            instance.install()
            instance.is_installed = True 
            instance.installation_event.set()
        return instance

    def submit_task(self, params):
        instance = self.get_available_instance()
        if not instance:
            instance = self.create_instance()
            self.install_instances()  # Ensure the instance are installed 
        future = self.executor.submit(instance.execute, params)
        instance.is_busy = True
        # print("Line 77, modulePool -> Test after the task is summited ... ")
        self.update_instances_info_file()
        return future

    def get_available_instance(self):
        # print(self.instances)
        for instance in self.instances:
            if not instance.is_busy and instance.is_installed:
                return instance
        return None

    def update_instances_info_file(self):
        installed_instances_info = {
            instance.PackID: {
                "is_installed": instance.is_installed,
                "is_busy": instance.is_busy
            } for instance in self.instances
        }
        data_to_save = {"installed_instances": installed_instances_info}
        with open(self.instances_info_file, 'w') as file:
            json.dump(data_to_save, file, indent=4)

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
                    installed_ids = installed_instances_info.get("installed_instances", [])
                    for pack_id, info in installed_ids.items():
                        self.reload_instance(pack_id=pack_id, pack_info=info)
                    self.logger.warning("Instance reloaded!")
            except FileNotFoundError:
                self.logger.info(f"{self.instances_info_file} not found. Assuming no instances are installed.")
        self._info_loaded = True

    def get_instances_info_file_path(self):
        return os.path.join(self.config['path'].replace("/@PackID", ""), f"{self.name}_instance_info.json")

    def init_instances_info_file(self):
        if not os.path.exists(os.path.dirname(self.instances_info_file)):
            # self.logger.warning("Init instance info file ")
            os.makedirs(os.path.dirname(self.instances_info_file))
        data_to_save = {"installed_instances": {}}
        with open(self.instances_info_file, 'w') as file:
            json.dump(data_to_save, file, indent=4)
