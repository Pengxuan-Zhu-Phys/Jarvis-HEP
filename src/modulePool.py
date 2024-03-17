#!/usr/bin/env python3

from concurrent.futures import ThreadPoolExecutor
from module import Parameters, CalculatorModule
from copy import deepcopy
from concurrent.futures import as_completed

class ModulePool:
    def __init__(self, module, max_workers=2):
        self.name = module.name 
        self.config = module.config 
        self.executor = ThreadPoolExecutor(max_workers=max_workers)
        self.module_name = module.name 
        self.type = module.type 
        self.instances = []
        self.id_counter = 0  # ID计数器


    def create_instance(self):
        self.id_counter += 1
        id = f"{self.id_counter:03}"
        config = deepcopy(self.config)
        self.logger.info(f"Create the {id} Instance for Module {self.name}")
        instance = CalculatorModule(self.name, config=config)
        instance.assign_ID(id)
        instance.logger = self.logger
        self.instances.append(instance)
        return instance


    def install_instances(self):
        futures = []
        for instance in self.instances:
            if not instance.is_installed:
                future = self.executor.submit(self.install_instance, instance)
                futures.append(future)
        from time import time
        for future in as_completed(futures):
            instance = future.result()
            instance.is_installed = True

    @staticmethod
    def install_instance(instance):
        instance.install()
        return instance

    def submit_task(self, params):
        instance = self.get_available_instance()
        if not instance:
            instance = self.create_instance()
            self.install_instances()  # Ensure the instance are installed 
        future = self.executor.submit(instance.execute, params)
        instance.is_busy = True
        return future

    def get_available_instance(self):
        for instance in self.instances:
            if not instance.is_busy and instance.is_installed:
                return instance
        return None