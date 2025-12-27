#!/usr/bin/env python3
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import threading
from modulePool import ModulePool
from pprint import pprint
from sample import Sample
import logging
import threading
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import time 

class WorkerFactory:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(WorkerFactory, cls).__new__(cls)
                cls._instance.module_manager = args[0] if args else None
                cls._instance.executor = None  # 延迟executor的创建
        return cls._instance

    def configure(self, module_manager=None, max_workers=4):
        if not hasattr(self, 'initialized'):
            self.module_manager = module_manager if module_manager else self.module_manager
            self.executor = ThreadPoolExecutor(max_workers=max_workers)
            self.task_count = 0  
            self.task_time = 0  
            self.last_time_mark = time.time()  
            self.initialized = True

    def get_executor(self):
        if self.executor is None:
            self.executor = ProcessPoolExecutor(max_workers=self._max_workers)
        return self.executor

    def set_config(self, config):
        self.config = config

    def add_layer(self, layer):
        self.layer.append(layer)

    def set_logger(self, logger):
        self.logger = logger
        self.logger.warning("Building the factory for workers ...")
        self.info   = {}

    def add_module(self, module, logger):
        if module.name not in self.module_pools:
            self.logger.warning(f"Adding ModulePool {module.name}. ")
            self.module_pools[module.name] = ModulePool(module, max_workers=self._max_workers)
            self.module_pools[module.name].set_logger(logger)
            self.module_pools[module.name].load_installed_instances()
        else:
            self.logger.warning(f"ModulePool for {module.name} already exists.")

    # def submit_task(self, params, sample_info):
    def submit_task(self, sample_info):
        # This method uses ModuleManager to execute a specific workflow
        # Asynchronous execution through ThreadPoolExecutor
        self.logger.info(f"Submit Task {sample_info['uuid']} into WorkerFactory ...")
        future = self.executor.submit(self.module_manager.execute_workflow, sample_info)
        self.task_count += 1
        self.print_status()
        return future

    def print_status(self):
        def format_duration(seconds):
            """Formating into HH:MM:SS.msc format"""
            hours = seconds // 3600
            minutes = (seconds % 3600) // 60
            seconds = seconds % 60
            millisec = seconds % 1
            return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}.{str(millisec)[2:5]}"


        if self.task_count % 100 == 0:
            end_time = time.time()  
            elapsed_time = format_duration(end_time - self.last_time_mark)  
            self.task_time += end_time - self.last_time_mark  
            self.logger.warning(f"Submitted {self.task_count} tasks, time for last 100 tasks: {elapsed_time}, total time: {format_duration(self.task_time)}")
            self.last_time_mark = end_time  


    def task_done(self, future, uuid):
        # Deal the thing after the task finished 
        del self.active_tasks[uuid]
        try:
            result = future.result()
            print(f"Task {uuid} result: {result}")
        except Exception as e:
            print(f"Task {uuid} failed with error: {e}")

    def setup_workflow(self, workflow):
        self.workflow = workflow  # Assign the workflow

    def add_module_to_pool(self, module):
        self.module_pool.add_module(module)

    def set_likelihood_func(self, func):
        # 设置用于计算的外部函数
        self.likelihood = func

    def submit_task_with_prior(self, params):
        # This method uses ModuleManager to execute a specific workflow
        # Asynchronous execution through ThreadPoolExecutor
        sample = Sample(params)
        sample.set_config(deepcopy(self.info['sample']))
        future = self.submit_task(sample.params, sample.info)
        return future.result()


