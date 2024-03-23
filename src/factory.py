#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor, as_completed
import threading
from modulePool import ModulePool
from pprint import pprint
from sample import Sample
import logging
import threading
from concurrent.futures import ThreadPoolExecutor

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
            self.initialized = True
            # 其他一次性初始化逻辑...


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

    def add_module(self, module, logger):
        if module.name not in self.module_pools:
            self.logger.warning(f"Adding ModulePool {module.name}. ")
            self.module_pools[module.name] = ModulePool(module, max_workers=self._max_workers)
            self.module_pools[module.name].set_logger(logger)
            self.module_pools[module.name].load_installed_instances()



        else:
            self.logger.warning(f"ModulePool for {module.name} already exists.")

    def submit_task(self, params, uuid):
        # 这个方法使用ModuleManager来执行一个特定的工作流
        # 通过ThreadPoolExecutor来异步执行
        future = self.executor.submit(self.module_manager.execute_workflow, params, uuid)
        return future

    def task_done(self, future, uuid):
        # 此处处理任务完成的逻辑，不论成功还是失败
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


    def compute_likelihood(self, sample):
        # 直接操作sample对象，包括其logger
        print("Testing the compute likelihood")
        # sample.logger.warning(f"Compute the sample {sample.params}")
        observables = sample.params
        for layer in self.workflow.values():
            for mod in layer:
                output = self.module_pools[mod].submit_task(observables)
                # 等待并处理每个模块的输出结果
                observables.update(output.result())  # 假设output.result()返回了所需的更新字典
        
        # 更新sample的状态或结果
        sample.likelihood = 0.5
        return sample