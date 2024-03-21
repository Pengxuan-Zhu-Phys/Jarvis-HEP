#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor
import threading
from modulePool import ModulePool
from pprint import pprint


class WorkerFactory:
    _instance = None
    _lock = threading.Lock()
    _max_workers = 4  # 默认的最大工作线程数

    def __new__(cls, max_workers=None, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(WorkerFactory, cls).__new__(cls)
                if max_workers is not None:
                    cls._max_workers = max_workers
                cls._instance.executor = ProcessPoolExecutor(max_workers=cls._max_workers)
        return cls._instance

    def __init__(self, max_workers=None):
        self.workflow = None  # Initialize with None
        self.active_task = 0
        self.module_pools = {}
        self.layer = []
        self.config = None

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

    def submit_task(self, params):
        # 提交任务前增加活跃任务计数
        self.active_tasks += 1
        future = self.executor.submit(self.likelihood, params)
        future.add_done_callback(self.task_done)
        return future

    def task_done(self, future):
        # 任务完成时减少活跃任务计数
        self.active_tasks -= 1
        

    def setup_workflow(self, workflow):
        self.workflow = workflow  # Assign the workflow

    def add_module_to_pool(self, module):
        self.module_pool.add_module(module)


    def set_likelihood_func(self, func):
        # 设置用于计算的外部函数
        self.likelihood = func

    def compute_likelihood(self, params):
        if not self.workflow:
            raise ValueError("Workflow is not set up.")
        # Use the workflow to compute the likelihood
        result = self.workflow.run(params)  # Assuming 'run' is the method that executes the workflow and computes the result
        return result

    # def compute_likelihood(self, params):
    #     result = simulate_external_likelihood_calculation(params)
    #     return result



def sampler_generator(factory, max_active_tasks, sampler):
    while True:
        # 检查活跃任务数是否小于最大允许并行任务数
        if factory.get_active_tasks_count() < max_active_tasks:
            # 获取下一个采样点
            params = next(sampler)
            # 提交新的任务
            yield factory.submit_task(params)
        else:
            # 暂时没有可用的工作线程，稍后再试
            time.sleep(0.1)
