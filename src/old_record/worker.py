#!/usr/bin/env python3

from base import Base

import uuid
from concurrent.futures import Future

class Worker(Base):
    def __init__(self, module_pool):
        self.id = uuid.uuid4()
        self.module_pool = module_pool
        self.future = None  # 用于存储异步执行结果

    def execute(self, module_name, params):
        # 根据模块名称从池中获取模块并执行计算任务
        module = self.module_pool.get_module(module_name)
        self.future = module.submit_task(params)
        return self.future