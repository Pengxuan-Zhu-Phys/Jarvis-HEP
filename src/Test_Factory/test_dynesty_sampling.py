#!/usr/bin/env python3

import numpy as np
import dynesty
from concurrent.futures import ProcessPoolExecutor
import time
from random import random

# 模拟的外部命令执行函数，计算多变量高斯的对数似然值
def simulate_external_likelihood_calculation(params):
    mean = np.zeros(len(params))
    cov = np.eye(len(params))
    # time.sleep(0.1+ 0.05 * random())
    return -0.5 * np.dot((params-mean).T, np.linalg.solve(cov, (params-mean)))

# WorkerFactory类，用于管理并行的Likelihood计算任务
# class WorkerFactory:
#     def __init__(self, max_workers=4):
#         self.executor = ProcessPoolExecutor(max_workers=max_workers)

#     def compute_likelihood_async(self, params):
#         return self.executor.submit(simulate_external_likelihood_calculation, params)

from concurrent.futures import ProcessPoolExecutor
import threading

class WorkerFactory:
    _instance = None
    _lock = threading.Lock()
    
    def __new__(cls, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(WorkerFactory, cls).__new__(cls)
                # 假设初始化逻辑在这里
                cls._instance.executor = ProcessPoolExecutor(max_workers=4)
        return cls._instance

    def compute_likelihood(self, params):
        # 这里调用外部程序计算likelihood
        # 假设这是外部程序的调用结果
        result = simulate_external_likelihood_calculation(params)
        return result


# 先验转换函数
def prior_transform(utheta):
    return 10. * utheta - 5.

# log_likelihood函数，通过WorkerFactory提交任务并等待结果
# def log_likelihood(theta):
#     # print("Start run {} - {}".format(theta, time.time()))
#     future = factory.compute_likelihood_async(theta)
#     result = future.result()  # 阻塞等待计算结果
#     # print("Finish run {} - {} - {}".format(theta, time.time(), result))
#     return result

# 全局作用域下初始化WorkerFactory单例
factory = WorkerFactory()

def log_likelihood(params):
    # 直接使用factory实例计算likelihood
    start = time.time()
    print("Start run {} - {}".format(params, start))

    result = factory.compute_likelihood(params)
    print("Finish run {} - {} - {}".format(params, time.time() - start, result))
    # result = future.result()
    return result
# def create_log_likelihood(factory):
#     def log_likelihood(params):
#         # 使用factory计算并返回likelihood
#         result = factory.compute_likelihood_async(params)
#         return result
#     return log_likelihood


if __name__ == "__main__":
    # 初始化WorkerFactory
    # factory = WorkerFactory(max_workers=2)
    # log_likelihood = create_log_likelihood(factory)
    # 使用ProcessPoolExecutor创建并行池
    with ProcessPoolExecutor(max_workers=2) as executor:
        # 包装executor以适配dynesty所需的接口
        # pool = dynesty.utils.Pool(processes=2, pool=executor)

        # 初始化并行的dynesty采样器
        ndim = 2
        sampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, ndim, pool=executor, queue_size=4)

        # sampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, ndim, pool=pool, queue_size=2)
        sampler.run_nested(dlogz_init=0.01, print_progress=False)
        results = sampler.results

    # 关闭工厂
    factory.executor.shutdown()

    print("Sampling completed.")
