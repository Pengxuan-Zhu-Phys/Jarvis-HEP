#!/usr/bin/env python3

import numpy as np
# import dynesty
from py.dynesty import dynesty
from concurrent.futures import ThreadPoolExecutor
import time
from random import random

# 模拟的外部命令执行函数，计算多变量高斯的对数似然值
def simulate_external_likelihood_calculation(params):
    mean = np.zeros(len(params))
    cov = np.eye(len(params))
    # time.sleep(0.0001+ 0.00005 * random())
    return -0.5 * np.dot((params-mean).T, np.linalg.solve(cov, (params-mean)))

# WorkerFactory类，用于管理并行的Likelihood计算任务
# class WorkerFactory:
#     def __init__(self, max_workers=4):
#         self.executor = ProcessPoolExecutor(max_workers=max_workers)

#     def compute_likelihood_async(self, params):
#         return self.executor.submit(simulate_external_likelihood_calculation, params)
from uuid import uuid4
import threading

class WorkerFactory:
    _instance = None
    _lock = threading.Lock()
    
    def __new__(cls, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(WorkerFactory, cls).__new__(cls)
                # 假设初始化逻辑在这里
                cls._instance.executor = ThreadPoolExecutor(max_workers=4)
        return cls._instance

    def compute_likelihood(self, params):
        # 这里调用外部程序计算likelihood
        # 假设这是外部程序的调用结果
        result = simulate_external_likelihood_calculation(params)
        return result


# 先验转换函数
def prior_transform(utheta):
    u = 10. * utheta - 5.
    ret = np.append(u, str(uuid4()))
    return ret

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
    param = params[:-1].astype(float)

    uuid = params[-1]
    start = time.time()
    # print("{} Start run {} - {}".format(uuid, param, start))

    result = factory.compute_likelihood(param)
    # print("{} Finish run {} - {} - {}".format(uuid, param, time.time() - start, result))
    # result = future.result()
    return result
# def create_log_likelihood(factory):
#     def log_likelihood(params):
#         # 使用factory计算并返回likelihood
#         result = factory.compute_likelihood_async(params)
#         return result
#     return log_likelihood

import pandas as pd

def save_dynesty_results_to_csv(results, csv_file_path):
    """
    将dynesty的结果保存为CSV文件。

    参数:
        results: dynesty的结果对象。
        csv_file_path: 保存CSV文件的路径。
    """
    # 创建一个字典来存储我们感兴趣的结果数据
    data = {
        'samples_uid': results['samples_uid'],  # 假设结果对象中已经有了sample_uid
        'logl': results['logl'],  # 对数似然值
        'logwt': results['logwt'],  # 样本权重
    }
    
    # 假设样本参数是多维的，我们需要为每个维度创建一个列
    for i in range(results['samples'].shape[1]):  # 假设samples是二维数组，形状为(niter, ndim)
        data[f'param_{i}'] = results['samples'][:, i]
    
    # 使用pandas创建数据帧
    df = pd.DataFrame(data)
    
    # 保存到CSV文件
    df.to_csv(csv_file_path, index=False)
    print(f"Results saved to {csv_file_path}")



if __name__ == "__main__":
    # 初始化WorkerFactory
    # factory = WorkerFactory(max_workers=2)
    # log_likelihood = create_log_likelihood(factory)
    # 使用ProcessPoolExecutor创建并行池
    with ThreadPoolExecutor(max_workers=16) as executor:
        # 包装executor以适配dynesty所需的接口
        # pool = dynesty.utils.Pool(processes=2, pool=executor)

        # 初始化并行的dynesty采样器
        ndim = 2
        sampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, ndim, pool=executor, queue_size=32)

        # sampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, ndim, pool=pool, queue_size=2)
        sampler.run_nested(dlogz_init=0.01, print_progress=False)
        results = sampler.results
        save_dynesty_results_to_csv(results, "test_dynesty.csv")
        # print(results)
        print("Information", results['information'])
        print("logZ length", len(results['logz']))
        print(results.summary())
        print(results.isdynamic())
        
        import matplotlib.pyplot as plt 
        fig = plt.figure(figsize=(10, 4))
        ax = fig.add_axes([0.15, 0.15, 0.84, 0.84])

        ax.plot(- results['logvol'], results['samples_n'], 'o', color='blue' )
        # plt.show()
        plt.savefig("test_dynesty_02.png")
        # print(""results['logvol']))


    # 关闭工厂
    factory.executor.shutdown()

    print("Sampling completed.")
