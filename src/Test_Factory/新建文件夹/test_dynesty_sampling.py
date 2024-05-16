#!/usr/bin/env python3

import numpy as np
# import dynesty
from py.dynesty import dynesty
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time
from random import random

tmax = 6.0 * np.pi
def loglike(x):
    t = 2.0 * tmax * x - tmax
    return (2.0 + np.cos(t[0] / 2.0) * np.cos(t[1] / 2.0)) ** 5.0

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
                cls._instance.executor = ProcessPoolExecutor(max_workers=16)
        return cls._instance

    def compute_likelihood(self, params):
        # 这里调用外部程序计算likelihood
        # 假设这是外部程序的调用结果
        result = loglike(params)
        return result


# 先验转换函数
def prior_transform(x):
    # return x
    ret = np.append(x, str(uuid4()))
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
    import pandas as pd 
    data = {
        "uuid":             results['samples_uid'],
        "log_weight":       results['logwt'],
        "log_Like":         results['logl'],
        "log_PriorVolume":  results['logvol'],
        "log_Evidence":     results['logz'],
        "log_Evidence_err": results['logzerr'],
        "samples_nlive":    results['samples_n'],
        "ncall":            results['ncall'],
        "samples_it":       results['samples_it'],
        "samples_id":       results['samples_id'],
        "information":      results['information']
    }
    for ii in range(results['samples'].shape[1]):
        data[f'samples_v[{ii}]'] = results['samples'][:, ii]
    for ii in range(results['samples_u'].shape[1]):
        data[f'samples_u[{ii}]'] = results['samples'][:, ii]
    df = pd.DataFrame(data)
    df.to_csv(csv_file_path, index=False)
    print(f"Results saved to {csv_file_path}")
    print(results.summary())



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
        sampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, ndim, pool=executor, queue_size=32, log_file_path="dynesty.log")

        # sampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, ndim, pool=pool, queue_size=2)
        # sampler.run_nested(dlogz_init=0.01, print_progress=False, nlive_init=50, wt_kwargs={'pfrac': 1.0})
        sampler.run_nested(print_progress=True, wt_kwargs={'pfrac': 0.8}, dlogz_init=0.5, nlive_init=360)
        results = sampler.results
        from dynesty import plotting as dyplot
        fig, axes = dyplot.runplot(results, kde=False) 
        fig.savefig("save01.png", dpi=300)
        fig, axes = dyplot.runplot(results, kde=True) 
        fig.savefig("save02.png", dpi=300)

        save_dynesty_results_to_csv(results, "test_dynesty_02.csv")
        print(results.keys())
        # print("Information", results['information'])
        # print("logZ length", len(results['logz']))
        # print("samples Number", len(results['samples_uid']))
        for kk, vv in results.items():
            if type(vv) == np.ndarray:
                print(f"Length {kk} \t->  {type(vv[0])}\t {vv.shape}")
        
        print("====================")
        for kk, vv in results.items():
            if type(vv) != np.ndarray:
                print(f"{kk} \t->\t{vv}")
        print(results.summary())
        print(results.isdynamic())
        
        import matplotlib.pyplot as plt 
        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_axes([0.15, 0.15, 0.84, 0.84])
        ax.plot(- results['logvol'], results['samples_n'], '.', color='blue', markersize=0.5, alpha=0.6 )
        plt.savefig("test_dynesty_d_01.png", dpi=300)
        print("figure 1 succeed")

        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_axes([0.15, 0.15, 0.84, 0.84])
        ax.plot(- results['logvol'], np.exp(results['logl'] - max(results['logl'])), '.', color='blue', markersize=0.5, alpha=0.6 )
        plt.savefig("test_dynesty_d_02.png", dpi=300)
        print("figure 2 succeed")


        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_axes([0.15, 0.15, 0.84, 0.84])
        ax.plot(- results['logvol'], np.exp(results['logwt'] - results['logz'][-1]), '.', color='blue', markersize=0.5, alpha=0.6 )
        plt.savefig("test_dynesty_d_03.png", dpi=300)
        print("figure 3 succeed")

        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_axes([0.15, 0.15, 0.84, 0.84])
        ax.plot(- results['logvol'], np.exp(results['logz']), '.', color='blue', markersize=0.5, alpha=0.6 )
        plt.savefig("test_dynesty_d_04.png", dpi=300)
        print("figure 4 succeed")

        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_axes([0.15, 0.15, 0.84, 0.84])
        ax.plot(- results['logvol'], results['samples_it'], '.', color='blue', markersize=0.5, alpha=0.6 )
        plt.savefig("test_dynesty_d_05.png", dpi=300)
        print("figure 5 succeed")

        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_axes([0.15, 0.15, 0.84, 0.84])
        ax.plot(- results['logvol'], results['ncall'], '.', color='blue', markersize=0.5, alpha=0.6 )
        plt.savefig("test_dynesty_d_06.png", dpi=300)
        print("figure 6 succeed")

        fig = plt.figure(figsize=(10, 6))
        ax1 = fig.add_axes([0.15, 0.2, 0.84, 0.4])
        ax2 = fig.add_axes([0.15, 0.6, 0.84, 0.4])
        ax1.plot(- results['logvol'], results['samples_u'][:, 0], '.', color='blue', markersize=0.5, alpha=0.6 )
        ax2.plot(- results['logvol'], results['samples_u'][:, 1], '.', color='blue', markersize=0.5, alpha=0.6 )
        plt.savefig("test_dynesty_d_07.png", dpi=300)
        print("figure 7 succeed")

        # print(""results['logvol']))


    # 关闭工厂
    factory.executor.shutdown()

    print("Sampling completed.")
