#!/usr/bin/env python3 

import numpy as np

# jDE控制参数
tau = 0.1

# gencrossover: 选择二进制或指数交叉来生成试验向量
def gencrossover(X, V, n, run_params):
    if run_params['jDE']:
        trialCr = newCr(X, n)
        U = bincrossover(X, V, n, run_params, trialCr)
    elif run_params['expon']:
        U = expcrossover(X, V, n, run_params)
    else:
        U = bincrossover(X, V, n, run_params)
    return U, trialCr

# 二进制交叉函数
def bincrossover(X, V, n, run_params, trialCr=None):
    D = run_params['D']
    U = np.copy(X['vectors'][n, :])  # 初始化为目标向量
    Cr = trialCr if trialCr is not None else run_params['Cr']

    # 确保至少有一个维度发生交叉
    jrand = np.random.randint(0, D)

    randj = np.random.rand(D)
    U[randj <= Cr] = V[randj <= Cr]  # 使用供体向量
    U[jrand] = V[jrand]  # 保证交叉点

    return U

# 指数交叉函数
def expcrossover(X, V, n, run_params):
    D = run_params['D']
    U = np.copy(X['vectors'][n, :])
    
    # 随机选择交叉起点
    j = np.random.randint(0, D)
    L = 0

    # 确定交叉的长度
    while L < D:
        L += 1
        if np.random.rand() > run_params['Cr']:
            break

    # 实现循环交叉
    for k in range(L):
        U[(j + k) % D] = V[(j + k) % D]

    return U

# jDE方法中的Cr自适应调整函数
def newCr(X, n):
    if np.random.rand() < tau:
        return np.random.rand()  # 随机生成Cr
    else:
        return X['CrjDE'][n]  # 使用上一代的Cr值

# 初始化CrjDE
def init_CrjDE(size):
    return np.random.rand(size)

# 生成随机整数
def random_int(min_value, max_value):
    return np.random.randint(min_value, max_value + 1)
