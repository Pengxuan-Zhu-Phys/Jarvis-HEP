#!/usr/bin/env python3 

import numpy as np

# 定义常量
meanimprovement = 0  # 用于收敛准则
checkpopres = False  # 是否检查种群的分辨率

# 初始化收敛准则
def init_convergence(run_params):
    if run_params['convergence_criterion'] == meanimprovement:
        run_params['meanlike'] = np.inf
        run_params['improvements'] = 1.0

# 检查证据是否收敛
def evidence_done(Z, Zerr, tol):
    return np.log(Z / (Z - Zerr)) <= tol

# 检查种群是否收敛
def converged(X, run_params):
    if run_params['verbose'] >= 3:
        print("  Checking convergence...")

    if run_params['convergence_criterion'] == meanimprovement:
        is_converged = check_SFIM(X, run_params)
    else:
        is_converged = False  # 未实现其他收敛准则

    if is_converged:
        if run_params['verbose'] >= 3:
            print("  Converged.")
    else:
        if run_params['verbose'] >= 3:
            print("  Not converged.")

    # 检查种群分辨率
    if checkpopres:
        check_population_resolution(X, run_params, is_converged)
    
    return is_converged

# 计算平滑的均值改进，判断是否收敛
def check_SFIM(X, run_params):
    curval = np.mean(X['values'])  # 当前种群的平均适应度
    inf_threshold = 0.001 * np.finfo(float).max  # 防止数值无限大的问题

    # 确保没有遇到无穷大问题
    if curval > inf_threshold:
        run_params['meanlike'] = inf_threshold
        fracdiff = 1.0
    else:
        fracdiff = 1.0 - curval / run_params['meanlike']  # 计算两代之间的改进
        run_params['meanlike'] = curval

    # 存储新的改进值并丢弃最早的改进值
    run_params['improvements'] = np.roll(run_params['improvements'], shift=-1)
    run_params['improvements'][-1] = fracdiff

    # 计算平滑的改进值
    sfim = np.mean(run_params['improvements'])

    if run_params['verbose'] >= 3:
        print(f"  Smoothed fractional improvement of the mean = {sfim}")

    # 比较改进值和收敛阈值
    return sfim < run_params['convthresh']

# 检查种群分辨率是否过小
def check_population_resolution(X, run_params, is_converged):
    avg_vector = np.mean(X['vectors'], axis=0)  # 计算种群的平均向量
    diff_vectors = np.abs(X['vectors'] - avg_vector)  # 计算每个个体与平均向量的差异

    if run_params['verbose'] >= 3:
        print("  Checking population resolution...")
        print("  Average vector:", avg_vector)

    for i in range(run_params['D']):
        resolution = 10.0 * np.spacing(avg_vector[i])  # 计算分辨率

        if np.any(diff_vectors[:, i] < resolution):
            if run_params['verbose'] >= 3:
                print(f"    WARNING: at least one vector within allowed resolution for dimension {i}")

            res_pt_count = run_params['DE']['NP'] // 4
            if np.count_nonzero(diff_vectors[:, i] < resolution) >= res_pt_count:
                is_converged = True
                if run_params['verbose'] >= 3:
                    print(f"WARNING: Points along dimension {i} cannot be resolved further. Ending civilization.")
        else:
            if run_params['verbose'] >= 3:
                print("    Population resolution okay.")
