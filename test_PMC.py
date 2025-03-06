#!/usr/bin/env python3 

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde

# ============ 1) 定义问题: Eggbox + 目标 z=0.4 ============

def eggbox_function(theta):
    """
    theta: shape (..., 2), 每行 [x, y]
    return: z = sin(x)*cos(y)
    """
    return np.sin(theta[..., 0]) * np.cos(theta[..., 1])

def prior_density(theta, bounds):
    """
    简单的均匀先验: 如果 (x,y) 落在 bounds 内, 返回 1/Area, 否则 0
    bounds: [(low_x, high_x), (low_y, high_y)]
    """
    (low_x, high_x), (low_y, high_y) = bounds
    in_bounds = (
        (theta[0] >= low_x) & (theta[0] <= high_x) &
        (theta[1] >= low_y) & (theta[1] <= high_y)
    )
    area = (high_x - low_x)*(high_y - low_y)
    return (1.0/area) if in_bounds else 0.0

def soft_likelihood(theta, target=0.4, sigma_lik=0.1, bounds=None):
    """
    软似然: L(theta) = exp( -(z-0.4)^2 / (2*sigma_lik^2) ) * prior(theta)
    如果超出先验边界, 则返回 0
    """
    # 先验检查
    if bounds is not None:
        pd = prior_density(theta, bounds)
        if pd == 0:
            return 0.0
    else:
        pd = 1.0

    z_val = eggbox_function(theta)
    # 计算高斯型惩罚, 距离 (z - target) 越近, 似然越大
    distance = z_val - target
    ll = np.exp(-0.5 * (distance**2) / (sigma_lik**2))

    return pd * ll

def proposal_pdf(theta_new, theta_old, Sigma):
    """
    q(theta_new | theta_old) = N(theta_new | mean=theta_old, cov=Sigma)
    返回其密度值
    """
    return multivariate_normal.pdf(theta_new, mean=theta_old, cov=Sigma)

def proposal_sample(theta_old, Sigma):
    """
    从 q(theta_new | theta_old) = N(theta_old, Sigma) 中采样一个新点
    """
    return np.random.multivariate_normal(mean=theta_old, cov=Sigma)

# ============ 2) 关键函数: 生成下一代样本 + 权重更新 ============

def generate_new_population(old_samples, old_weights, N, Sigma, pi_function, q_function):
    """
    使用公式:
      w_{i,t} = pi(theta_{i,t}) / sum_j [ w_{j,t-1} * q(theta_{i,t} | theta_{j,t-1}) ]
    来更新新一代的样本和权重
    """
    d = old_samples.shape[1]  # 维度
    new_samples = np.zeros((N, d))
    new_weights = np.zeros(N)

    # 用于离散抽样的累积权重
    cum_weights = np.cumsum(old_weights)

    for i in range(N):
        # 1) 在上一代样本中, 根据权重选出一个点 theta_j
        u = np.random.rand()
        j = np.searchsorted(cum_weights, u)
        theta_old = old_samples[j]

        # 2) 用核函数(高斯) 从 theta_old 出发采样新点
        theta_new = proposal_sample(theta_old, Sigma)

        # 3) 计算 pi(theta_new)
        pi_val = pi_function(theta_new)

        # 4) 计算分母: sum_j [ w_{j,t-1} * q(theta_new | theta_{j,t-1}) ]
        denom = 0.0
        for k in range(N):
            denom += old_weights[k] * q_function(theta_new, old_samples[k], Sigma)

        # 5) 得到 w_new = pi_val / denom
        w_new = 0.0 if denom == 0 else pi_val/denom

        new_samples[i] = theta_new
        new_weights[i] = w_new

    # 6) 归一化
    sum_w = np.sum(new_weights)
    if sum_w > 0:
        new_weights /= sum_w
    else:
        new_weights[:] = 1.0/N

    return new_samples, new_weights

# ============ 3) 主循环: ABC-PMC ============

def abc_pmc(
    N=500,            # 粒子数
    T=5,              # 迭代代数
    bounds=[(-2*np.pi,2*np.pi),(-2*np.pi,2*np.pi)],
    sigma_lik=0.02,    # soft_likelihood 的标准差
    sigma_prop=0.3  # proposal 高斯的标准差(对角)
):
    """
    在 [(-5π,5π),(-5π,5π)] 范围内定义均匀先验,
    pi(theta) = prior(theta)* exp(-(sin(x)*cos(y)-0.4)^2/(2*sigma_lik^2)),
    proposal = N(theta_old, sigma_prop^2 * I).
    使用 ABC-PMC 迭代 T 代, 返回最终样本与权重, 以及每代记录.
    """
    d = 2  # 维度

    # 1) 初始化: 从先验均匀采样 N 个点
    old_samples = np.zeros((N,d))
    for i in range(N):
        old_samples[i,0] = np.random.uniform(bounds[0][0], bounds[0][1])
        old_samples[i,1] = np.random.uniform(bounds[1][0], bounds[1][1])

    # 计算初始权重: w_i = pi(theta_i)
    old_weights = np.zeros(N)
    for i in range(N):
        old_weights[i] = soft_likelihood(old_samples[i], 0.4, sigma_lik, bounds)
    sum_w = np.sum(old_weights)
    if sum_w > 0:
        old_weights /= sum_w
    else:
        old_weights[:] = 1.0/N

    # 2) 定义 pi_function, q_function, Sigma
    def pi_function(theta):
        return soft_likelihood(theta, 0.4, sigma_lik, bounds)

    def q_function(theta_new, theta_old, Sigma):
        return proposal_pdf(theta_new, theta_old, Sigma)

    Sigma = (sigma_prop**2)*np.eye(d)

    # 保存每一代的 (samples, weights)
    all_generations = []
    all_generations.append((old_samples.copy(), old_weights.copy()))

    # 3) 迭代 T 代
    for t in range(1, T+1):
        new_samples, new_weights = generate_new_population(
            old_samples, old_weights, N, Sigma, pi_function, q_function
        )
        old_samples, old_weights = new_samples, new_weights

        all_generations.append((old_samples.copy(), old_weights.copy()))
        print(f"Iteration {t}: sum(weights) = {np.sum(old_weights):.5f}, max(weights) = {np.max(old_weights):.5f}")

    return old_samples, old_weights, all_generations

# ============ 4) 主函数: 可视化所有代的点 + 最终点 ============

if __name__ == "__main__":
    N = 400
    T = 20
    final_samples, final_weights, all_gens = abc_pmc(N=N, T=T)

    # 1) 先在一张图中绘制所有代 (generations) 的采样点
    plt.figure(figsize=(8,6))
    # 为不同代分配不同颜色
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_gens)))
    for i, (samples_i, weights_i) in enumerate(all_gens):
        plt.scatter(samples_i[:,0], samples_i[:,1],
                    s=10, color=colors[i], alpha=0.7, label=f"Gen {i}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("All generations sample points (Eggbox ~ 0.4)")
    # plt.legend()
    plt.show()

    # 2) 再单独绘制最终一代的结果, 用“(sin(x)*cos(y)-0.4)^2”或软似然值做颜色
    final_z = eggbox_function(final_samples)
    # 这里使用 soft_likelihood 作为颜色映射, 以体现“谁更符合 z=0.4”
    final_lik = np.array([soft_likelihood(th, 0.4, 0.1, [(-2*np.pi,2*np.pi),(-2*np.pi,2*np.pi)]) 
                          for th in final_samples])

    plt.figure(figsize=(8,6))
    sc = plt.scatter(final_samples[:,0], final_samples[:,1],
                     c=final_lik, cmap='viridis', s=10, alpha=0.7, label="Final Samples")
    plt.colorbar(sc, label="Soft-likelihood")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Final Generation Samples (Color by Soft-likelihood)")
    plt.legend()
    plt.show()