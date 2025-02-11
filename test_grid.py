#!/usr/bin/env python3 

import numpy as np

# import numpy as np
import matplotlib.pyplot as plt

def generate_random_points():
    """ 在单位圆上随机撒四个点，并计算它们的极角 """
    angles = np.sort(np.random.uniform(0, 2 * np.pi, 4))  # 生成 4 个随机角度并排序
    return angles

def check_points_in_half_circle(angles):
    """
    判断 4 个点是否在同一个半圆上：
    如果 4 个点中最大的间隔 >= π，则其余 3 个点落入某个半圆
    """
    angle0 = angles.min()
    angle_gaps0 = angles - angle0
    angle_gaps1 = 2 * np.pi - angles - angle0  
    # 计算角度间隔
    max_gap = np.max(angle_gaps)  # 找到最大间隔
    return max_gap <= np.pi  # 最大间隔 <= π，说明可以找到一个半圆包含 4 个点

def plot_points(angles, result):
    """ 绘制单位圆并标注 4 个点 """
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

    # 画出单位圆
    circle = plt.Circle((0, 0), 1, color='black', fill=False, linewidth=1.2)
    ax.add_patch(circle)

    # 计算点的位置并绘制
    x, y = np.cos(angles), np.sin(angles)
    ax.scatter(x, y, color='red', s=80, label="Random Points")

    # 标注点
    for i, (xi, yi) in enumerate(zip(x, y)):
        ax.text(xi, yi, f"P{i+1}", fontsize=12, ha='right', color='blue')

    # 标题
    ax.set_title(f"All Points in a Half-Circle: {'Yes' if result else 'No'}")

    plt.legend()
    plt.show()

# 运行程序
n = 0
for ii in range(1000):  
    angles = generate_random_points()
    result = check_points_in_half_circle(angles)
    if result: 
        n += 1
# plot_points(angles, result)

# 输出判断结果
print("Prob -> {}".format(n/1000.))
print("Angles (radians):", angles)
print("All points in a half-circle:", result)