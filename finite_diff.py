#!/usr/bin/env python3 

import matplotlib.pyplot as plt
import time

# 数据列表，每个图形的数据
data_list = [
    ([1, 2, 3, 4, 5], [2, 3, 5, 7, 11], 'Plot 1'),
    ([1, 2, 3, 4, 5], [1, 4, 9, 16, 25], 'Plot 2'),
    ([1, 2, 3, 4, 5], [5, 4, 3, 2, 1], 'Plot 3')
]

for x, y, title in data_list:
    # 创建一个图形
    plt.plot(x, y)
    plt.title(title)
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')

    # 显示图形 (非阻塞)
    # plt.show(block=False)
    plt.show()
    time.sleep(1)


    # 关闭当前图形
    plt.close()

# 程序结束时关闭所有图形窗口
plt.close('all')
