#!/usr/bin/env python3 

from base import Base

import psutil
import time
import threading
import os
class Monitor(Base):
    def __init__(self, interval=0.5):
        self.interval = interval  # Monitor interval 
        self.keep_running = True
        self.logger = None  

    def start(self):
        # 启动资源监视线程
        monitor_thread = threading.Thread(target=self.monitor)
        monitor_thread.daemon = True
        monitor_thread.start()

    def monitor(self):
        # 持续监视资源，直到 self.keep_running 为 False
        while self.keep_running:
            self.log_resources()
            time.sleep(self.interval)

    def log_resources(self):
        # 监控当前进程的资源使用情况
        process = psutil.Process(os.getpid())
        
        # CPU 和内存信息
        cpu_usage = process.cpu_percent(interval=None)  # 获取CPU使用率
        memory_info = process.memory_info()  # 获取内存使用情况
        
        # 文件描述符和I/O操作
        open_files = process.num_fds()  # 当前打开的文件描述符数量
        io_counters = process.io_counters()  # I/O操作计数

        # 打印或记录监视数据
        if self.logger:
            self.logger.info(f"CPU Usage: {cpu_usage}%")
            self.logger.info(f"Memory Usage: {memory_info.rss / 1024 / 1024:.2f} MB")
            self.logger.info(f"Open Files: {open_files}")
            self.logger.info(f"I/O Counters: {io_counters}")
        else:
            print(f"CPU Usage: {cpu_usage}%")
            print(f"Memory Usage: {memory_info.rss / 1024 / 1024:.2f} MB")
            print(f"Open Files: {open_files}")
            print(f"I/O Counters: {io_counters}")

    def stop(self):
        # 停止资源监控
        self.keep_running = False

# 在 Jarvis-HEP 中集成监视器
monitor = Monitor(interval=10)
monitor.start()

# 主任务执行代码
try:
    # 运行主并发任务
    pass
finally:
    monitor.stop()  # 停止监视器
