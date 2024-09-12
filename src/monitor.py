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

    def get_total_resources():
        # 获取主进程和所有子进程
        process = psutil.Process(os.getpid())
        children = process.children(recursive=True)  # 获取所有子进程，递归查找
        all_processes = [process] + children
    
        # 累加每个进程的资源使用
        total_cpu = 0.0
        total_memory = 0
        total_open_files = 0
    
        for p in all_processes:
            try:
                total_cpu += p.cpu_percent(interval=None)
                total_memory += p.memory_info().rss  # 以字节为单位
                total_open_files += p.num_fds()
            except psutil.NoSuchProcess:
                # 如果子进程已经终止，跳过该进程
                continue
            
        return total_cpu, total_memory, total_open_files

    def log_resources(self):
        # 获取主进程和所有子进程
        process = psutil.Process(os.getpid())
        # children = process.children(recursive=True)  # 获取所有子进程，递归查找
        # all_processes = [process] + children

        # # 累加每个进程的资源使用情况
        # total_cpu = 0.0
        # total_memory = 0
        # total_open_files = 0
        # total_io_counters = psutil._common.siocnt(0, 0, 0, 0)  # 初始化I/O计数器

        # for p in all_processes:
        #     try:
        #         total_cpu += p.cpu_percent(interval=None)
        #         total_memory += p.memory_info().rss  # 以字节为单位
        #         total_open_files += p.num_fds()

        #         # I/O 计数器 (读取与写入操作)
        #         io_counters = p.io_counters()
        #         total_io_counters = total_io_counters._replace(
        #             read_count=total_io_counters.read_count + io_counters.read_count,
        #             write_count=total_io_counters.write_count + io_counters.write_count,
        #             read_bytes=total_io_counters.read_bytes + io_counters.read_bytes,
        #             write_bytes=total_io_counters.write_bytes + io_counters.write_bytes,
        #         )
        #     except psutil.NoSuchProcess:
        #         # 如果子进程已经终止，跳过该进程
        #         continue

        # # 将资源使用情况记录到日志或输出到控制台
        # if self.logger:
        #     self.logger.info(f"Total CPU Usage: {total_cpu}%")
        #     self.logger.info(f"Total Memory Usage: {total_memory / 1024 / 1024:.2f} MB")
        #     self.logger.info(f"Total Open Files: {total_open_files}")
        #     self.logger.info(f"Total I/O Counters: {total_io_counters}")
        # else:
        #     print(f"Total CPU Usage: {total_cpu}%")
        #     print(f"Total Memory Usage: {total_memory / 1024 / 1024:.2f} MB")
        #     print(f"Total Open Files: {total_open_files}")
        #     print(f"Total I/O Counters: {total_io_counters}")


        # CPU 和内存信息
        cpu_usage = process.cpu_percent(interval=None)  # 获取CPU使用率
        memory_info = process.memory_info()  # 获取内存使用情况
        
        # 文件描述符和I/O操作
        open_files = process.num_fds()  # 当前打开的文件描述符数量
        io_counters = process.io_counters()  # I/O操作计数

        open_files = process.open_files()
    
        if open_files:
            print(f"Currently open files (total {len(open_files)}):")
            for file in open_files:
                print(f"File: {file.path} - Mode: {file.mode}")
        else:
            print("No files are currently open.")



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


    # def log_resources(self):
    #     # 监控当前进程的资源使用情况
    #     process = psutil.Process(os.getpid())
        
    #     # CPU 和内存信息
    #     cpu_usage = process.cpu_percent(interval=None)  # 获取CPU使用率
    #     memory_info = process.memory_info()  # 获取内存使用情况
        
    #     # 文件描述符和I/O操作
    #     open_files = process.num_fds()  # 当前打开的文件描述符数量
    #     io_counters = process.io_counters()  # I/O操作计数

    #     # 打印或记录监视数据
    #     if self.logger:
    #         self.logger.info(f"CPU Usage: {cpu_usage}%")
    #         self.logger.info(f"Memory Usage: {memory_info.rss / 1024 / 1024:.2f} MB")
    #         self.logger.info(f"Open Files: {open_files}")
    #         self.logger.info(f"I/O Counters: {io_counters}")
    #     else:
    #         print(f"CPU Usage: {cpu_usage}%")
    #         print(f"Memory Usage: {memory_info.rss / 1024 / 1024:.2f} MB")
    #         print(f"Open Files: {open_files}")
    #         print(f"I/O Counters: {io_counters}")

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
