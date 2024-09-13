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
        self.log_file = None

    def start(self):
        # Start a monitor thread 
        monitor_thread = threading.Thread(target=self.monitor)
        monitor_thread.daemon = True
        monitor_thread.start()

    def monitor(self):
        # Keep monitoring until keep runing is False 
        while self.keep_running:
            self.log_resources()
            time.sleep(self.interval)

    def get_total_resources():
        # Get main process and all subprocesses 
        process = psutil.Process(os.getpid())
        children = process.children(recursive=True)  
        all_processes = [process] + children
    
        # Accumulate all source used 
        total_cpu = 0.0
        total_memory = 0
        total_open_files = 0
        total_io = 0
    
        for p in all_processes:
            try:
                total_cpu += p.cpu_percent(interval=None)
                total_memory += p.memory_info().rss  # memory in byte 
                total_open_files += p.num_fds()
                total_io += p.io_counters()
            except psutil.NoSuchProcess:
                # Skip the finished subprocess
                continue
            
        return (total_cpu, total_memory, total_open_files, total_io)

    def log_resources(self):
        usage = self.get_total_resources()
        if self.log_file is not None:
            
    
        if open_files:
            print(f"Currently open files (total {len(open_files)}):")
            for file in open_files:
                print(f"File: {file.path} - Mode: {file.mode}")
        else:
            print("No files are currently open.")



        if self.logger:
            self.logger.info(f"CPU Usage: {cpu_usage}%")
            self.logger.info(f"Memory Usage: {memory_info.rss / 1024 / 1024:.2f} MB")
            self.logger.info(f"Open Files: {open_files}")
            self.logger.info(f"I/O Counters: {io_counters}")

    def stop(self):
        # Stop monitor
        self.keep_running = False

# 在 Jarvis-HEP 中集成监视器
monitor = Monitor(interval=10)
monitor.start()

try:
    pass
finally:
    monitor.stop() 
