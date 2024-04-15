#!/usr/bin/env python3 

import psutil
import time

def monitor_cpu_usage(duration=10, interval=0.1):
    """ Monitor CPU usage percentage over a duration of time. """
    start_time = time.time()
    while time.time() - start_time < duration:
        cpu_percent = psutil.cpu_percent(interval=interval)
        print(f"CPU Usage: {cpu_percent}%")
        # time.sleep(interval - 0.01)  # adjust sleep if necessary

monitor_cpu_usage()
