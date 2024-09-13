#!/usr/bin/env python3 

import curses
import time
import psutil

def monitor_resources(stdscr, monitor_file):
    # 隐藏光标
    curses.curs_set(0)

    # 清屏
    stdscr.clear()

    try:
        # 读取PID
        with open(monitor_file, 'r') as f:
            pid = int(f.read().strip())

        # 获取进程
        process = psutil.Process(pid)

        while True:
            # 获取系统资源使用情况
            cpu_usage = process.cpu_percent(interval=1)
            memory_info = process.memory_info().rss / 1024 / 1024  # 以MB为单位
            open_files = process.num_fds()

            # 清空屏幕
            stdscr.clear()

            # 打印监控信息
            stdscr.addstr(0, 0, f"Monitoring PID: {pid}")
            stdscr.addstr(1, 0, f"CPU Usage: {cpu_usage:.2f}%")
            stdscr.addstr(2, 0, f"Memory Usage: {memory_info:.2f} MB")
            stdscr.addstr(3, 0, f"Open Files: {open_files}")
            stdscr.addstr(4, 0, "-----------------------------")
            stdscr.addstr(5, 0, "Press 'q' to quit.")

            # 刷新屏幕
            stdscr.refresh()

            # 检查是否按下了退出键
            if stdscr.getch() == ord('q'):
                break

            time.sleep(1)  # 每秒刷新一次

    except FileNotFoundError:
        stdscr.addstr(0, 0, "Monitor file not found.")
        stdscr.refresh()
        stdscr.getch()
    except psutil.NoSuchProcess:
        stdscr.addstr(0, 0, f"Process with PID {pid} not found.")
        stdscr.refresh()
        stdscr.getch()
    except KeyboardInterrupt:
        stdscr.addstr(0, 0, "Exiting monitor.")
        stdscr.refresh()
        stdscr.getch()

if __name__ == "__main__":
    monitor_file = '/path/to/your/monitor_file.txt'
    curses.wrapper(monitor_resources, monitor_file)
