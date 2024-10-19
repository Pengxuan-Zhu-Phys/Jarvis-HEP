#!/usr/bin/env python3 

import curses
import time
import psutil
import os 
from datetime import datetime
import random


def monitor_resources(stdscr, monitor_file):
    # 隐藏光标
    curses.curs_set(0)

    # 清屏
    stdscr.clear()

    try:
        # read the pid of the process 
        print("Just a inform!")
        with open(monitor_file, 'r') as f:
            pid = int(f.read().strip())
        pid = psutil.Process(os.getpid())
        print(pid)
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




def draw_section(stdscr, start_y, start_x, height, width, title, content):
    # Draw box with title
    draw_box(stdscr, start_y, start_x, height, width, title)
    
    # Populate box content
    y_offset = start_y + 1
    for line in content:
        stdscr.addstr(y_offset, start_x + 2, line)
        y_offset += 1




def init_colors():
    curses.start_color()
    curses.use_default_colors()  
    for i in range(1, 8):
        curses.init_pair(i, i, -1) 

def draw_progress_bar(stdscr, start_y, start_x, width, percentage):
    stdscr.addstr(start_y, start_x, '[')
    stdscr.addstr(start_y, start_x + width + 1, ']')

    fill_width = int(width * percentage / 100)

    stdscr.attron(curses.color_pair(2))
    for x in range(fill_width):
        stdscr.addstr(start_y, start_x + 1 + x, '█')
    stdscr.attroff(curses.color_pair(2))

    stdscr.attron(curses.color_pair(3))
    for x in range(fill_width, width):
        stdscr.addstr(start_y, start_x + 1 + x, '.')
    stdscr.attroff(curses.color_pair(3))



def init_screen(stdscr):
    """ 初始化屏幕设置 """
    curses.curs_set(0)  
    stdscr.nodelay(True) 
    stdscr.clear()  


def draw_box(stdscr, start_y, start_x, end_y, end_x, title="", text_color_pair=3, border_color_pair=3):
    stdscr.attron(curses.color_pair(border_color_pair))
    # Draw corner
    stdscr.addch(start_y, start_x, curses.ACS_ULCORNER)
    stdscr.addch(start_y, end_x, curses.ACS_URCORNER)
    stdscr.addch(end_y, start_x, curses.ACS_LLCORNER)
    stdscr.addch(end_y, end_x, curses.ACS_LRCORNER) 
    # draw top and bottom boarders
    for x in range(start_x + 1, end_x):
        stdscr.addch(start_y, x, curses.ACS_HLINE)
        stdscr.addch(end_y, x, curses.ACS_HLINE)    
    # draw left and right boarders
    for y in range(start_y + 1, end_y):
        stdscr.addch(y, start_x, curses.ACS_VLINE)
        stdscr.addch(y, end_x, curses.ACS_VLINE)
    # add titles 
    if title:
        stdscr.addstr(start_y, start_x + 2, title)

    # add texts 
    stdscr.attroff(curses.color_pair(border_color_pair))
    stdscr.attron(curses.color_pair(text_color_pair))
    stdscr.addstr(start_y + 1, start_x + 1, " " * (end_x - start_x - 1))
    stdscr.attroff(curses.color_pair(text_color_pair))

def draw_summary(stdscr, start_y, start_x, cpu, mem, files, procs):
    """ 绘制汇总信息包括进度条和数字 """
    draw_progress_bar(stdscr, start_y, start_x+14, 60, cpu)
    stdscr.addstr(start_y, 2, "{}  : {:.1f}%".format("CPU", cpu))

    draw_progress_bar(stdscr, start_y + 1, start_x+14, 60, mem)
    stdscr.addstr(start_y+1, 2, "{}  : {:.1f} MB".format("MEM", mem))

    stdscr.addstr(start_y + 2, start_x, f"Files: {files}")
    stdscr.addstr(start_y + 3, start_x, f"Procs: {procs}")

def draw_file_box(stdscr, start_y, start_x, height, width, content_list, current_index):
    """ 绘制带有动态内容的 box """

    max_content = height  # 计算可显示的最大文件数
    # max_content = height - 2  # 计算可显示的最大文件数
    for i in range(max_content):
        idx = current_index + i
        if idx < len(content_list):
            file_path = content_list[idx]
            line = f"{idx+1:3}: {file_path}"  # 格式化编号和路径
            stdscr.addstr(start_y + 1 + i, start_x + 1, line[:width-2])  # 防止字符串超出 box 宽度

def draw_subprocess_box(stdscr, start_y, start_x, height, width):
    """ 绘制子进程信息的盒子 """
    procs = []
    cpu_usage = 0.
    mem_total = 0. 

    for p in psutil.process_iter(['pid', 'cpu_percent', 'memory_info', 'name']):
        try:
            pid = p.info['pid']
            cpu_percent = p.info['cpu_percent']
            # 如果 cpu_percent 是 None 或无法转换为 float，设置为 0
            cpu_usage += cpu_percent if cpu_percent is not None else 0
            cpu_percent_str = f"{cpu_percent if cpu_percent is not None else 0}%"
                        # 如果 memory_info 为 None 则设置为 "N/A"
            if p.info['memory_info']:
                mem_total += p.info['memory_info'].rss / 1024 / 1024 / 1024
                mem_usage = f"{p.info['memory_info'].rss / 1024 / 1024:.1f} MB"  # 内存使用量（MB）
            else:
                mem_usage = "N/A"
            name = p.info['name']  # 进程名称
            procs.append((pid, cpu_percent_str, mem_usage, name))
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            # 忽略无法访问的进程
            continue

    procs.sort(key=lambda x: float(x[1][:-1]), reverse=True)

    max_procs = height
    # for i, proc in enumerate(procs[:max_procs]):
        # stdscr.addstr(start_y + 1 + i, start_x + 1, f"{proc[0]:<6} {proc[1]:<6} {proc[2]:<6}  {proc[3]}")
    stdscr.attron(curses.color_pair(3))
    stdscr.addstr(start_y, start_x + 1, f"{'PID':<8} {'CPU':<7} {'MEM':<10} {'COMMAND':<20}")
    stdscr.attroff(curses.color_pair(3))
    # files_count = 0


    for i, proc in enumerate(procs[:max_procs]):
        line = f"{proc[0]:<8} {proc[1]:<7} {proc[2]:<10} {proc[3]:<20}"
        stdscr.addstr(start_y + 1 + i, start_x + 1, line[:width - 2])
        # if proc[2] != "N/A":
        #     mem_usage += proc[2]
        # files_count += proc[2]
    return cpu_usage, mem_total

def format_runtime(seconds):
    days = seconds // (24 * 3600)
    seconds = seconds % (24 * 3600)
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60

    return f"{days:02d}/{hours:02d}:{minutes:02d}:{seconds:02d}"




def main(stdscr):
    init_screen(stdscr)
    init_colors()

    current_view = "Summary"

    files = [
        "/home/user/documents/file1.txt",
        "/home/user/documents/file2.txt",
        "/home/user/documents/file3.txt",
        "/var/log/syslog",
        "/etc/hosts",
        "/usr/local/bin/script.sh",
        "/usr/local/bin/script.sh",
        "/home/user/downloads/image1.png",
        "/home/user/downloads/image2.png",
        "/home/user/downloads/image3.png",
        "/home/user/downloads/image4.png",
        "/home/user/downloads/image5.png",
        "/home/user/downloads/image6.png",
        "/home/user/downloads/image7.png",
        "/home/user/downloads/image8.png",
        "/home/user/downloads/image9.png",
        "/home/user/downloads/image10.png",
        "/home/user/downloads/image11.png",
        "/home/user/downloads/image12.png",
        "/home/user/downloads/image13.png",
        "/home/user/downloads/image14.png",
        "/home/user/downloads/image15.png",
        "/home/user/downloads/image16.png",
        "/home/user/downloads/image17.png",
        "/home/user/downloads/image18.png",
        "/home/user/downloads/image19.png",
        "/home/user/downloads/image20.png",
        "/home/user/downloads/image21.png",
        "/home/user/downloads/image22.png",
        "/home/user/downloads/image23.png",
        "/home/user/downloads/image24.png",
        "/home/user/downloads/image25.png",
        "/home/user/downloads/image26.png",
        "/home/user/downloads/image27.png",
        "/home/user/downloads/image28.png",
        "/home/user/downloads/image29.png",
        "/home/user/downloads/image30.png",
        "/home/user/downloads/image31.png",
        "/home/user/downloads/image32.png",
        "/home/user/downloads/image33.png",
        "/home/user/downloads/image34.png",
        "/home/user/music/song.mp3"
    ]

    try:
        current_index = 0
        while True:
            process = psutil.Process(os.getpid())
            # process = psutil.Process(8221)
            stdscr.clear()
            title = " Jarvis-HEP Monitor System - {} ".format(current_view)
            height, width = 24, 80
            if current_view == "Summary":

                draw_box(stdscr, 0, 0, 23, width - 2, title, 1, 1)
                footer = "[s] Summary  [f] Files  [p] Proc  [q] Quit"

                stdscr.addstr(24, 2, footer)
                # 模拟数据更新
                # cpu_usage = random.randint(10, 90)
                cpu_usage = process.cpu_percent()
                # print(process.cpu_percent())
                # time.sleep(5)
                mem_usage = process.memory_info().rss / 1024 / 1024  # 以MB为单位

                # mem_usage = random.randint(10, 90)
                files_count = random.randint(10, 90)
                procs_count = random.randint(10, 90)

                stdscr.attron(curses.color_pair(2))
                stdscr.addstr(7, 2, "> Opened File List (Short)")
                stdscr.attroff(curses.color_pair(2))
                stdscr.attron(curses.color_pair(3))
                stdscr.addstr(8, 2, "  No. Path")
                stdscr.attroff(curses.color_pair(3))
                fbh = 5
                draw_file_box(stdscr, 8, 2, fbh, 80, files, current_index)
                current_index += fbh
                if current_index > len(files) + fbh:
                    current_index = 0  

                stdscr.attron(curses.color_pair(2))
                stdscr.addstr(15, 2, "> Subprocess List (Short)")
                stdscr.attroff(curses.color_pair(2))

                cpu_usage, mem_usage = draw_subprocess_box(stdscr, 16, 2, 5, width - 4)
                # stdscr.addstr(15, 50, "{}".format(cpu_usage))
                draw_summary(stdscr, 2, 2, cpu_usage, mem_usage, files_count, procs_count)

            elif current_view == "Files":
                draw_box(stdscr, 0, 0, 23, width - 2, title, 1, 1)
                stdscr.addstr(24, 2, footer)
                ttt = "> Opened Files List: \t{} in total".format(len(files))
                stdscr.attron(curses.color_pair(2))
                # stdscr.border(0)
                stdscr.addstr(1, 2, ttt)
                stdscr.attroff(curses.color_pair(2))
                # current_index = 0 
                draw_file_box(stdscr, 1, 2, 20, 80, files, current_index)

                current_index += 10
                if current_index > len(files) + 10:
                    current_index = 0  
                stdscr.refresh()

            elif current_view == "Subprocess":
                runtime_in_seconds = 90061  # 假设运行了90061秒
                formatted_runtime = format_runtime(runtime_in_seconds)
                draw_box(stdscr, 0, 0, 23, width-2, title, 1, 1)
                stdscr.addstr(24, 2, footer)
                ttt = "Processes: 50 total; CPU usage: 1180.0%; Memery usage: 1000 GB   {}".format(formatted_runtime)
                stdscr.attron(curses.color_pair(2))
                # stdscr.border(0)
                stdscr.addstr(1, 2, ttt)
                stdscr.attroff(curses.color_pair(2))
                draw_subprocess_box(stdscr, 2, 2, 20, width - 4)


            key = stdscr.getch()
            if key == ord('q'):  # Press 'q' to quit
                break
            elif key == ord('s'):  # Press 'b' to go back to summary
                current_view = "Summary"
            elif key == ord('c'):  # Show detailed CPU usage
                current_view = "CPU"
            elif key == ord('f'):  # Show detailed file information
                current_view = "Files"
            elif key == ord('p'):  # Show detailed subprocess information
                current_view = "Subprocess"
            elif key == -1:  # No input
                curses.napms(1000)  # Sleep briefly to reduce CPU usage
    finally:
        # Set nodelay back to False when exiting to clean up properly
        stdscr.nodelay(False)


    stdscr.nodelay(False)
    curses.endwin()


if __name__ == "__main__":
    monitor_file = '/path/to/your/monitor_file.txt'
    curses.wrapper(main)
