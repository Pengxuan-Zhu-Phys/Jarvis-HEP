#!/usr/bin/env python3 

import curses
import asyncio
import time
import psutil
import os 
# from datetime import datetime
# import random
# import setproctitle 

class JarvisMonitor:
    def __init__(self, process_name):
        self.process_name = process_name
        self.pid = self.get_pid_by_name(process_name)
        
        self.current_view = "Summary"
        self.runtime_in_seconds = 0
        self.scroll_offset = 0
        self.horizontal_scroll = False
        self.current_index = 0
        self.start_time = time.time()
         

    def get_pid_by_name(self, name):
        """ Get PID by process name """
        for proc in psutil.process_iter(['cmdline', 'pid', 'name']):
            if proc.info['cmdline']:
            # if proc.info['pid'] == self.pid:
                if name in proc.info['cmdline']:
                # if proc.info['name'] == name:
                    # print(proc.info)
                    return proc.info['pid']
        return None

    def get_open_files(self):
        """ Get list of open files for the specified process and its children """
        open_files = []
        try:
            for process in [psutil.Process(self.pid)] + psutil.Process(self.pid).children(recursive=True):
                try:
                    for file in process.open_files():
                        open_files.append(file.path.replace(os.path.expanduser('~'), '~'))
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        return open_files

    def init_screen(self, stdscr):
        """ Initialize screen settings """
        curses.curs_set(0)  
        stdscr.nodelay(True) 
        stdscr.clear()  

    def init_colors(self):
        curses.start_color()
        curses.use_default_colors()  
        for i in range(1, 8):
            curses.init_pair(i, i, -1) 

    def draw_box(self, stdscr, start_y, start_x, end_y, end_x, title="", text_color_pair=3, border_color_pair=3):
        stdscr.attron(curses.color_pair(border_color_pair))
        # Draw corner
        stdscr.addch(start_y, start_x, curses.ACS_ULCORNER)
        stdscr.addch(start_y, end_x, curses.ACS_URCORNER)
        stdscr.addch(end_y, start_x, curses.ACS_LLCORNER)
        stdscr.addch(end_y, end_x, curses.ACS_LRCORNER) 
        # Draw top and bottom borders
        for x in range(start_x + 1, end_x):
            stdscr.addch(start_y, x, curses.ACS_HLINE)
            stdscr.addch(end_y, x, curses.ACS_HLINE)    
        # Draw left and right borders
        for y in range(start_y + 1, end_y):
            stdscr.addch(y, start_x, curses.ACS_VLINE)
            stdscr.addch(y, end_x, curses.ACS_VLINE)
        # Add title 
        if title:
            stdscr.addstr(start_y, start_x + 2, title)

        stdscr.attroff(curses.color_pair(border_color_pair))
        stdscr.attron(curses.color_pair(text_color_pair))
        stdscr.addstr(start_y + 1, start_x + 1, " " * (end_x - start_x - 1))
        stdscr.attroff(curses.color_pair(text_color_pair))

    def draw_progress_bar(self, stdscr, start_y, start_x, width, percentage):
        max_y, max_x = stdscr.getmaxyx()
        if start_x + width + 2 > max_x:
            width = max_x - start_x - 2

        stdscr.addstr(start_y, start_x, '[')
        stdscr.addstr(start_y, start_x + width + 1, ']')

        fill_width = int(width * percentage / 100)

        stdscr.attron(curses.color_pair(2))
        for x in range(fill_width):
            if start_x + 1 + x < max_x:
                stdscr.addstr(start_y, start_x + 1 + x, '█')
        stdscr.attroff(curses.color_pair(2))

        stdscr.attron(curses.color_pair(3))
        for x in range(fill_width, width):
            if start_x + 1 + x < max_x:
                stdscr.addstr(start_y, start_x + 1 + x, '.')
        stdscr.attroff(curses.color_pair(3))

    def draw_summary(self, stdscr, start_y, start_x, cpu, mem, files, procs):
        """ Draw summary information including progress bars and numbers """
        self.draw_progress_bar(stdscr, start_y, start_x + 13, 61, cpu / psutil.cpu_count())
        stdscr.addstr(start_y, 2, "{}  : {:.1f}%".format("CPU", cpu))

        # Determine whether to show memory in MB or GB
        if mem < 1024 * 1024 * 1024:
            mem_display = mem / 1024 / 1024
            mem_unit = "MB"
        else:
            mem_display = mem / 1024 / 1024 / 1024
            mem_unit = "GB"

        self.draw_progress_bar(stdscr, start_y + 1, start_x + 13, 61, mem * 100 / psutil.virtual_memory().total)
        stdscr.addstr(start_y + 1, 2, "{}  : {} {}".format("MEM", int(mem_display), mem_unit))

        stdscr.addstr(start_y + 2, start_x, f"Files: {files}")
        stdscr.addstr(start_y + 3, start_x, f"Procs: {procs}")

    def draw_file_box(self, stdscr, start_y, start_x, height, width, content_list, scroll_offset=0):
        """ Draw a box with dynamic content """
        max_content = height  # Calculate the maximum number of files that can be displayed
        for i in range(max_content):
            idx = self.current_index + i
            if idx < len(content_list):
                file_path = content_list[idx]
                # Scroll the display if the path length exceeds the width
                if len(file_path) > (width - 10):
                    file_path = file_path[scroll_offset:scroll_offset + (width - 10)]
                line = f"{idx + 1:3}: {file_path.replace(os.path.expanduser('~'), '~')}"  # Format the number and path
                stdscr.addstr(start_y + 1 + i, start_x + 1, line[:width - 2])  # Prevent the string from exceeding the box width

    async def draw_subprocess_box(self, stdscr, start_y, start_x, height, width):
        """ Draw a box with subprocess information """
        procs = []
        cpu_usage = 0.
        mem_total = 0.

        for child in psutil.Process(self.pid).children(recursive=True) + [psutil.Process(self.pid)]:
            try:
                thread_count = child.num_threads()
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                thread_count = 'N/A'
            try:
                cpu_percent = await self.loop.run_in_executor(None, child.cpu_percent, 0.01)
                
                cpu_usage += cpu_percent if cpu_percent is not None else 0
                cpu_percent_str = f"{cpu_percent if cpu_percent is not None else 0}%"
                if child.memory_info():
                    mem_total += child.memory_info().rss
                    mem_usage = f"{child.memory_info().rss / 1024 / 1024:.1f} MB"
                else:
                    mem_usage = "N/A"
                name = child.name()
                procs.append((child.pid, cpu_percent_str, mem_usage, name, thread_count))
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                continue

        procs.sort(key=lambda x: float(x[1][:-1]), reverse=True)

        max_procs = height
        stdscr.attron(curses.color_pair(3))
        stdscr.addstr(start_y, start_x + 1, f"{'PID':<8} {'CPU':<7} {'MEM':<10} {'THREADS':<8} {'COMMAND':<20}")
        stdscr.attroff(curses.color_pair(3))

        for i, proc in enumerate(procs[:max_procs]):
            line = f"{proc[0]:<8} {proc[1]:<7} {proc[2]:<10} {proc[4]:<8} {proc[3]:<20}"
            stdscr.addstr(start_y + 1 + i, start_x + 1, line[:width - 2])

        return cpu_usage, mem_total

    def format_runtime(self, seconds):
        days = seconds // (24 * 3600)
        seconds = seconds % (24 * 3600)
        hours = seconds // 3600
        seconds %= 3600
        minutes = seconds // 60
        seconds %= 60

        return f"{days:02d}/{hours:02d}:{minutes:02d}:{seconds:02d}"

    async def main(self, stdscr):
        self.loop = asyncio.get_running_loop()
        self.loop = asyncio.get_running_loop()
        self.init_screen(stdscr)
        self.init_colors()

        files = self.get_open_files()
        try:
            while True:
                stdscr.clear()
                title = " Jarvis-HEP Monitor System - {} ".format(self.current_view)
                height, width = 24, 80

                if self.current_view == "Summary":
                    self.draw_box(stdscr, 0, 0, 23, width - 2, title, 1, 1)
                    footer = "[s] Summary  [f] Files  [p] Proc  [q] Quit"
                    stdscr.addstr(24, 2, footer)

                    files_count = len(files)
                    procs_count = len(psutil.Process(self.pid).children(recursive=True))

                    stdscr.attron(curses.color_pair(2))
                    stdscr.addstr(7, 2, "> Opened File List (Short)")
                    stdscr.attroff(curses.color_pair(2))
                    stdscr.attron(curses.color_pair(3))
                    stdscr.addstr(8, 2, "  No. Path")
                    stdscr.attroff(curses.color_pair(3))
                    fbh = 5
                    self.draw_file_box(stdscr, 8, 2, fbh, 80, files, scroll_offset=self.scroll_offset)

                    # Horizontal scroll
                    if self.horizontal_scroll:
                        self.scroll_offset += 10
                        if self.scroll_offset > max(len(files[self.current_index + i]) for i in range(fbh) if self.current_index + i < len(files)) - (width - 10):
                            self.scroll_offset = 0
                            self.horizontal_scroll = False
                    else:
                        self.current_index += fbh
                        if self.current_index >= len(files):
                            self.current_index = 0
                        self.horizontal_scroll = True

                    stdscr.attron(curses.color_pair(2))
                    stdscr.addstr(15, 2, "> Subprocess List (Short)")
                    stdscr.attroff(curses.color_pair(2))

                    cpu_usage, mem_usage = await self.draw_subprocess_box(stdscr, 16, 2, 5, width - 4)
                    self.draw_summary(stdscr, 2, 2, cpu_usage, mem_usage, files_count, procs_count)

                elif self.current_view == "Files":
                    self.draw_box(stdscr, 0, 0, 23, width - 2, title, 1, 1)
                    stdscr.addstr(24, 2, footer)
                    ttt = "> Opened Files List: \t{} in total".format(len(files))
                    stdscr.attron(curses.color_pair(2))
                    stdscr.addstr(1, 2, ttt)
                    stdscr.attroff(curses.color_pair(2))
                    self.draw_file_box(stdscr, 1, 2, 20, 80, files, scroll_offset=self.scroll_offset)

                    # Horizontal scroll
                    if self.horizontal_scroll:
                        self.scroll_offset += 10
                        if self.scroll_offset > max(len(files[self.current_index + i]) for i in range(20) if self.current_index + i < len(files)) - (width - 10):
                            self.scroll_offset = 0
                            self.horizontal_scroll = False
                    else:
                        if len(files) > 20:
                            self.current_index += 10
                            if self.current_index >= len(files):
                                self.current_index = 0
                            self.horizontal_scroll = True
                    stdscr.refresh()

                elif self.current_view == "Subprocess":
                    formatted_runtime = self.format_runtime(int(time.time() - psutil.Process(self.pid).create_time()))
                    self.draw_box(stdscr, 0, 0, 23, width - 2, title, 1, 1)
                    stdscr.addstr(24, 2, footer)
                    cpu_usage, mem = await self.draw_subprocess_box(stdscr, 2, 2, 20, width - 4)
                    if mem < 1024 * 1024 * 1024:
                        mem_display = mem / 1024 / 1024
                        mem_unit = "MB"
                    else:
                        mem_display = mem / 1024 / 1024 / 1024
                        mem_unit = "GB"
                    ttt = "Processes {}: {} total; CPU: {:.1f}%; Memory: {:.1f} {}   {}".format(self.pid, len(psutil.Process(self.pid).children(recursive=True)), cpu_usage, mem_display, mem_unit, formatted_runtime)
                    stdscr.attron(curses.color_pair(2))
                    stdscr.addstr(1, 2, ttt)
                    stdscr.attroff(curses.color_pair(2))

                key = stdscr.getch()
                if key == ord('q'):
                    break
                elif key == ord('s'):
                    self.current_view = "Summary"
                elif key == ord('f'):
                    self.current_view = "Files"
                elif key == ord('p'):
                    self.current_view = "Subprocess"
                elif key == -1:
                    curses.napms(1000)
        finally:
            stdscr.nodelay(False)

if __name__ == "__main__":
    # setproctitle.setproctitle("Jarvis-HEP")
    monitor = JarvisMonitor("Safari")
    asyncio.run(curses.wrapper(monitor.main))
