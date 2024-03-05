#!/usr/bin/env python3
import subprocess

def get_terminal_output(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    

    