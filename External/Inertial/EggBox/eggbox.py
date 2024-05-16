#!/usr/bin/env python3 
# Author: Pengxuan Zhu 
# Email: zhupx99@icloud.com

import json 
import time 
from numpy import sin, cos 
# read input from inp.dat
with open("input.json", 'r') as f1: 
    data = json.loads(f1.read())

from random import random

f = {
    "z" :   sin(data["xx"]) * cos(data["yy"]),
    "Time": random()
}


# time.sleep(f['Time'] * 0.1)

print("TestFunction: input is \n{}\nOutput is \n{}".format(data, f))

with open("output.json",'w') as f1:
    json.dump(f, f1) 