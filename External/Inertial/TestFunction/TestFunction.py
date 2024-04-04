#!/usr/bin/env python3
# Test function
# Used by the bin/example.ini
# Yang, 2019.01.23
from numpy import loadtxt, sin, cos
import time 
import json

# read input from inp.dat
data = loadtxt(r"TestFunction_input.dat")
from random import random

f = {
    "Z" :   sin(data[0])**2 + cos(data[1])**2,
    "Time": random()
}


time.sleep(f['Time'])

print("TestFunction: input is \n{}\nOutput is \n{}".format(data, f))

with open("output.json",'w') as f1:
    json.dump(f, f1) 