#!/usr/bin/env python3
# Test Multi-mode Programs
# Used by the bin/test_bridson.ini
# Xuan, 2019.01.23
from os import path
from numpy import loadtxt, sin, cos
import time 
import json

pwd = path.abspath(path.dirname(__file__))
# read input from inp.dat
with open(path.join(pwd, "mode.json"), 'r') as f1:
    card = json.loads(f1.read())

def func(fm):
    if fm == "f1":
        return sin(data[0])**2 + cos(data[1])**2
    elif fm == "f2":
        return sin(data["x"])**2 - cos(data["y"])**2


if card['mode'] == "numpy":
    data = loadtxt(path.join(pwd, "test_numpy_input.dat"))
    from random import random

    f = {
        "Z" :   func(card['func']),
        "Time": random()
    }

    time.sleep(f['Time'])

    print("TestFunction: input is \n{}\nOutput is \n{}".format(data, f))

    with open("output.json",'w') as f1:
        json.dump(f, f1) 
elif card['mode'] == "json":
    with open(path.join(pwd, "test_json_input.json"), 'r') as f1:
        data = json.loads(f1.read())
    from random import random

    f = {
        "Z" :   func(card['func']),
        "Time": random()
    }

    time.sleep(f['Time'])
    print("TestFunction: input is \n{}\nOutput is \n{}".format(data, f))

    with open("output.json",'w') as f1:
        json.dump(f, f1) 