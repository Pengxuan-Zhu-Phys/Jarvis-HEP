#!/usr/bin/env python3 
from calendar import day_name
from lib2to3.pgen2.token import RPAR
import logging
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
from mpmath.functions.functions import re
import numpy as np
import pandas as pd 
from numpy.lib.function_base import meshgrid
from pandas.core.series import Series
from sympy import sympify
from sympy.geometry import parabola
import time 
from Func_lib import decode_path_from_file
from random import randint
from multiprocessing import Manager, Lock
import json

pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(pwd, "Sampling"))) 

class SManager():
    """
    The interface of sampler and workers, to deal with the conflicts 
    when the things runs in the multi process manners.
    """
    def __init__(self) -> None:
        self.ruid = Manager().dict()
        self.pack = {}
        self.sinf = Manager().dict()
        # self.lock = Lock()
        self.info = {}
        self.path = {}
        self.logger = None

    def make_factory(self):
        from IOs.IOs import to_file_woblk
        # print(self.path)
        print(self.info['pack']['incl'])
        # sys.exit()

        for pkg in self.info['pack']['incl']:
            self.pack[pkg] = {
                "workers":   Manager().dict(),
                "status":   "Init",
                "lock":     "OPEN",
                "runweb":   os.path.join(os.path.dirname(self.info['pack']['include'][pkg]["run info"]), "runweb.json")
            }
        to_file_woblk(self.pack.copy(), self.path['prwb'], method='to_json')
        # print(self.pack)


