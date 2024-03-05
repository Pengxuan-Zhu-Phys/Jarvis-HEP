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
from sampling import SamplingVirtial
import json

class Birdson(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()

    def load_schema_file(self):
        self.schema = self.path['BirdsonSchema']

