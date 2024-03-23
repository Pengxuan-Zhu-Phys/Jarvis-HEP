#!/usr/bin/env python3

from copy import deepcopy
import imp
import logging
import os
from plistlib import FMT_XML
import sys
from re import L
import json
from matplotlib.pyplot import axis
import numpy
import xslha
import pyslha
import xmltodict
import pandas as pd 
from base import Base

class Parameter:
    def __init__(self, name, description, distribution):
        self.name = name
        self.description = description
        self.distribution = distribution

    def generate_value(self):
        # 根据self.distribution的类型和参数生成参数值
        # 示例代码，具体逻辑需根据分布类型实现
        if self.distribution['type'] == 'Flat':
            return (self.distribution['parameters']['min'] + self.distribution['parameters']['max']) / 2
        elif self.distribution['type'] == 'Log':
            # Log分布的值生成逻辑，这里仅为示例
            return (self.distribution['parameters']['min'] + self.distribution['parameters']['max']) / 2
        else:
            raise ValueError("Unsupported distribution type")
