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

from IOs.IOs import InputFile, IOfile

class SLHAInputFile(InputFile):
    def write(self, param_values):
        """Writing the input variables using the SLHA format"""

        self.path = self.decode_path(self.path)
        self.logger.info(f"Loading {self.path}")
        if not self.path:
            self.logger.error(f"Input file {self.path} not found")
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

        self.logger.debug(f"Writing the {param_values} -> \n\t{self.path}")
        try:
            with open(self.path, 'w') as slha_file:
                for var in self.variables:
                    if var['method'] == "Replace":
                        value = param_values.get(var['name'], "MISSING_VALUE")
                        placeholder = var['placeholder']
                        # 在这里添加具体的SLHA格式处理逻辑
                        content = f"{placeholder} {value}\n"
                        slha_file.write(content)
        except Exception as e:
            self.logger.error(f"Error writing SLHA input file '{self.name}': {e}")

class JsonInputFile(InputFile):
    def write(self, param_values):
        """将参数值以Json格式写入文件"""
        if not self.path:
            raise ValueError("File path not generated. Call 'generate_path' first.")
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

        data_to_write = {}
        for var in self.variables:
            # 假设每个变量都需要写入Json文件，且param_values包含所有必需的值
            value = param_values.get(var['name'], "MISSING_VALUE")
            data_to_write[var['name']] = value

        try:
            with open(self.path, 'w') as json_file:
                json.dump(data_to_write, json_file, indent=4)
        except Exception as e:
            self.logger.error(f"Error writing Json input file '{self.name}': {e}")
