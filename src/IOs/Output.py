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

from IOs.IOs import IOfile

class OutputFile(IOfile):
    def read_variables(self):
        """读取并处理文件中的变量"""
        if not self.path or not os.path.exists(self.path):
            raise ValueError("File not found. Ensure the path is correct and file exists.")
        
        try:
            with open(self.path, 'r') as file:
                content = file.read()
                # 示例：根据SLHA文件类型和变量配置处理文件内容
                extracted_variables = {}
                for var in self.variables:
                    if var['method'] == "SLHA":
                        # 示例：从内容中提取SLHA块中的变量值
                        block_name = var['block']
                        # 这里需要实现具体的内容提取逻辑
                        extracted_variables[var['name']] = "extracted_value"  # 示例值
                return extracted_variables
        except Exception as e:
            self.logger.error(f"Error reading output file '{self.name}': {e}")
