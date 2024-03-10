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

class IOfile(Base):
    def __init__(self) -> None:
        self.filepath   = None
        self.vars       = None
        self.type       = None
        self.logger     = None
        self.PackID     = None 

    def decode_path(self, path, PackID=None) -> None:
        """
        Resolves special markers in the provided path.

        Parameters:
        - path: The path to be resolved.
        - jarvis_root: The root directory path of Jarvis.

        Returns:
        - The resolved full path.
        """
        # Replace the Jarvis root directory marker &J
        if "&J" in path:
            path = path.replace("&J", self.path['jpath'])

        # Replace the user home directory marker ~
        if "~" in path:
            path = os.path.expanduser(path)

        if PackID is not None and self.PackID is not None:
            path = path.replace("@PackID", self.PackID)
        
        return path  

    def set_logger(self, logger):
        self.logger = logger

    def set_packID(self, packID):
        self.PackID = packID

    def set_filepath(self, filepath):
        fpath = self.decode_path(fpath)
        if os.path.exists(fpath):
            self.filepath = fpath
        else:
            self.logger.error(f"File not found: {filepath}")
            raise FileNotFoundError

class InputFile(IOfile):
    def __init__(self) -> None:
        super().__init__()
        self.type = 'input'
    
    def load_config(self, config):
        print(config)


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

# class Parameters(Module):
#     def __init__(self, name, parameter_definitions):
#         super().__init__(name)
#         # 将参数字典转化为Parameter实例的列表
#         self.parameters = [Parameter(**param_def) for param_def in parameter_definitions]

#     def generate_parameters(self):
#         # 生成并返回所有参数的值
#         return {param.name: param.generate_value() for param in self.parameters}

# # 使用示例
# parameter_definitions = [
#     {'name': 'mn1', 'description': 'Variable X following a Flat distribution for sampling.', 'distribution': {'type': 'Flat', 'parameters': {'min': 2.0, 'max': 5.0, 'length': 60}}},
#     {'name': 'msmuR', 'description': 'Variable Y following a Logarithmic distribution for sampling.', 'distribution': {'type': 'Log', 'parameters': {'min': 0.01, 'max': 10.0, 'length': 50}}}
# ]

# parameters_module = Parameters("Sampling Parameters", parameter_definitions)
# generated_parameters = parameters_module.generate_parameters()
# print(generated_parameters)
