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
    def __init__(self, name, path, file_type, variables, save, logger, PackID):
        self.logger     = logger
        self.PackID     = PackID
        self.name       = name
        self.path_template = None
        self.file_type = file_type
        self.variables = variables
        self.save = save
        self.path = path  

    def decode_path(self, path) -> None:
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

        if self.PackID is not None:
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

    @classmethod
    def create(cls, name, path, file_type, variables, save, logger, PackID):
        print(file_type)
        if file_type == "Json":
            from IOs.Input import JsonInputFile
            logger.debug(f"Adding the file {name} as 'JsonInputFile' type")
            return JsonInputFile(name, path, file_type, variables, save, logger, PackID)
        elif file_type == "SLHA":
            from IOs.Input import SLHAInputFile
            logger.debug(f"Adding the file {name} as 'JsonInputFile' type")
            return SLHAInputFile(name, path, file_type, variables, save, logger, PackID)
        else:
            raise ValueError(f"Unsupported file type: {file_type}")

class InputFile(IOfile):
    def write(self, param_values):
        raise NotImplementedError("This method should be implemented by subclasses.")

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
