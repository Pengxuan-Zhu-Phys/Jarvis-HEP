#!/usr/bin/env python3

from copy import deepcopy
import imp
import logging
import os
from plistlib import FMT_XML
import sys
from re import L
import json
import numpy
import xslha
import pyslha
import xmltodict
import pandas as pd 
from jarvishep.base import Base
import asyncio
import aiofiles


class IOfile(Base):
    def __init__(self, name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs):
        self.logger     = logger
        self.PackID     = PackID
        self.name       = name
        self.path_template = None
        self.file_type = file_type
        self.variables = variables
        self.save = save
        self.path = path  
        self.sample_save_dir = sample_save_dir
        self.module = module
        self.funcs = funcs

    def decode_path(self, path) -> None:
        """
        Resolves special markers in the provided path.

        Parameters:
        - path: The path to be resolved.
        - jarvis_root: The root directory path of Jarvis.

        Returns:
        - The resolved full path.
        """
        if path is None or not isinstance(path, str):
            return path

        normalized = path.replace("\\", "/")
        if (normalized.startswith("&SRC/") or normalized.startswith("&J/")) and "/src/card/" in normalized:
            raise ValueError(
                "Legacy card path prefix is no longer supported: "
                f"{path}. Use '&SRC/card/...' in packaged runtime."
            )

        runtime_root = None
        if isinstance(getattr(self, "path", None), dict):
            runtime_root = self.path.get("task_root") or self.path.get("jpath")
        if not runtime_root:
            runtime_root = os.getenv("JARVIS_HEP_TASK_ROOT") or os.getenv("JHEP_TASK_ROOT") or os.getcwd()

        source_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

        if "&SRC" in path:
            path = path.replace("&SRC", source_root)
        if "&J" in path:
            path = path.replace("&J", runtime_root)

        # Replace the user home directory marker ~
        if "~" in path:
            path = os.path.expanduser(path)

        if "://" not in path and not os.path.isabs(path):
            path = os.path.abspath(os.path.join(runtime_root, path))

        if self.PackID is not None:
            path = path.replace("@PackID", self.PackID)
        
        return path  

    def set_logger(self, logger):
        self.logger = logger

    def set_packID(self, packID):
        self.PackID = packID

    def set_filepath(self, filepath):
        fpath = self.decode_path(filepath)
        if os.path.exists(fpath):
            self.filepath = fpath
        else:
            self.logger.error(f"File not found: {filepath}")
            raise FileNotFoundError

    @classmethod
    def create(cls, name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs):
        if file_type == "Json":
            from jarvishep.IOs.Input import JsonInputFile
            logger.debug(f"Adding the file {name} as 'JsonInputFile' type")
            return JsonInputFile(name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs)
        elif file_type == "SLHA":
            from jarvishep.IOs.Input import SLHAInputFile
            logger.debug(f"Adding the file {name} as 'SLHAInputFile' type")
            return SLHAInputFile(name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs)
        else:
            raise ValueError(f"Unsupported file type: {file_type}")

    @classmethod
    def load(cls, name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs):
        if file_type == "SLHA":
            from jarvishep.IOs.Output import SLHAOutputFile
            logger.debug(f"Loading the file {name} as 'SLHAOutputFile' type")
            return SLHAOutputFile(name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs)
        elif file_type == "Json":
            from jarvishep.IOs.Output import JsonOutputFile
            logger.debug(f"Loading the file {name} as 'JsonOutputFile' type")
            return JsonOutputFile(name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs)
        elif file_type == "xSLHA":
            from jarvishep.IOs.Output import xSLHAOutputFile
            logger.debug(f"Loading the file {name} as 'xSLHAOutputFile' type")
            return xSLHAOutputFile(name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs)
        elif file_type == "File":
            from jarvishep.IOs.Output import FileOutput
            logger.debug(f"Loading the file {name} as 'fileOutput' type")
            return FileOutput(name, path, file_type, variables, save, logger, PackID, sample_save_dir, module, funcs)

class InputFile(IOfile):
    async def write(self, param_values):
        raise NotImplementedError("This method should be implemented by subclasses.")



class OutputFile(IOfile):
    async def read(self):
        """Read the varaible values in the output files """
        raise NotImplementedError("This method should be implemented by subclasses.")




class Parameter:
    def __init__(self, name, description, distribution, ptype=None):
        self.name = name
        self.description    = description
        self.distribution   = distribution
        self.type           = "Param"
        if ptype is not None: 
            self.type           = ptype
             

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
