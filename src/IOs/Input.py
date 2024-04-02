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
import asyncio 
import aiofiles
from sympy.utilities.lambdify import lambdify
import sympy as sp 
from IOs.IOs import InputFile

class SLHAInputFile(InputFile):
    async def write(self, param_values):
        """
        Asynchronously writes the input variables to the file in the SLHA format.
        
        This method supports different actions for updating the file, including:
        - Replace: Directly replaces placeholders in the file content with new values.
        - SLHA: Updates specific entries within SLHA blocks using pyslha library.
        - File: Copies the content from another source file to this file.
        
        Args:
            param_values (dict): A dictionary containing the new values for variables.
                                 The keys should match the variable names expected by the actions.
        
        The method first reads the entire content of the target file, then processes
        each action specified in `self.variables`. Depending on the action type, it updates
        the file content accordingly. After processing all actions, it writes the updated
        content back to the file. If `self.save` is True, it also saves a copy of the updated
        file in a specified directory with a modified name indicating the module name.
        
        Raises:
            Exception: If an error occurs during file reading or writing.
        """
        
        # Decode the file path to handle any special characters or template strings
        self.path = self.decode_path(self.path)
        self.logger.info(f"Start writing the input file -> {self.path}")
        
        # Ensure the file path is not empty
        if not self.path:
            self.logger.error(f"Input file {self.path} not found")
        
        # Create the directory for the file if it doesn't exist
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        observables = {}

        self.logger.debug(f"Writing the {param_values} -> \n\t{self.path}")
        try:
            content = ""
            async with aiofiles.open(self.path, 'r') as slha_file:
                content = await slha_file.read()

            for action in self.variables:
                # Direct replacement of placeholders with new values
                if action['type'] == "Replace":
                    for var in action['variables']:
                        if not "expression" in var:
                            value = param_values.get(var['name'], "MISSING_VALUE")
                        else: 
                            # self.logger.warning(f"{var} -> {param_values} <- {self.funcs}")
                            expr = sp.sympify(var['expression'], locals=self.funcs)
                            para = set(expr.free_symbols)
                            num_expr = lambdify([str(par) for par in expr.free_symbols], expr, modules=[self.funcs, "numpy"])
                            symbol_values_strs = {str(key): value for key, value in param_values.items() if key in {str(par) for par in para}}
                            value = num_expr(**symbol_values_strs)
                            observables[var['name']] = value
                            # self.logger.warning(f"{expr} -> {param_values} -> {value}")
                            # value = f"{float(value):.1E}"
                        placeholder = var['placeholder']
                        if value != "MISSING_VALUE":
                            value = f"{float(value):.8E}"
                        content = content.replace(placeholder, value)
                
                # Update SLHA blocks using pyslha
                elif action['type'] == "SLHA":
                    content = pyslha.readSLHA(content)
                    # self.logger.info(content.blocks)
                    for var in action['variables']:
                        if "block" in var.keys():
                            if isinstance(var['entry'], int):
                                if not "expression" in var:
                                    value = param_values.get(var['name'], "MISSING_VALUE")
                                else:
                                    # self.logger.warning(f"{var} -> {param_values} <- {self.funcs}")
                                    expr = sp.sympify(var['expression'], locals=self.funcs)
                                    para = set(expr.free_symbols)
                                    num_expr = lambdify([str(par) for par in expr.free_symbols], expr, modules=[self.funcs, "numpy"])
                                    symbol_values_strs = {str(key): value for key, value in param_values.items() if key in {str(par) for par in para}}
                                    value = num_expr(**symbol_values_strs)
                                    observables[var['name']] = value
                                    # self.logger.warning(f"{expr} -> {param_values} -> {value}")
                                    # value = f"{float(value):.1E}"
                                if value != "MISSING_VALUE":
                                    value = f"{float(value):.8E}"
                                content.blocks[var['block']][var['entry']] = value 
                            elif isinstance(var['entry'], tuple):
                                if not "expression" in var:
                                    value = param_values.get(var['name'], "MISSING_VALUE")
                                else:
                                    expr = sp.sympify(var['expression'], locals=self.funcs)
                                    para = set(expr.free_symbols)
                                    num_expr = lambdify([str(par) for par in expr.free_symbols], expr, modules=[self.funcs, "numpy"])
                                    symbol_values_strs = {str(key): value for key, value in param_values.items() if key in {str(par) for par in para}}
                                    value = num_expr(**symbol_values_strs)
                                    observables[var['name']] = value
                                if value != "MISSING_VALUE":
                                    value = f"{float(value):.8E}"
                                content.blocks[var['block']][tuple(var['entry'])] = value 
                            else: 
                                self.logger.warning(f"Invalid Entry type: {type(var['entry'])}. Entry must be an integer or a tuple of integers.")
                    content = pyslha.writeSLHA(content, ignorenobr=True)
                
                
                # Copy content from another file
                elif action['type'] == "File":
                    source_path = param_values.get(action['source'], None)
                    if os.path.exists(source_path):
                        async with aiofiles.open(source_path, 'r') as source_file: 
                            content = await source_file.read()
                        break
                    else: 
                        self.logger.warning(f"Input source file is not found: -> \n\t{action['source']} \n\t{source_path}")
            
            # Write the updated content back to the file
            async with aiofiles.open(self.path, 'w') as slha_file:
                await slha_file.write(content)
            
            # Save a copy of the file if required
            if self.save:
                target = os.path.join(self.sample_save_dir, f"{os.path.basename(self.path)}@{self.module}")
                async with aiofiles.open(target, "w") as dst_file:
                    await dst_file.write(content)
                observables[self.name] = target

            self.logger.info(f"Finish writing the input file -> {self.path}")
            return observables
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
