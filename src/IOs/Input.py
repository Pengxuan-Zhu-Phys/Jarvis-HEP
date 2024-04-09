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
from inner_func import update_const

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
                            expr = sp.sympify(var['expression'], locals=update_const(self.funcs))
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
                                    expr = sp.sympify(var['expression'], locals=update_const(self.funcs))
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
                                    expr = sp.sympify(var['expression'], locals=update_const(self.funcs))
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
    async def write(self, param_values):
        """
        Asynchronously updates and writes data to a specified JSON file based on actions and variable expressions defined in self.variables.

        This method handles:
        - Reading existing JSON data from the specified file path. If the file doesn't exist or contains invalid JSON, starts with an empty dictionary.
        - Processing actions for each variable defined in self.variables:
            * For variables without expressions, fetches values directly from param_values.
            * For variables with expressions, evaluates the expressions using provided parameters and updates the variable with the result.
            * If a variable specifies an 'entry', updates the JSON data at the specified nested path; otherwise, updates at the root level.
        - Asynchronously writing the updated JSON data back to the file, formatted for readability.

        Parameters:
        - param_values (dict): A dictionary where keys are names of parameters to be used in variable expressions, and values are their corresponding values.

        Returns:
        - observables (dict): A dictionary containing the results of evaluated expressions for variables designated with "Dump" action, useful for monitoring or further processing.

        The method logs key actions and errors for debugging and monitoring purposes, ensuring transparency in the update and write processes.
        """
        
        self.path = self.decode_path(self.path)
        self.logger.warning(f"Start writing the input file -> {self.path}")
        observables = {}

        try:
            async with aiofiles.open(self.path, 'r') as f1:
                content = await f1.read()
                data_to_write = json.loads(content)
        except FileNotFoundError:
            self.logger.error(f"File not found: {self.path}")
            data_to_write = {}
        except json.JSONDecodeError:
            self.logger.error(f"Error decoding JSON from file: {self.path}")
            data_to_write = {}

        for action in self.variables:
            if action['type'] == "Dump":
                for var in action['variables']:
                    if not "expression" in var:
                        value = param_values.get(var['name'], "MISSING_VALUE")
                    else: 
                        # self.logger.warning(f"{var} -> {param_values} <- {self.funcs}")
                        expr = sp.sympify(var['expression'], locals=update_const(self.funcs))
                        para = set(expr.free_symbols)
                        num_expr = lambdify([str(par) for par in expr.free_symbols], expr, modules=[self.funcs, "numpy"])
                        symbol_values_strs = {str(key): value for key, value in param_values.items() if key in {str(par) for par in para}}
                        value = num_expr(**symbol_values_strs)
                        observables[var['name']] = value
                        self.logger.info(f"{expr} -> {param_values} -> {value}")

                    if "entry" not in var: 
                        data_to_write.update({var['name']: value})
                    else: 
                        self.update_json_by_entry(data_to_write, var['entry'], value)

            self.logger.warning(data_to_write)

        try:
            # 将字典转换为格式化的字符串
            json_str = json.dumps(data_to_write, indent=4)
            # 异步写入文件
            async with aiofiles.open(self.path, 'w') as json_file:
                await json_file.write(json_str)
        except Exception as e:
            self.logger.error(f"Error writing Json input file '{self.name}': {e}")
        return observables


    def update_json_by_entry(self, json_dict, entry, new_value):
        parts = entry.split('.')
        current = json_dict
        for part in parts[:-1]: 
            if part not in current:
                current[part] = {} 
            current = current[part]

        current[parts[-1]] = new_value

