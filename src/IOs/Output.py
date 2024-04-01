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
from IOs.IOs import OutputFile

class SLHAOutputFile(OutputFile):
    """
        A class designed to asynchronously read, process, and optionally save SLHA files based on specified observables.

        Inherits from OutputFile and utilizes asynchronous file operations to efficiently process SLHA files without blocking the event loop. This class is particularly useful in environments where responsiveness and non-blocking I/O operations are critical.

        Attributes:
            Inherits all attributes from the OutputFile class such as path, logger, and save flag.
            variables (list of dicts): Each dictionary specifies an observable to extract from the SLHA file. The keys in the dictionary include:
                - name: A string representing the unique identifier of the observable.
                - block: The SLHA block from which the value should be extracted.
                - entry: The entry (or entries) within the block identifying the specific value or values to be extracted. Can be an integer for single values or a list for complex entries such as decay branching ratios.
    """

    async def read(self):
        """
            Methods:
                Asynchronously reads the SLHA file specified by the path attribute, extracts the values of the observables defined in the variables list, and optionally saves a copy of the file.
                
                Processing Steps:
                1. Decodes the file path to resolve any placeholders.
                2. Asynchronously reads the entire SLHA file content into memory.
                3. Parses the SLHA content using the pyslha library to facilitate data extraction.
                4. Iterates through the specified variables, extracting values from the corresponding blocks and entries.
                5. If the save flag is True, saves a copy of the original SLHA file to a specified directory or a temporary location, appending the module's name to the file name for identification.
                6. Returns a dictionary containing the extracted observables, keyed by their names as specified in the variables configuration.
                
                Parameters:
                    None, uses class attributes for configuration.
                
                Returns:
                    observables (dict): A dictionary with keys corresponding to the names of the observables specified in the variables attribute, and values being the extracted data from the SLHA file.
                
                Raises:
                    Exception: If an error occurs during file reading or parsing, logs an error message with the file name.
                
                Usage:
                    This method should be called within an asynchronous context using the 'await' keyword to ensure non-blocking operation. The returned dictionary can be used for further analysis or reporting.

            Note:
                This class requires an asynchronous environment and the pyslha library for parsing SLHA files. Ensure aiofiles and pyslha are installed and properly configured in your environment.
        """

        self.path = self.decode_path(self.path)
        self.logger.warning(f"Start reading the output file -> {self.path}")

        observables = {}
        content = None
        source = None
        try: 
            async with aiofiles.open(self.path, 'r') as slha_file: 
                source = await slha_file.read()
            content = pyslha.readSLHA(source)

            for var in self.variables:
                # self.logger.warning(f"{var}, {type(var['entry'])}")
                if var['block'] == "DECAY":
                    if isinstance(var['entry'], int):
                        value = float(content.decays[var['entry']].__dict__["totalwidth"])
                        observables[var['name']] = value
                    elif isinstance(var['entry'], list): 
                        decays = content.decays[var['entry'][0]].__dict__["decays"]
                        value = 0. 
                        for decay in decays:
                            if set(var['entry'][1:]) == set(decay.ids):
                                value = decay.br 
                                break
                        observables[var['name']] = value 
                    else:
                        self.logger.error(f"Unsupport decay entry {var['entry']}")
                elif content.blocks[var['block']]:
                    if isinstance(var['entry'], int):
                        value = content.blocks[var['block']][var['entry']]
                        observables[var['name']] = value 
                    elif isinstance(var['entry'], list):
                        value = content.blocks[var['block']][var['entry']]
                        observables[var['name']] = value 
                    else: 
                        self.logger.error(f"Unsupported block entry {var['entry']}")
                else: 
                    self.logger.error(f"Unsupported SLHA read item {var}")

            if self.save:
                target = os.path.join(self.sample_save_dir, f"{os.path.basename(self.path)}@{self.module}")
                async with aiofiles.open(target, "w") as dst_file:
                    await dst_file.write(source)
                observables[self.name] = target
            else: 
                target_path = os.path.join(self.sample_save_dir, ".temp")
                if not os.path.exists(target_path):
                    os.makedirs(target_path)
                target = os.path.join(target_path, f"{os.path.basename(self.path)}@{self.module}")
                async with aiofiles.open(target, "w") as dst_file: 
                    await dst_file.write(source)
                observables[self.name] = target
            observables['TestFunc'] = self.funcs['XenonSD2019'](100)


            self.logger.warning(f"Finish reading the input file -> {self.path}")
            return observables
        except Exception as e: 
            self.logger.error(f"Error reading SLHA input file '{self.name}': {e}")


