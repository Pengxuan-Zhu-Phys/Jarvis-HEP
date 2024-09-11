#!/usr/bin/env python3 
from base import Base
from IOs.parameter import Parameter
import uuid
from pprint import pprint
import yaml 
# import logging
from loguru import logger
import subprocess
import threading
import os
from time import sleep
from Module.module import Module
import asyncio
import sympy as sp 
import json
class CalculatorModule(Module):
    def __init__(self, name, config):
        super().__init__(name)
        if config["modes"]:
            self.modes = True
            self.analyze_config_multi()
        else:
            self.config             = config 
            self.type               = "Calculator"
            self.required_modules   = config["required_modules"]
            self.clone_shadow       = config["clone_shadow"]
            self.installation       = config["installation"]
            self.initialization     = config["initialization"]
            self.execution          = config["execution"]
            self.input              = config["execution"].get("input", [])
            self.output             = config["execution"].get("output", [])
            self.basepath           = config['path']
            self.handlers           = {}
            self.modes              = False
            self.is_installed       = False
            self.installation_event = None
            self.is_busy            = False
            self.PackID             = None
            self.run_log_file       = None
            self.sample_info        = {}
            self._funcs             = {}
            self.analyze_config()

    def assign_ID(self, PackID):
        self.PackID = PackID 

    def analyze_config(self):
        # get the variables for inputs and outputs, prepare for the workflow chart 
        from inner_func import update_funcs, update_const
        for ipf in self.input:
            if ipf['type'] == "SLHA":
                ipf['variables'] = {}
                for act in ipf['actions']:
                    if act['type'] == "Replace":
                        for var in act['variables']:
                            if "expression" in var.keys():
                                expr = sp.sympify(var['expression'], locals=update_funcs(update_const({})))
                                varis = set(expr.free_symbols)
                                for vv in varis:
                                    self.inputs[str(vv)] = None
                                    # ipf['variables'][str(vv)] = var 
                                var.update({"inc": varis})
                                ipf['variables'][var['name']] = var 
                            else:
                                self.inputs[var['name']] = None
                                ipf['variables'][var['name']] = var
                    elif act['type'] == "SLHA":
                        for var in act['variables']:
                            if "expression" in var.keys():
                                expr = sp.sympify(var['expression'], locals=update_funcs(update_const({})))
                                varis = set(expr.free_symbols)
                                for vv in varis:
                                    self.inputs[str(vv)] = None
                                    # ipf['variables'][str(vv)] = var 
                                var.update({"inc": varis})
                                ipf['variables'][var['name']] = var 
                            else:
                                self.inputs[var['name']] = None
                                ipf['variables'][var['name']] = var
                    elif act['type'] == "File":
                        self.inputs[act['source']] = None
                        ipf['variables'][act['source']] = {"name": act['source']}
                # print(ipf['variables'])

            elif ipf['type'] == "Json":
                ipf['variables'] = {}
                for act in ipf['actions']:
                    if act['type'] == "Dump":
                        for var in act['variables']:
                            if "expression" in var.keys():
                                expr = sp.sympify(var['expression'], locals=update_funcs(update_const({})))
                                varis = set(expr.free_symbols)
                                for vv in varis:
                                    self.inputs[str(vv)] = None
                                    # ipf['variables'][str(vv)] = var 
                                var.update({"inc": varis})
                                ipf['variables'][var['name']] = var 
                            # print(var.keys())
                            else:
                                ipf['variables'][var['name']] = var 
                                self.inputs[var['name']] = None
                # pprint(ipf['variables'])

            else:
                for ipv in ipf['variables']:
                    self.inputs[ipv['name']] = None
        # pprint(self.inputs)
        for opf in self.output:
            # self.outputs[opf['name']] = None
            for opv in opf['variables']:
                self.outputs[opv['name']] = None

    @property
    def funcs(self):
        return self._funcs

    def analyze_config_multi(self):
        pass 

    def add_handler(self, name, logpath, level, formatter):
        file_handler = logging.FileHandler(logpath, mode="a")
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)

        self.handlers[name] = file_handler

    def remove_handler(self, name):
        self.handlers[name].close()
        del self.handlers[name]

    def custom_format(record):
        module = record["extra"].get("module", "No module")
        if "raw" in record["extra"]:
            return "{message}"
        else:
            return f"\n·•· <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{record['message']}</level>"
    
    def create_basic_logger(self):
        logger_name = f"{self.name}-{self.PackID}"
        def filte_func(record):
            return record['extra']['module'] == logger_name
        
        self.basepath = self.decode_shadow_path(self.basepath)
        if not os.path.exists(self.basepath):
            os.makedirs(self.basepath)

        install_file_log = os.path.join(self.basepath, f"Installation_{self.name}-{self.PackID}.log")
        self.logger = logger.bind(module=logger_name, to_console=True, Jarvis=True)
        install_handler = self.logger.add(install_file_log, format=CalculatorModule.custom_format, level="DEBUG", rotation=None, retention=None, filter=filte_func )
        self.handlers['install'] = install_handler

    def update_sample_logger(self, sample_info):
        logger_name = f"Sample@{sample_info['uuid']} ({self.name}-No.{self.PackID})"
        def filte_func(record):
            return record['extra']['module'] == logger_name
        
        if not os.path.exists(sample_info['save_dir']):
            os.makedirs(sample_info['save_dir'])
        
        self.logger = logger.bind(module=logger_name, to_console=True, Jarvis=True)
        sample_handler = self.logger.add(sample_info['run_log'], format=CalculatorModule.custom_format, level="DEBUG", rotation=None, retention=None, filter=filte_func)
        self.handlers['sample'] = sample_handler
        self.logger.info("Sample created into the Disk")

    async def install(self):
        self.create_basic_logger()
        self.logger.warning(f"Start install {self.name}-{self.PackID}")

        for cmd in self.installation:
            if self.clone_shadow:
                command = self.decode_shadow_commands(cmd)
                self.logger.info(f" Run command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> ")
                # asyncio.run(self.run_command(command=command))
                await self.run_command(command=command)
            else: 
                asyncio.run(self.run_command(command=cmd))
                await self.run_command(command=cmd)
        
        logger.remove(self.handlers['install'])
        self.logger = None
        self.is_installed = True

    async def run_command(self, command):
        process = await asyncio.create_subprocess_shell(
            command['cmd'],
            stdout=asyncio.subprocess.PIPE, 
            stderr=asyncio.subprocess.PIPE, 
            cwd=command['cwd']
        )
        # stdout, stderr = await asyncio.gather(
        #     self.log_stream_info(process.stdout),
        #     self.log_stream_error(process.stderr)
        # )

        # await process.wait()
        try:
            stdout, stderr = await asyncio.gather(
                self.log_stream_info(process.stdout),
                self.log_stream_error(process.stderr)
            )
            await process.wait()
        finally:
            process.stdout.close()
            process.stderr.close()


    def log_stream(self, stream, level, logger):
        # Using the child logger to handle the logging 
        for line in iter(stream.readline, ''):
            logger.log(level, "\t{}".format(line.strip()))

    async def log_stream_info(self, stream):
        async for line in stream:
            self.logger.bind(raw=True).info(f"\t{line.decode()}")

    async def log_stream_error(self, stream):
        async for line in stream:
            self.logger.bind(raw=True).error(f"\t{line.decode()}")

    def initialize(self):
        for command in self.initialization:
            if self.clone_shadow: 
                command = self.decode_shadow_commands(command)
                self.logger.info(f" Run initialize command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> \n")
                asyncio.run(self.run_command(command=command))

    def execute(self, input_data, sample_info):
        self.sample_info = sample_info
        self.update_sample_logger(sample_info)

        self.initialize()
        # self.logger.info(f"Executing sample {sample_info['uuid']} in {self.name} on inputs: {json.dumps(input_data)}")
        
        result = {}
        input_obs = asyncio.run(self.load_input(input_data=input_data))

        if isinstance(input_obs, dict):
            result.update(input_obs)
        
        for command in self.execution['commands']:
            if self.clone_shadow:
                command = self.decode_shadow_commands(command)
                self.logger.info(f" Run execution command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> ")
                asyncio.run(self.run_command(command=command))

        output_obs = asyncio.run(self.read_output())
        if isinstance(output_obs, dict):
            result.update(output_obs)

        logger.remove(self.handlers['sample'])
        self.logger = None
        return result

    async def read_output(self):
        from IOs.IOs import IOfile
        read_coroutines = [
            IOfile.load(
                ffile['name'],
                path=ffile['path'],
                file_type=ffile['type'],
                variables=ffile['variables'],
                save=ffile['save'],
                logger=self.logger,
                PackID=self.PackID,
                sample_save_dir=self.sample_info['save_dir'],
                module=self.name,
                funcs=self.funcs
            ).read()
            for ffile in self.output
        ]
        observables = await asyncio.gather(*read_coroutines) 
        # print("Line 261 ->", observables)
        try:
            merged_observables = {key: val for d in observables for key, val in d.items()}
        except:
            merged_observables = {}
        # print("Line 262->", merged_observables)
        return merged_observables

    async def load_input(self, input_data):
        """
            Asynchronously loads input data into SLHA files based on the specified configuration.

            This method reads the configuration for each file from the `self.input` list, creates 
            instances for handling the files, and then concurrently writes the input data to these 
            files using their respective `write` methods. The operation is performed asynchronously 
            to improve performance when dealing with I/O operations and multiple files.

            Args:
            input_data (dict): The input data to be written into the files. This dictionary should 
                               contain the necessary information that matches the expected structure 
                               for each file type being written.

            The method uses `asyncio.gather` to concurrently execute all write operations for the 
            files defined in `self.input`. Each file is handled based on its configuration, including 
            the path, type, actions (variables), and whether it should be saved, along with other 
            metadata like `PackID`, the directory to save the file (`sample_save_dir`), and the 
            module name (`self.name`). The logger is used for logging purposes, and it's passed to 
            each file handler instance for consistent logging throughout the operation.

            After all files have been processed and the input data written, a log message is generated 
            to indicate completion of the loading process.
        """
        from IOs.IOs import IOfile
        write_coroutines = [
            IOfile.create(
                ffile["name"], 
                path=ffile['path'],
                file_type=ffile["type"],
                variables=ffile["actions"],
                save=ffile['save'],
                logger=self.logger, 
                PackID=self.PackID,
                sample_save_dir=self.sample_info['save_dir'],
                module=self.name,
                funcs=self.funcs
            ).write(input_data)  # Return directly to the coroutine
            for ffile in self.input
        ]
        # Perform all write operations concurrently and wait for them all to complete
        observables = await asyncio.gather(*write_coroutines)
        merged_observables = {key: val for d in observables for key, val in d.items()}
        return merged_observables
        # self.logger.warn(self.funcs)

    def decode_shadow_commands(self, cmd):
        command = {
            "cmd":  cmd['cmd'].replace("@PackID", self.PackID),
            "cwd":  cmd['cwd'].replace("@PackID", self.PackID)
        }
        return command

    def decode_shadow_path(self, path):
        path = self.decode_path(path)
        if "@PackID" in path:
            path = path.replace("@PackID", self.PackID)
        print(path)
        return path

