#!/usr/bin/env python3 
from base import Base
from IOs.parameter import Parameter
import uuid
from pprint import pprint
import yaml 
import logging
import subprocess
import threading
import os
from time import sleep
from Module.module import Module

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
            self.modes = False
            self.is_installed = False
            self.is_busy = False
            self.PackID  = None
            self.analyze_config()
            self.run_log_file   = None

    def assign_ID(self, PackID):
        self.PackID = PackID 

    def analyze_config(self):
        # get the variables for inputs and outputs, prepare for the workflow chart 
        for ipf in self.input:
            # self.inputs[ipf['name']] = None
            for ipv in ipf['variables']:
                self.inputs[ipv['name']] = None
        
        for opf in self.output:
            # self.outputs[opf['name']] = None
            for opv in opf['variables']:
                self.outputs[opv['name']] = None


    def analyze_config_multi(self):
        pass 

    def create_child_logger(self, logger_name):
        child_installation_logger_name = f"{self.logger.name}.{logger_name}"
        self.child_logger = logging.getLogger(child_installation_logger_name)
        self.child_logger.propagate = False

        simple_formatter    = logging.Formatter('%(message)s')
        file_handler        = logging.FileHandler(self.path['run_log_file'], mode="a")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(simple_formatter)
        self.child_logger.addHandler(file_handler)


    def update_sample_logger(self, sample_info):
        logger_name = f"Sample@{sample_info['uuid']} <{self.name}-No.{self.PackID}>"
        if not os.path.exists(sample_info['save_dir']):
            os.makedirs(sample_info['save_dir'])
        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(logging.DEBUG)
        self.logger.propagate = False

        self.path['run_log_file'] = os.path.join(sample_info['save_dir'], "Sample_running.log")
        formatter = logging.Formatter("\n·•· %(name)s\n\t -> %(asctime)s - %(levelname)s >>> \n%(message)s")
        file_handler = logging.FileHandler(self.path['run_log_file'], mode='a')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)

        jlog_handler = logging.FileHandler(sample_info['jarvis_log'])
        jlog_handler.setFormatter(formatter)
        jlog_handler.setLevel(logging.WARNING)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.WARNING)
        console_handler.setFormatter(formatter)

        self.logger.addHandler(console_handler)
        self.logger.addHandler(file_handler)


        self.logger.warning("Sample created into the Disk")
        self.logger.warning(self.logger.handlers)


    def install(self):
        self.basepath = self.decode_shadow_path(self.basepath)
        if not os.path.exists(self.basepath):
            os.makedirs(self.basepath)

        self.logger.setLevel(logging.DEBUG) 
        install_file_log = os.path.join(self.basepath, f"Installation_{self.name}-{self.PackID}.log")
        formatter = logging.Formatter(" %(name)s - %(asctime)s - %(levelname)s >>> \n%(message)s")
        Afile_handler = logging.FileHandler(install_file_log, mode='a')
        Afile_handler.setLevel(0)
        Afile_handler.setFormatter(formatter)
        self.logger.addHandler(Afile_handler)
        self.logger.info(f"Installing {self.name}-{self.PackID} ....")
        self.create_child_logger("Installation")
        
        simple_formatter = logging.Formatter('%(message)s')
        file_handler = logging.FileHandler(install_file_log, mode="a")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(simple_formatter)
        self.child_logger.addHandler(file_handler)


        for cmd in self.installation:
            if self.clone_shadow:
                command = self.decode_shadow_commands(cmd)
                self.logger.info(f" Run command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> \n")
                self.run_command(command=command, child_logger=self.child_logger)
            else: 
                self.run_command(command=cmd, child_logger=self.child_logger)

        self.logger.removeHandler(Afile_handler)
        Afile_handler.close()
        self.child_logger = None

    def remove_running_file_log(self):
        # self.logger.removeHandler(self.run_log_handler)
        # self.run_log_handler.close()
        self.child_logger = None

    def run_command(self, command, child_logger):
        with subprocess.Popen(
                                command['cmd'], 
                                shell=True, 
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, 
                                text=True, 
                                cwd=command['cwd']) as process:
            stdout_thread = threading.Thread(target=self.log_stream, args=(process.stdout, logging.INFO, child_logger))
            stdout_thread.start()
            stderr_thread = threading.Thread(target=self.log_stream, args=(process.stderr, logging.ERROR, child_logger)) 
            stderr_thread.start()

            stdout_thread.join()
            stderr_thread.join()

            process.wait()

    def log_stream(self, stream, level, logger):
        # Using the child logger to handle the logging 
        for line in iter(stream.readline, ''):
            logger.log(level, "\t{}".format(line.strip()))

    def initialize(self, sample_info):
        self.create_child_logger("initialization")
        for command in self.initialization:
            if self.clone_shadow: 
                command = self.decode_shadow_commands(command)
                self.logger.info(f" Run initialize command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> ")
                self.run_command(command=command, child_logger=self.child_logger)
        


    def execute(self, input_data, sample_info):
        self.update_sample_logger(sample_info)
        self.initialize(sample_info['uuid'])
        self.logger.info(f"Executing sample {sample_info['uuid']} in {self.name} on inputs: {input_data}")
        self.load_input(input_data=input_data)


        self.remove_running_file_log()
        return {"TestOutPutVar": 50.}

    def load_input(self, input_data):
        from IOs.IOs import IOfile
        for ffile in self.input:
            input_file = IOfile.create(
                ffile["name"], 
                path=ffile['path'],
                file_type=ffile["type"],
                variables=ffile["variables"],
                save=ffile['save'],
                logger=self.logger, 
                PackID=self.PackID
            )
            input_file.write(input_data)




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

