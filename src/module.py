#!/usr/bin/env python3 
from IOs import Parameter
import uuid
from pprint import pprint
import yaml 
import logging
import subprocess
import threading
import os
from base import Base
from time import sleep

class Module(Base): 
    def __init__(self, name, inputs=None, outputs=None):
        super().__init__()
        self.logger = None
        self.name = name
        self.input = inputs if inputs else []
        self.output = outputs if outputs else []
        self.inputs = {}
        self.outputs = {}
        self.required_modules = []

    def execute(self):
        raise NotImplementedError("This method should be implemented by subclasses.")


class Parameters(Module):
    def __init__(self, name, parameter_definitions=None):
        super().__init__(name)
        # Create Parameter object list via the information 
        if parameter_definitions is not None:
            self.parameters = [Parameter(**param_def) for param_def in parameter_definitions]
        self.unique_id = None
        self.type       = "Parameter"

    def generate_unique_id(self):
        return str(uuid.uuid4())

    def analyze_ios(self):
        for pars in self.parameters:
            self.outputs[pars.name] = None

    def generate_parameters(self):
        # Generate a unique uuid series if no uuid assigned 
        if self.unique_id is None:
            self.unique_id = self.generate_unique_id()        
        
        output_values = {}
        for param in self.parameters:
            output_values[param.name] = param.generate_value()
        
        # 将生成的参数值设置为本模块的输出
        self.outputs = output_values  # 此处的outputs是字典形式，根据需要可以调整格式

    def execute(self):
        # 调用generate_parameters方法来生成参数，并设置outputs
        self.generate_parameters()
        # 假设执行方法不需要返回值，因为输出已经被设置在self.outputs上

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

    def assign_ID(self, PackID):
        self.PackID = PackID 
        print(self.config.keys())

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

    def install(self):
        self.basepath = self.decode_shadow_path(self.basepath)
        if not os.path.exists(self.basepath):
            os.makedirs(self.basepath)

        self.logger.setLevel(logging.DEBUG)  # 或更低的级别

        install_file_log = os.path.join(self.basepath, f"Installation_{self.name}-{self.PackID}.log")
        formatter = logging.Formatter(" %(name)s - %(asctime)s - %(levelname)s >>> \n%(message)s")
        Afile_handler = logging.FileHandler(install_file_log, mode='a')
        Afile_handler.setLevel(0)
        Afile_handler.setFormatter(formatter)
        self.logger.addHandler(Afile_handler)
        self.logger.info(f"Installing {self.name}-{self.PackID} ....")


        child_installation_logger_name = f"{self.logger.name}.command_logger"
        self.child_installation_logger = logging.getLogger(child_installation_logger_name)
        self.child_installation_logger.propagate = False
        
        # setting a simple formmater for each command 
        simple_formatter = logging.Formatter('%(message)s')
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(simple_formatter)
        stream_handler.setLevel(logging.DEBUG)
        # self.child_installation_logger.addHandler(stream_handler)

        file_handler = logging.FileHandler(install_file_log, mode="a")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(simple_formatter)
        self.child_installation_logger.addHandler(file_handler)


        for cmd in self.installation:
            if self.clone_shadow:
                command = self.decode_shadow_commands(cmd)
                self.logger.info(f" Run command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> \n")
                self.run_command(command=command, child_logger=self.child_installation_logger)
            else: 
                self.run_command(command=cmd, child_logger=self.child_installation_logger)

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


    def initialize(self):
        # 模拟初始化逻辑
        print(f"Initializing {self.name} with: {self.initialization}")

    def execute(self, input_data):
        # 执行模块计算逻辑，处理输入和产生输出
        print(f"Executing {self.name} on inputs: {input_data}")
        # 此处应根据模块具体逻辑处理输入数据和执行命令

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


class LibraryModule(Module):
    def __init__(self, name, required_modules, installed, installation):
        super().__init__(name)
        self.required_modules = required_modules
        self.installed = installed
        self.installation = installation
        self.logger = None
        self._skip_library = False

    def set_logger(self, logger):
        self.logger = logger
        self.logger.warning("Initializating Library -> {}".format(self.name))

    def check_installed(self):
        # Checking the installation configuration file
        if os.path.exists(self.path['config_file_path']):
            with open(self.path['config_file_path'], 'r') as file:
                existing_config = yaml.safe_load(file)
            # Comparing the settings
            if existing_config['installed'] != True:
                self.logger.info(f"{self.name} has found, but not be correctly installed.")
                self.installed = False  # Setting the installation to false due to config change
            else: 
                self.logger.info(f"{self.name} has found and correctly installed.")
                existing_config['installed'] = False
                if existing_config == self.config:
                    self.logger.info(f"{self.name} is already installed with matching configuration.")
                    self.installed = True
                else: 
                    self.logger.info(f"Configuration for {self.name} has changed.")
        else:
            self.logger.info(f"No existing installation found for {self.name}.")
            self.installed = False


    def set_config(self, config):
        self.config = config 

    def set_library_card(self, path):
        self.path['log_file_path'] = os.path.join(path, f"Library_{self.name}_installation.log")
        self.path['config_file_path'] = os.path.join(path, f'{self.name}_config.yaml')

    def install(self):
        # self.check_installed()
        if not self._skip_library:
            if self.installed:
                user_input = input(f"Configuration for {self.name} has not changed. Reinstall? (y/n): ")
                if user_input.lower() != 'y':
                    self.logger.warning(f"Skipping installation of {self.name}.")
                else:
                    self.run_install_commands()
            if not self.installed:
                self.run_install_commands()

            self.config['installed'] = True
            self.installed = True

            with open(self.path['config_file_path'], 'w') as file:
                yaml.dump(self.config, file)
                self.logger.info(f"Configuration for {self.name} has been written to {self.path['config_file_path']}")
        else: 
            self.logger.warning(f"Skipping the installation of library -> {self.name}")

    def run_install_commands(self):
        child_logger_name = f"{self.logger.name}.command_logger"
        self.child_logger = logging.getLogger(child_logger_name)
        self.child_logger.propagate = False
        
        # setting a simple formmater for each command 
        simple_formatter = logging.Formatter('%(message)s')
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(simple_formatter)
        stream_handler.setLevel(logging.INFO)
        self.child_logger.addHandler(stream_handler)

        file_handler = logging.FileHandler(self.path['log_file_path'], mode="a")
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(simple_formatter)
        self.child_logger.addHandler(file_handler)

        # self.child_logger.setLevel(self.logger.level)

        for command in self.installation['commands']:
            self.logger.info(f" Run command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> \n")
            self.run_command(command)

    def run_command(self, command):
        with subprocess.Popen(
                                command['cmd'], 
                                shell=True, 
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, 
                                text=True, 
                                cwd=command['cwd']) as process:
            stdout_thread = threading.Thread(target=self.log_stream, args=(process.stdout, logging.INFO, self.child_logger))
            stdout_thread.start()
            stderr_thread = threading.Thread(target=self.log_stream, args=(process.stderr, logging.ERROR, self.child_logger)) 
            stderr_thread.start()

            stdout_thread.join()
            stderr_thread.join()

            process.wait()
            # if process.returncode != 0:
                # raise subprocess.CalledProcessError(process.returncode, command)

            # child_logger.removeHandler(handler)
        # sleep(4)

    def log_stream(self, stream, level, logger):
        # Using the child logger to handle the logging 
        for line in iter(stream.readline, ''):
            logger.log(level, "\t{}".format(line.strip()))
        