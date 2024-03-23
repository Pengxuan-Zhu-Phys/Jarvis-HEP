#!/usr/bin/env python3 
from IOs.IOs import Parameter
import uuid
from pprint import pprint
import yaml 
import logging
import subprocess
import threading
import os
from base import Base
from time import sleep
from Module.module import Module

class LibraryModule(Module):
    def __init__(self, name, required_modules, installed, installation):
        super().__init__(name)
        self.required_modules = required_modules
        self.installed = installed
        self.installation = installation
        self.logger = None
        self._skip_library = False

    def set_logger(self, logger):
        """
            Assigns a logger instance to the module and logs an initialization message.

            This method allows an externally created logger instance to be associated with the current module,
            enabling the logging of messages within this module to be handled by the shared logger. This approach
            facilitates sharing the same logging configuration across different parts of the application and centralizes
            the management of how and where log messages are output.

            Upon assigning the logger, this method immediately uses the provided logger to record a warning level log message,
            indicating that the corresponding library module is being initialized. This helps in tracking the execution flow
            and status of the program.

            Parameters:
            - logger: A configured logging.Logger object that will be used for logging within this module.

            Returns:
            None
        """
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
        