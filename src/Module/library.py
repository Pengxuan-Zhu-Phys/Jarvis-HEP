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
from loguru import logger
import asyncio
class LibraryModule(Module):
    def __init__(self, name, required_modules, installed, installation):
        super().__init__(name)
        self.required_modules = required_modules
        self.installed = installed
        self.installation = installation
        self.logger = None
        self._skip_library = False
        self.loggerID = None

    def set_logger(self):
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
                # def get_sample_logger(self, sample_id, module):
                    #  sample_log_file = os.path.join(self.log_directory, f"{sample_id}_RUNNING.log")
                    #  custom_format = "{time} {level} - {name} | {message}"
                    #  sample_logger = logger.bind(sample_id=sample_id, module=module)

                    #  sample_logger.add(sample_log_file, format=custom_format, level="INFO", rotation=None, retention=None)
                    #  return sample_logger
        def custom_format(record):
            module = record["extra"].get("module", "No module")
            if "raw" in record["extra"]:
                return "{message}"
            else:
                return f"\n <cyan>{module}</cyan> \n\t- <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{record['message']}</level> "

        def filte_func(record):
            # print(record['extra'].keys(), record['extra']['module'], f"Library.{self.name}", record['extra']['module'] == f"Library.{self.name}")
            return record['extra']['module'] == f"Library.{self.name}"

        slogger = logger.bind(module=f"Library.{self.name}", to_console=True, Jarvis=True)
        self.loggerID = slogger.add(self.path['log_file_path'], format=custom_format, level="DEBUG", rotation=None, retention=None, filter=filte_func)
        # slogger.add(self.path['log_file_path'], format=custom_format, level="DEBUG", rotation=None, retention=None)
        self.logger = slogger
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
                self.logger.warning(f"Configuration for {self.name} has not changed. ")
                sleep(0.01)
                user_input = input(f"\nReinstall? (y/n): ")
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
        self.logger.remove(self.loggerID)

    def run_install_commands(self):
        for command in self.installation['commands']:
            self.logger.info(f" Run command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> \n")
            asyncio.run(self.run_command(command))

    async def run_command(self, command):
        # Create subprocess
        process = await asyncio.create_subprocess_shell(
            command['cmd'],
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=command['cwd']
        )

        # Gather both stdout and stderr concurrently
        stdout, stderr = await asyncio.gather(
            self.log_stream_info(process.stdout),
            self.log_stream_error(process.stderr)
        )

        # Wait for the process to finish
        await process.wait()


    async def log_stream_info(self, stream):
        async for line in stream:
            self.logger.bind(raw=True).info(f"\t{line.decode()}")

    async def log_stream_error(self, stream):
        async for line in stream:
            self.logger.bind(raw=True).error(f"\t{line.decode()}")

