#!/usr/bin/env python3

import configparser
from copy import deepcopy
import json
import os
from subprocess import Popen, run
import sys
import logging 
import logging.config
import time 
import yaml
from base import Base

pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(pwd, "Sampling"))) 


from loguru import logger

class AppLogger(Base):
    def __init__(self, logger_name='Jarvis', log_file_name="Jarvis.log", debug_mode=False):
        super().__init__()
        self.logger = logger.bind(module=logger_name)
        self.debug_mode = debug_mode

        def global_log_filter(record):
            # 只有包含 'global_log' 标记的日志消息才被写入
            return record["extra"].get("Jarvis", False)

        logger.add(f"{log_file_name}_{time}", 
                   rotation="5 MB", 
                   encoding='utf-8', 
                   format="\n {name} \n\t- {time} - [{level}] >>> \n{message}", 
                   enqueue=True, 
                   level="DEBUG" if self.debug_mode else "WARNING",
                   filter=global_log_filter
                   )


    def create_dynamic_logger(self, logger_name, level="WARNING", log_file=None):
        # Generate a dynamic logger with a specific log file if provided
        level = "DEBUG" if self.debug_mode else level
        log_file = log_file or self.log_file_name
        format_str = "\n{name} \n\t- {time} - [{level}] >>> \n{message}"

        # Add a new sink if a specific log file is requested for the dynamic logger
        if log_file:
            logger.add(log_file, rotation="1 day", level=level, format=format_str, enqueue=True)

        dlogger = logger.bind(name=f"{self.logger_name}.{logger_name}")
        return 

    def delete_child_logger(self, logger_name):
        # In loguru, you cannot 'delete' loggers as you do in Python's logging,
        # because all logs are handled by the same logger object.
        # However, you can remove specific handler if needed, but usually,
        # you control this via configurations not by adding/removing loggers.
        pass

    def print_logo(self):
        with open(self.path['logo'], 'r') as f:
            self.logger.warning(f"\n{f.read()}")

# Usage
app_logger = AppLogger()
dynamic_logger = app_logger.create_dynamic_logger("Module1", level="DEBUG")
dynamic_logger.info("This is a test log from Module1.")

# class AppLogger(Base):
#     def __init__(self, config_path=os.path.join(pwd, 'card/jarvis_logging_config.yaml'), logger_name='Jarvis', log_file_name="Jarvis.log"):
#         super().__init__()
#         self.info = {
#             "parent_logger_name":   logger_name,
#             "parent_logger_handler": None,
#             "parent_log_file_name": log_file_name,
#             "parent_config_path":   config_path,
#             "debug_mode":   False
#         }
#         # self.configure_logging(config_path, log_file_name)

#     def configure_logging(self):
#         with open(self.info['parent_config_path'], 'r') as f:
#             config = yaml.safe_load(f.read())
#             config['handlers']['file_parent']['class'] = 'logging.handlers.TimedRotatingFileHandler'
#             config['handlers']['file_parent']['filename'] = self.info['parent_log_file_name']
#             config['handlers']['file_parent']['when'] = 'D'  
#             config['handlers']['file_parent']['interval'] = 1  
#             config['handlers']['file_parent']['backupCount'] = 0  
#             config['handlers']['file_parent']['encoding'] = 'utf-8'
#             if self.info['debug_mode']:
#                 config['handlers']['console']['level'] = logging.DEBUG
#             logging.config.dictConfig(config)
#         self.logger = logging.getLogger(self.info['parent_logger_name'])

#     def get_parent_file_handlers(self):
#         for handler in self.logger.handlers:
#             if isinstance(handler, (logging.FileHandler, logging.handlers.RotatingFileHandler, logging.handlers.TimedRotatingFileHandler)):
#                 return handler

#     def get_logger(self):
#         return self.logger
    
#     def print_logo(self):
#         with open(self.path['logo'], 'r') as f1:
#             self.logger.warning(f"\n{f1.read()}")

#     def create_dynamic_logger(self, logger_name, level=logging.WARNING, log_file=None):
#         """
#         Dynamically creates a logger.

#         :param logger_name: The name of the new logger.
#         :param level: The logging level for the new logger.
#         :param propagate: Whether to allow the log messages to propagate to the parent logger.
#         """
#         if self.info['debug_mode']:
#             level = logging.DEBUG
#         full_logger_name = f"{self.info['parent_logger_name']}.{logger_name}"
#         logger = logging.getLogger(full_logger_name)
#         logger.setLevel(level)
#         logger.propagate = False

#         if not logger.handlers:
#             console_handler = logging.StreamHandler()
#             console_handler.setLevel(level)
#             formatter = logging.Formatter("\n %(name)s \n\t- %(asctime)s - [%(levelname)s] >>> \n%(message)s")
#             console_handler.setFormatter(formatter)
#             logger.addHandler(console_handler)
#         if log_file is not None: 
#             file_handler = logging.FileHandler(log_file)
#             file_handler.setLevel(logging.DEBUG)
#             formatter = logging.Formatter("\n %(name)s \n\t- %(asctime)s - [%(levelname)s] >>> \n%(message)s")
#             file_handler.setFormatter(formatter)
#             logger.addHandler(file_handler)

#         parent_file_handler = logging.handlers.TimedRotatingFileHandler(
#             filename=self.info['parent_log_file_name'],
#             when="D",
#             interval=1,
#             backupCount=0,
#             encoding='utf-8'
#         )
#         formatter = logging.Formatter("\n %(name)s \n\t- %(asctime)s - [%(levelname)s] >>> \n%(message)s")
#         parent_file_handler.setFormatter(formatter)
#         parent_file_handler.setLevel(level)
#         logger.addHandler(parent_file_handler)

#         return logger
    
#     def delete_child_logger(self, logger_name, clogger):
#         logger_name = f"{self.info['parent_logger_name']}.{logger_name}"
#         if logger_name in logging.root.manager.loggerDict:
#             del logging.root.manager.loggerDict[logger_name]
#         del clogger
#         self.logger.warning(f"Logger is deleted-> {logger_name}")



if __name__ == "__main__":
    app_logger = AppLogger(logger_name="MainLogger")
    dynamic_logger = app_logger.create_dynamic_logger("DynamicWorker", level=logging.INFO, propagate=False)
    dynamic_logger.info("This is a warning message from DynamicWorker.")

