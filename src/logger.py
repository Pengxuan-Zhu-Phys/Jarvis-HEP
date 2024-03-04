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
from program import Pack
import yaml

jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(pwd, "Sampling"))) 

class AppLogger:
    def __init__(self, config_path=os.path.join(pwd, 'card/jarvis_logging_config.yaml'), logger_name='Jarvis', log_file_name="Jarvis.log"):
        self.info = {
            "parent_logger_name":   logger_name,
            "parent_logger_handler": None,
            "parent_log_file_name": log_file_name,
            "parent_config_path":   config_path,
            "debug_mode":   False
        }
        # self.configure_logging(config_path, log_file_name)

    def configure_logging(self):
        with open(self.info['parent_config_path'], 'r') as f:
            config = yaml.safe_load(f.read())
            config['handlers']['file_parent']['class'] = 'logging.handlers.TimedRotatingFileHandler'
            config['handlers']['file_parent']['filename'] = self.info['parent_log_file_name']
            config['handlers']['file_parent']['when'] = 'D'  
            config['handlers']['file_parent']['interval'] = 1  
            config['handlers']['file_parent']['backupCount'] = 0  
            config['handlers']['file_parent']['encoding'] = 'utf-8'
            if self.info['debug_mode']:
                config['handlers']['console']['level'] = logging.DEBUG
            logging.config.dictConfig(config)
        self.logger = logging.getLogger(self.info['parent_logger_name'])

    def get_parent_file_handlers(self):
        for handler in self.logger.handlers:
            if isinstance(handler, (logging.FileHandler, logging.handlers.RotatingFileHandler, logging.handlers.TimedRotatingFileHandler)):
                return handler

    def get_logger(self):
        return self.logger
    
    def print_logo(self, logo_file=os.path.join(pwd, "card/logo")):
        with open(logo_file, 'r') as f1:
            self.logger.warning(f"\n{f1.read()}")

    def create_dynamic_logger(self, logger_name, level=logging.DEBUG, log_file=None):
        """
        Dynamically creates a logger.

        :param logger_name: The name of the new logger.
        :param level: The logging level for the new logger.
        :param propagate: Whether to allow the log messages to propagate to the parent logger.
        """
        full_logger_name = f"{self.info['parent_logger_name']}.{logger_name}"
        logger = logging.getLogger(full_logger_name)
        logger.setLevel(level)
        logger.propagate = False

        if not logger.handlers:
            console_handler = logging.StreamHandler()
            console_handler.setLevel(logging.WARNING)
            formatter = logging.Formatter(" %(name)s - %(asctime)s - %(levelname)s - \n%(message)s")
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)
        if log_file is not None: 
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter(" %(name)s - %(asctime)s - %(levelname)s - \n%(message)s")
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        parent_file_handler = logging.handlers.TimedRotatingFileHandler(
            filename=self.info['parent_log_file_name'],
            when="D",
            interval=1,
            backupCount=0,
            encoding='utf-8'
        )
        formatter = logging.Formatter(" %(name)s - %(asctime)s - %(levelname)s - \n%(message)s")
        parent_file_handler.setFormatter(formatter)
        parent_file_handler.setLevel(logging.WARNING)
        logger.addHandler(parent_file_handler)

        print(self.logger.handlers)
        print(logger.handlers)
        return logger
    
    def delete_child_logger(self, logger_name, clogger):
        logger_name = f"{self.info['parent_logger_name']}.{logger_name}"
        if logger_name in logging.root.manager.loggerDict:
            del logging.root.manager.loggerDict[logger_name]
        del clogger
        print(logging.root.manager.loggerDict)

        self.logger.warning(f"Logger is deleted-> {logger_name}")



if __name__ == "__main__":
    app_logger = AppLogger(logger_name="MainLogger")
    dynamic_logger = app_logger.create_dynamic_logger("DynamicWorker", level=logging.INFO, propagate=False)
    dynamic_logger.info("This is a warning message from DynamicWorker.")

