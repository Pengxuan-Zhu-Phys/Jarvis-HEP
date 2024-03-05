#!/usr/bin/env python3

import configparser
from copy import deepcopy
import json
import os
from subprocess import Popen, run
import sys
import time 
from program import Pack
import argparse
from config import ConfigLoader
from logger import AppLogger
import logging
from base import Base
from sampler import Distributor

class Core(Base):
    def __init__(self) -> None:
        # print("Testing the first logging line")
        super().__init__()
        self.argparser: Any             = None
        self.info: Dict[str, Any]       = {}
        self.yaml: ConfigLoader         = ConfigLoader()
        self.sampler: Any               = None
        self.factory: Any               = None 
        self.likelihood: Any            = None 
        self.logger: AppLogger          = None

        
    def init_argparser(self) -> None:
        self.argparser = argparse.ArgumentParser(description="Jarvis Program Help Center")
        self.info['args'] = load_args_config(self.path['args_info'])
        
        for pos_arg in self.info['args'].get('positionals', []):
            self.argparser.add_argument(pos_arg['name'], help=pos_arg['help'])
        
        for opt in self.info['args'].get('options', []):
            kwargs = {
                'help': opt['help'],
                'action': opt.get('action', 'store'),
                'dest': opt.get('dest')
            }
            if 'default' in opt:
                kwargs['default'] = opt['default']
            if "type" in opt:
                if opt['type'] == 'int':
                    kwargs['type'] = int
                elif opt['type'] == 'float':
                    kwargs['type'] = float
                else:
                    kwargs['type'] = str
            if 'short' in opt and 'long' in opt:
                self.argparser.add_argument(
                    opt['short'],
                    opt['long'],
                    **kwargs
                )
            elif 'long' in opt:
                self.argparser.add_argument(
                    opt['long'], **kwargs
                )
        self.check_init_args()

    def init_logger(self) -> None:
        self.info["project_name"] = os.path.splitext(os.path.basename(self.args.file))[0]
        self.logger = AppLogger(
            config_path=self.path['logger_config_path'],
            logger_name="Jarvis-HEP",
            log_file_name=f"{self.info['project_name']}.log"
        )
        self.logger.info['debug_mode'] = self.args.debug
        self.logger.configure_logging()
        self.logger.print_logo()
        self.logger.logger.info("Jarvis-HEP logging system initialized successful!")
        if self.args.debug:
            self.logger.logger.info("Jarvis-HEP in debug mode currently!")
        
        # child_logger = self.logger.create_dynamic_logger("Test_Child_Logger", log_file="child.log")
        # self.logger.delete_child_logger("Test_Child_Logger", child_logger)

    def init_configparser(self) -> None: 
        self.yaml.logger = self.logger.create_dynamic_logger("ConfigParser")
        from copy import deepcopy
        self.yaml.path = deepcopy(self.path)
        self.yaml.load_config(os.path.abspath(self.args.file))
        self.yaml.check_dependency_installed()
        self.sampler = Distributor.set_method(self.yaml.get_sampling_method()) 
        self.yaml.set_schema(self.sampler.schema)
        self.yaml.validate_config()


    def check_init_args(self) -> None:
        try:
            self.args = self.argparser.parse_args()
        except argparse.ArgumentError as e:
            print(str(e))
            self.argparser.print_help()
            sys.exit(2)

    def initialization(self) -> None:
        self.init_argparser()
        self.init_logger()
        self.init_configparser()





def load_args_config(json_file):
    with open(json_file, 'r') as file:
        config = json.load(file)
    return config
