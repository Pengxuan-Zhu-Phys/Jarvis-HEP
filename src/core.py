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
import argparse

class Core():
    def __init__(self) -> None:
        self.argparser  = None
        self.info       = {}
        self.path = {
            "jpath":    os.path.abspath(os.path.dirname(__file__))            
        }
        
        self.init_argparser()

    def init_argparser(self) -> None:
        self.argparser = argparse.ArgumentParser(description="Jarvis Program Help Center")
        self.path['args_info'] = os.path.join(self.path['jpath'], "card/argparser.json")
        self.info['args'] = load_config(self.path['args_info'])
        
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

        # print(self.path)
        self.check_init_args()

    def check_init_args(self) -> None:
        try:
            self.args = self.argparser.parse_args()
            print(self.args)
        except argparse.ArgumentError as e: 
            print(str(e))
            self.argparser.print_help()
            sys.exit(2)





def load_config(json_file):
    with open(json_file, 'r') as file:
        config = json.load(file)
    return config