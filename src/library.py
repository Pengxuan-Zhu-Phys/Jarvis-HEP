#!/usr/bin/env python3 

from copy import deepcopy
import json
import logging
import os, sys
from base import Base
from pprint import pprint
jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
import yaml 
from Module.library import LibraryModule


class Library(Base):
    def __init__(self) -> None:
        super().__init__()
        self.modules = {}
        self._skip_library = False

    def set_logger(self, logger) -> None:
        self.logger = logger 
        self.logger.warning("Ready for preparing the Supporting Libraries ...")
        
    def set_config(self, config):
        self.config = config
        if self.config:
            self.config['SupportingLibrary']["path"] = self.decode_path(self.config['SupportingLibrary']['path'])
            self.path['library'] = self.config['SupportingLibrary']["path"]
            for module in self.config['SupportingLibrary']['Modules']:
                mod = LibraryModule(module['name'], module['required_modules'], installed=False, installation=module['installation'])
                mod._skip_library = self._skip_library
                mod.set_library_card(self.path['library'])
                mod.set_config(module)
                mod.set_logger()
                self.modules[module['name']] = mod

    def display_installation_summary(self):
        self.logger.info("Preparing to install libraries. Installation summary:")
        installed_modules = []
        not_installed_modules = []
        
        for module in self.modules.values():
            module.check_installed()  
            if module.installed:
                installed_modules.append(module.name)
            else:
                not_installed_modules.append(module.name)
        
        from prettytable import PrettyTable

        table = PrettyTable()
        table.field_names = ["Module Name", "Status"]

        for name in installed_modules:
            table.add_row([name, "Already Installed"])
        for name in not_installed_modules:
            table.add_row([name, "Will be Installed"])
        if not installed_modules and not not_installed_modules:
            table.add_row(["No modules", "N/A"])

        self.logger.warning(table)


