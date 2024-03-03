#!/usr/bin/env python3

import yaml 
import importlib
import math 
import os, sys 

class ConfigLoader():
    def __init__(self) -> None:
        self.filepath   = None
        self.config     = None
        self.schema     = None

    def set_schema(self, schema) -> None:
        self.schema = schema

    def load_yaml(self, filepath):
        try:
            with open(filepath, 'r') as file:
                self.config = yaml.safe_load(file)
                self.filepath = filepath
        except Exception as e: 
            print(f"Error: loading config file -> {e}")
            sys.exit(2)

    def validate_config(self) -> None: 
        validator = ConfigValidator(self.config, self.schema)

class ConfigValidator():
    def __init__(self) -> None:
        self.config = None
        self.schema = None
        self.passcheck = False

    def set_config(self, config):
        self.config = config

    def set_schema(self, schema):
        self.schema = schema

    def validate(self):
        if self.config is None or self.schema is None:
            print("Error: Config or schema not set.")
            sys.exit(2)
        
        missing_sections = [section for section in self.schema if section not in self.config]
        if missing_sections:
            print(f"配置文件缺少必需的部分: {missing_sections}")
            return False
        
        for section, params in self.schema.items():
            missing_params = [param for param in params if param not in self.config.get(section, {})]
            if missing_params:
                print(f"在 '{section}' 部分缺少参数: {missing_params}")
                return False
        
        self.passcheck = True

class ConfigTemplateGenerator():
    def __init__(self) -> None:
        self.schema = None
        self.output_path = None

    def set_schema(self, schema):
        self.schema = schema

    def set_output_path(self, output_path):
        self.output_path = output_path

    def generate_template(self):
        if self.schema is None or self.output_path is None:
            print("Error: Schema or output path not set.")
            sys.exit(2)
        
        template = {}
        for section, params in self.schema.items():
            template[section] = {param: '' for param in params}
        
        try:
            with open(self.output_path, 'w') as file:
                yaml.dump(template, file, default_flow_style=False)
            print(f"Template file is generated: {self.output_path}")
        except Exception as e:
            print(f"An error happened in generating template file: {e}")
            sys.exit(2)
