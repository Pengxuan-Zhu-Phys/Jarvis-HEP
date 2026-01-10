#!/usr/bin/env python3
from __future__ import annotations

import yaml 
import platform
import os, sys 
import re 
import subprocess
import pkg_resources
from base import Base
import json
from jsonschema import validate
from jsonschema.exceptions import ValidationError
from Module.module import Module

class ConfigLoader(Base):
    def __init__(self) -> None:
        super().__init__()
        self.filepath   = None
        self.config     = None
        self.schema     = None
        self.logger     = None 
        self.legal      = False 
        self.summary    = {}
        self.path       = None 

    def get_cmd_output(self, cmd):
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.stderr:
            self.logger.warning("An error happend when excute: \n -> {} \n\t{}".format(cmd, result.stderr))    
        return result.stdout.strip()

    def set_schema(self, schema) -> None:
        self.schema = schema

    def load_config(self, filepath):
        try:
            self.logger.info(f"Start configuring the input file -> {filepath}")
            with open(filepath, 'r') as file:
                self.config = yaml.safe_load(file)
                self.filepath = filepath
        except Exception as e: 
            self.logger.error(f"Error: loading config file -> {e}")
            sys.exit(2)

    def check_dependency_installed(self) -> None:
        dependencies = self.config["EnvReqs"]
        # Update the environment requirement by default setting 
        if "Check_default_dependences" in dependencies:
            self.update_dependences()
        # Check the system version
        if "OS" in dependencies:
            self.check_OS_requirement(dependencies["OS"])
        # Check the CERN-ROOT setting
        if "CERN_ROOT" in dependencies:
            self.check_ROOT()
        # Check the Python Dependences 
        if "Python" in dependencies:
            self.check_PYTHON_env()
        
    
    def get_sampling_method(self) -> str:
        try: 
            return self.config['Sampling']['Method']
        except:
            self.logger.error(f"Iegal configparser file, No Sampling method founded -> {self.filepath}")
            sys.exit(2)

    def validate_config(self) -> None: 
        validator = ConfigValidator()
        validator.logger = self.logger
        validator.set_config(self.config)
        validator.set_schema(self.schema)
        validator.validate_yaml()

    def check_PYTHON_env(self) -> None: 
        py_env = self.config['EnvReqs']['Python']
        def compare_versions(version1, version2):
            v1 = tuple(map(int, version1.split('.')))
            v2 = tuple(map(int, version2.split('.')))
            return v1 >= v2
        
        def check_package_requirement(name, required, min_version):
            try:
                dist = pkg_resources.get_distribution(name)
                if pkg_resources.parse_version(dist.version) >= pkg_resources.parse_version(min_version):
                    self.logger.info(f"{name} is found, version: {dist.version} meets the requirement")
                    self.summary[name] = dist.version
                else:
                    self.logger.info(f"{name} version: {dist.version} is installed but does not meet the requirement")
                    subprocess.run(f"python -m pip install --upgrade {name} --user", shell=True)
                    dist = pkg_resources.get_distribution(name)
                    self.logger.warning(f"Jarvis-HEP is trying to upgrade {name}, please check the current version")
                    self.summary[name] = dist.version
            except pkg_resources.DistributionNotFound:
                if required:
                    self.logger.warning(f"{name} is required but not installed, Jarvis-HEP is trying to install it via pip")
                    subprocess.run(f"python3 -m pip install {name}", shell=True)
                    self.summary[name] = dist.version
                else:
                    self.logger.warning(f"{name} is optional and not installed")
                    self.summary[name] = None
       
        vf = sys.version_info
        py_version = f"{vf.major}.{vf.minor}.{vf.micro}"
        version_matched = re.match(r">=(\d+(?:\.\d+)?)", py_env['version'])
        if version_matched:
            min_version = version_matched.group(1)
            if compare_versions(py_version, min_version):
                self.summary["Python version"] = py_version
                self.logger.info(f"Python-{py_version} found.")
            else:
                self.logger.error(f"Python-{py_version} version does not meet the requirement >= {min_version}.")
                sys.exit(2)
        
        if "Dependencies" in py_env:
            for lib in py_env["Dependencies"]:
                try:
                    min_version = re.match(r">=(\d+(?:\.\d+)?)", lib['version']).group(1)
                    check_package_requirement(lib['name'], lib['required'], min_version)
                except:
                    self.logger.warning(f"Python library {lib['name']} checking un-successful!")

    def check_ROOT(self) -> None:
        def compare_versions(version1, version2):
            version1 = version1.replace("/", '.')
            v1 = tuple(map(int, version1.split('.')))
            v2 = tuple(map(int, version2.split('.')))
            return v1 >= v2
        
        root_req = self.config['EnvReqs']['CERN_ROOT']
        if root_req['required']:
            try:
                result = subprocess.run("root-config --version", shell=True, capture_output=True, text=True)
                self.summary['ROOT'] = True
            except subprocess.CalledProcessError:
                self.logger.error("CERN ROOT is not installed or root-config is not in the PATH")
                self.summary['ROOT'] = False
            except FileNotFoundError:
                self.logger.error("CERN ROOT is not installed or root-config is not in the PATH")
                self.summary['ROOT'] = False

            if self.summary['ROOT']:
                if "get_path_command" in root_req:
                    root_prefix = self.get_cmd_output(root_req['get_path_command'])
                elif "path" in root_req:
                    root_prefix = self.decode_path(root_req['path'])
                self.summary['ROOT path'] = root_prefix

                if "version" in root_req:
                    root_version = self.get_cmd_output("root-config --version")
                    version_matched = re.match(r">=(\d+(?:\.\d+)?)", root_req['version'])
                    if version_matched:
                        min_version = version_matched.group(1)
                        if compare_versions(root_version, min_version):
                            self.summary["ROOT version"] = root_version
                            self.logger.info(f"CERN ROOT-{root_version} found.")
                        else:
                            self.logger.error(f"ROOT Version-{root_version} version does not meet the requirement >= {min_version}.")
                            sys.exit(2)
                
                if "Dependencies" in root_req:
                    for dep in root_req['Dependencies']:
                        if dep['required']:
                            result = subprocess.run(dep['check_command'], shell=True, capture_output=True, text=True)
                            if not result.stdout.strip() == dep['expected_output']:
                                self.logger.error("ROOT-{} not founded, please check the CERN ROOT installation".format(dep['name']))
                                self.summary['ROOT-{}'.format(dep['name'])] = False
                                self.summary['ROOT'] = False
                            else:
                                self.summary['ROOT-{}'.format(dep['name'])] = True 

            if not self.summary['ROOT']:
                self.logger.error("ROOT installation is not meeting the requirements!")
                sys.exit(2)
            else:
                self.logger.info("CERN ROOT meets the requirements!")

    def check_OS_requirement(self, os_requirement) -> None: 
        # def compare_versions(version1, version2):
        #     v1 = tuple(map(int, version1.split('.')))
        #     v2 = tuple(map(int, version2.split('.')))
        #     return v1 >= v2

        def compare_versions(version1, version2):
            # Extract numeric parts of the version strings
            numeric_version1 = re.findall(r'\d+', version1)
            numeric_version2 = re.findall(r'\d+', version2)

            # Convert numeric strings to integers
            v1 = tuple(map(int, numeric_version1))
            v2 = tuple(map(int, numeric_version2))

            # Compare version tuples
            return v1 >= v2



        current_os = platform.system()
        current_version = platform.release()
        self.summary['OS'] = "{}-{}".format(current_os, current_version)

        os_matched = False
        for os_req in os_requirement:
            if os_req['name'].lower() == current_os.lower(): 
                version_matched = re.match(r">=(\d+(?:\.\d+)?)", os_req['version'])
                if version_matched:
                    min_version = version_matched.group(1)
                    if compare_versions(current_version, min_version):
                        os_matched = True
                        self.logger.info(f"{current_os} version {current_version} meets the requirement >= {min_version}.")
                    else:
                        os_matched = False
                        self.logger.critical(f"{current_os} version {current_version} does not meet the requirement >= {min_version}.")
                break
        if not os_matched:
            sys.exit(2)

    def update_dependences(self) -> None:
        if "required" in self.config['EnvReqs']['Check_default_dependences'] and "default_yaml_path" in self.config['EnvReqs']['Check_default_dependences']: 
            if self.config['EnvReqs']['Check_default_dependences']['required']:
                try:
                    self.config['EnvReqs']['Check_default_dependences']['default_yaml_path'] = self.decode_path(self.config['EnvReqs']['Check_default_dependences']['default_yaml_path'])
                    with open(self.decode_path(self.config['EnvReqs']['Check_default_dependences']['default_yaml_path']), 'r') as file:
                        default_env = yaml.safe_load(file)
                        self.config['EnvReqs'].update(default_env['EnvReqs'])
                except FileExistsError:
                    self.logger.error("Jarvis-HEP load file error: {} not found".format(self.config['EnvReqs']['default_yaml_path']))
        self.logger.info("Updating the Environment Requirements from default setting file \n\t{}".format(self.config['EnvReqs']['Check_default_dependences']['default_yaml_path']))

    def get_modules(self):
        modules = {
            "Parameter": self.config['Sampling']['Variables']
        }
        
        if self.config['Sampling'].get("Nuisance", False): 
            modules['Nuisance'] = self.config['Sampling']['Nuisance']['Variables']
        
        if hasattr(self.config, "LibDeps"):
            if self.config["LibDeps"]["Modules"]:
                self.analysis_Library()
                modules['Library'] = self.config["LibDeps"]["Modules"]
        if self.config["Calculators"]["Modules"]:
            self.analysis_calculator()
            modules['Calculator'] = self.config["Calculators"]["Modules"]
        return modules
        
    def decode_lib_command_via_config(self, command, pwd, pack):
        def parse_cd_command(command, cwd):
            parts = command.split()
            if len(parts) == 2 and parts[0] == "cd":
                return parts[1]
            else:
                return cwd 
        cmd = self.resolve_placeholders_config(self.config, command)
        cmd = self.resolve_placeholders_config(pack, cmd)
        cmd = self.resolve_placeholders_summary(cmd)
        cmdd =  {"cmd": cmd, "cwd": pwd}
        cwd = parse_cd_command(cmd, pwd)
        return cmdd, cwd  

    def analysis_Library(self):
        def analysis_path(lib):
            lib['installation']['path'] = self.decode_path(lib['installation']['path'])
            lib['installation']['source'] = self.decode_path(lib['installation']['source'])
            return lib 

        def analysis_commands(lib):
            cmds = lib['installation']['commands']
            cwd  = self.config['LibDeps']['path']
            commands = []
            for command in cmds:
                cmd, cwd = self.decode_lib_command_via_config(command, cwd, lib['installation'])
                commands.append(cmd)
            lib['installation']['commands'] = commands
            return lib 

        self.config["LibDeps"]['path'] = self.decode_path(self.config["LibDeps"]['path'])
        libs_config = self.config["LibDeps"]["Modules"]
        for lib in libs_config:
            lib.update(analysis_path(lib))
            lib.update(analysis_commands(lib))

    def decode_calc_command_via_config(self, command, pwd, pack):
        def parse_cd_command(command, cwd):
            parts = command.split()
            if len(parts) == 2 and parts[0] == "cd":
                return parts[1]
            else:
                return cwd 
        # print("0 ->", command)
        cmd = self.resolve_placeholders_config(self.config, command)
        # print("1 ->", cmd)
        cmd = self.resolve_placeholders_config(pack, cmd)
        # print("2 ->", cmd)
        cmd = self.resolve_placeholders_summary(cmd)
        # print("3 ->", cmd)
        cmdd =  {"cmd": cmd, "cwd": pwd}
        cwd = parse_cd_command(cmd, pwd)
        return cmdd, cwd  

    def analysis_calculator(self):
        def analysis_path(calc):
            if 'path' in calc:
                calc['path'] = self.decode_path(calc['path'])
            if "source" in calc:
                calc['source'] = self.decode_path(calc['source'])
            return calc 

        def analysis_install_commands(calc):
            cmds = calc['installation']
            cwd  = self.config['Calculators']['path']
            commands = []
            for command in cmds:
                cmd, cwd = self.decode_calc_command_via_config(command, cwd, calc)
                commands.append(cmd)
            calc['installation'] = commands
            return calc 
        
        def analysis_initialization_commands_single(calc):
            cmds = calc['initialization']
            cwd  = calc['path']
            commands = []
            for command in cmds:
                cmd, cwd = self.decode_calc_command_via_config(command, cwd, calc)
                commands.append(cmd)
            calc['initialization'] = commands
            return calc
        
        def analysis_execution_single(calc):
            calc_setting = calc['execution']

            calc_setting['path'] = self.decode_path(calc_setting['path'])

            commands = []
            for command in calc_setting['commands']:
                cmd, cwd = self.decode_calc_command_via_config(command, calc_setting['path'], calc_setting)
                commands.append(cmd)
            calc_setting['commands'] = commands

            for ipf in calc_setting['input']:
                ipf['path'] = self.decode_path(ipf['path'])

            for opf in calc_setting['output']:
                opf['path'] = self.decode_path(opf['path'])

            return calc

        self.config['Calculators']['path'] = self.decode_path(self.config['Calculators']['path'])
        calc_config = self.config["Calculators"]['Modules']
        for calc in calc_config:
            calc.update(analysis_path(calc))
            calc.update(analysis_install_commands(calc))
            if "modes" in calc:
                self.logger.info(f"Loading configuration of Module {calc['name']} in multiple modes")
            else:
                self.logger.info(f"Loading configuration of Module {calc['name']} in single mode")
                calc['modes'] = False
                calc.update(analysis_initialization_commands_single(calc))
                calc.update(analysis_execution_single(calc))

    def decode_path(self, path) -> None:
        return super().decode_path(path)

    def resolve_placeholders_config(self, data, text, separator=':'):
        pattern = re.compile(r'\$\{([^}]+)\}')
        def find_value(d, keys):
            if keys and keys[0] in d:
                return find_value(d[keys[0]], keys[1:]) if len(keys) > 1 else d[keys[0]]
            return None
        for match in pattern.finditer(text):
            key_path = match.group(1)
            keys = key_path.split(separator)
            value = find_value(data, keys)
            if value is not None:
                text = text.replace(match.group(0), str(value))
        return text
    
    def resolve_placeholders_summary(self, text):
        pattern = re.compile(r'\@\{([^}]+)\}')
        def find_value(d, keys):
            if keys and keys[0] in d:
                return find_value(d[keys[0]], keys[1:]) if len(keys) > 1 else d[keys[0]]
            return None

        for match in pattern.finditer(text):
            key_path = match.group(1)
            keys = key_path.split(":")
            value = find_value(self.summary, keys)
            if value is not None:
                text = text.replace(match.group(0), str(value).strip())
        return text





class ConfigValidator(Base):
    def __init__(self) -> None:
        super().__init__()
        self.logger = None
        self.config = None
        self.schema = None
        self.passcheck = False

    def set_config(self, config):
        self.config = config

    def set_schema(self, schema):
        with open(schema, 'r') as file:
            self.schema = json.load(file)

    def validate_yaml(self):
        if self.config is None or self.schema is None:
            self.logger.error("Error: Config or schema not set.")
            sys.exit(2)
        try:
            # print(self.schema['schemaBlock'].keys())
            for kk, vv in self.schemablock.items():
                self.schema['schemaBlock'][kk]['$ref'] = self.schemablock[kk]
            
            validate(instance=self.config, schema=self.schema)
            self.logger.warning("Validation successful. The input YAML file meets the schema requirement.")
            self.passcheck = True
        except ValidationError as e:
            self.logger.error(f"Validation error: {e.message},\n the problematic module data: {e.instance}")


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
