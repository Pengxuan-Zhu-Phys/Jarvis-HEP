#!/usr/bin/env python3

import yaml 
import platform
import importlib
import math 
import os, sys 
import re 
import subprocess
import pkg_resources

class ConfigLoader():
    def __init__(self) -> None:
        self.filepath   = None
        self.config     = None
        self.schema     = None
        self.logger     = None 
        self.legal      = False 
        self.summary    = {}
        self.path       = None 

    def decode_path(self, path) -> None: 
        """
        Resolves special markers in the provided path.

        Parameters:
        - path: The path to be resolved.
        - jarvis_root: The root directory path of Jarvis.

        Returns:
        - The resolved full path.
        """
        # Replace the Jarvis root directory marker &J
        if "&J" in path:
            path = path.replace("&J", self.path['jpath'])

        # Replace the user home directory marker ~
        if "~" in path:
            path = os.path.expanduser(path)

        return path        

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
            print(f"Error: loading config file -> {e}")
            sys.exit(2)

    def check_dependency_installed(self) -> None:
        dependencies = self.config["EnvironmentRequirements"]
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

        print(self.summary)
    

        

    def validate_config(self) -> None: 
        validator = ConfigValidator(self.config, self.schema)

    def check_PYTHON_env(self) -> None: 
        py_env = self.config['EnvironmentRequirements']['Python']
        def compare_versions(version1, version2):
            v1 = tuple(map(int, version1.split('.')))
            v2 = tuple(map(int, version2.split('.')))
            return v1 >= v2
        
        def check_package_requirement(name, required, min_version):
            try:
                dist = pkg_resources.get_distribution(name)
                print(dist, min_version)
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
                    subprocess.run(f"python -m pip install {name}", shell=True)
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
                    self.logger.warning("Python library {} checking un-successful!")


    def check_ROOT(self) -> None:
        def compare_versions(version1, version2):
            v1 = tuple(map(int, version1.split('.')))
            v2 = tuple(map(int, version2.split('.')))
            return v1 >= v2
        
        root_req = self.config['EnvironmentRequirements']['CERN_ROOT']
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
            from pprint import pprint
            pprint(self.summary)

    def check_OS_requirement(self, os_requirement) -> None: 
        def compare_versions(version1, version2):
            v1 = tuple(map(int, version1.split('.')))
            v2 = tuple(map(int, version2.split('.')))
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
        if "required" in self.config['EnvironmentRequirements']['Check_default_dependences'] and "default_yaml_path" in self.config['EnvironmentRequirements']['Check_default_dependences']: 
            if self.config['EnvironmentRequirements']['Check_default_dependences']['required']:
                try:
                    self.config['EnvironmentRequirements']['Check_default_dependences']['default_yaml_path'] = self.decode_path(self.config['EnvironmentRequirements']['Check_default_dependences']['default_yaml_path'])
                    with open(self.decode_path(self.config['EnvironmentRequirements']['Check_default_dependences']['default_yaml_path']), 'r') as file:
                        default_env = yaml.safe_load(file)
                        self.config['EnvironmentRequirements'].update(default_env['EnvironmentRequirements'])
                        from pprint import pprint
                except FileExistsError:
                    self.logger.error("Jarvis-HEP load file error: {} not found".format(self.config['EnvironmentRequirements']['default_yaml_path']))
        self.logger.info("Updating the Environment Requirements from default setting file \n\t{}".format(self.config['EnvironmentRequirements']['Check_default_dependences']['default_yaml_path']))

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
