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

class Parameters(Module):
    def __init__(self, name, parameter_definitions=None):
        super().__init__(name)
        # Create Parameter object list via the information 
        if parameter_definitions is not None:
            self.parameters = [Parameter(**param_def) for param_def in parameter_definitions]
        self.unique_id = None
        self.type       = "Parameter"

    def generate_unique_id(self):
        return str(uuid.uuid4())

    def analyze_ios(self):
        for pars in self.parameters:
            self.outputs[pars.name] = None

    def generate_parameters(self):
        # Generate a unique uuid series if no uuid assigned 
        if self.unique_id is None:
            self.unique_id = self.generate_unique_id()        
        
        output_values = {}
        for param in self.parameters:
            output_values[param.name] = param.generate_value()
        
        # 将生成的参数值设置为本模块的输出
        self.outputs = output_values  # 此处的outputs是字典形式，根据需要可以调整格式

    def execute(self):
        # 调用generate_parameters方法来生成参数，并设置outputs
        self.generate_parameters()
        # 假设执行方法不需要返回值，因为输出已经被设置在self.outputs上

