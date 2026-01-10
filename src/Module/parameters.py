#!/usr/bin/env python3 
from IOs.IOs import Parameter
import uuid
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
        self.nuisance   = None

    def generate_unique_id(self):
        return str(uuid.uuid4())

    def add_nuisance(self, nuisances): 
        for nvar in nuisances: 
            nvar['ptype'] = "Nuisa"
            self.parameters.append(Parameter(**nvar))

    def analyze_ios(self):
        for pars in self.parameters:
            self.outputs[pars.name] = pars.type

    def generate_parameters(self):
        # Generate a unique uuid series if no uuid assigned 
        if self.unique_id is None:
            self.unique_id = self.generate_unique_id()        
        
        output_values = {}
        for param in self.parameters:
            output_values[param.name] = param.generate_value()
        
        self.outputs = output_values  

    def execute(self):
        self.generate_parameters()

