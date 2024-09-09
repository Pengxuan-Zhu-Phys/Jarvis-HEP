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

class Module(Base): 
    def __init__(self, name, inputs=None, outputs=None):
        super().__init__()
        self.logger = None
        self.name = name
        self.input = inputs if inputs else []
        self.output = outputs if outputs else []
        self.inputs = {}
        self.outputs = {}
        self.required_modules = []

    def execute(self):
        raise NotImplementedError("This method should be implemented by subclasses.")
