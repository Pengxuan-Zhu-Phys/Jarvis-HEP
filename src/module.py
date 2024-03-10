#!/usr/bin/env python3 
from IOs import Parameter
import uuid
from pprint import pprint

class Module: 
    def __init__(self, name, inputs=None, outputs=None):
        self.logger = None
        self.name = name
        self.input = inputs if inputs else []
        self.output = outputs if outputs else []
        self.inputs = {}
        self.outputs = {}
        self.required_modules = []

    def execute(self):
        raise NotImplementedError("This method should be implemented by subclasses.")


class Parameters(Module):
    def __init__(self, name, parameter_definitions):
        super().__init__(name)
        # Create Parameter object list via the information 
        self.parameters = [Parameter(**param_def) for param_def in parameter_definitions]
        self.unique_id = None

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

class LibraryModule(Module):
    def __init__(self, name, required_modules, installed, installation):
        super().__init__(name)
        self.required_modules = required_modules
        self.installed = installed
        self.installation = installation
        

    def execute(self):
        if not self.installed:
            # 安装逻辑
            pass

class CalculatorModule(Module):
    def __init__(self, name, config):
        super().__init__(name)
        if config["modes"]:
            self.modes = True
            self.analyze_config_multi()
        else:
            self.required_modules   = config["required_modules"]
            self.clone_shadow       = config["clone_shadow"]
            self.installation       = config["installation"]
            self.initialization     = config["initialization"]
            self.execution          = config["execution"]
            self.input              = config["execution"].get("input", [])
            self.output             = config["execution"].get("output", [])
            self.modes = False
            self.analyze_config()

    def analyze_config(self):
        # get the variables for inputs and outputs, prepare for the workflow chart 
        for ipf in self.input:
            # self.inputs[ipf['name']] = None
            for ipv in ipf['variables']:
                self.inputs[ipv['name']] = None
        
        for opf in self.output:
            # self.outputs[opf['name']] = None
            for opv in opf['variables']:
                self.outputs[opv['name']] = None


    def analyze_config_multi(self):
        pass 

    def install(self):
        # 模拟安装逻辑
        # 实际应用中可能需要执行安装命令
        print(f"Installing {self.name} using commands: {self.installation['commands']}")

    def initialize(self):
        # 模拟初始化逻辑
        print(f"Initializing {self.name} with: {self.initialization}")

    def execute(self, input_data):
        # 执行模块计算逻辑，处理输入和产生输出
        print(f"Executing {self.name} on inputs: {input_data}")
        # 此处应根据模块具体逻辑处理输入数据和执行命令

