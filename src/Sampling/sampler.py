#!/usr/bin/env python3 
from base import Base
from abc import ABCMeta, abstractmethod
from variables import Variable

class SamplingVirtial(Base):
    __metaclass__ = ABCMeta
    def __init__(self) -> None:
        super().__init__()
        self.schema                     = None
        self.info: Dict[str, Any]       = {}
        self.vars: Tuple[Variable]      = None
        self.method                     = None

    @abstractmethod
    def load_schema_file(self) -> None:
        pass 

    @abstractmethod
    def set_config(self, config_info) -> None:
        pass 

    @abstractmethod
    def init_generator(self) -> None:
        pass 

    @abstractmethod
    def load_variable(self) -> None:
        variables = []
        for var_config in self.config['Sampling'].get("Variables", []):  # 使用get以防'Variables'键不存在
            var = Variable(
                name=var_config['name'],
                description=var_config['description'],
                distribution=var_config['distribution']['type'],
                parameters=var_config['distribution']['parameters']
            )
            variables.append(var)
        self.vars = tuple(variables)

    @abstractmethod
    def set_logger(self, logger) -> None:
        self.logger = logger 
        