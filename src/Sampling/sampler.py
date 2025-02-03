#!/usr/bin/env python3 
from base import Base
from abc import ABCMeta, abstractmethod
from variables import Variable
import sympy as sp 
from inner_func import update_const, update_funcs
import numpy as np 
import sys, os 

class BoolConversionError(Exception):
    """Nothing but raise bool error"""
    pass
class SamplingVirtial(Base):
    __metaclass__ = ABCMeta
    def __init__(self) -> None:
        super().__init__()
        self.schema                     = None
        self.info: Dict[str, Any]       = {}
        self.vars: Tuple[Variable]      = None
        self.method                     = None
        self.max_workers                = 4
        self._selectionexp              = None

    def set_max_workers(self, nn):
        self.max_workers = nn 

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
    def set_factory(self) -> None: 
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

    @abstractmethod 
    def combine_data(self, df_full) -> None:
        pass 
        
    @abstractmethod
    def evaluate_selection(self, expression, variables) -> bool: 
        """
        Evaluates a selection condition.

        Args:
            expression (str): A sympy-compatible condition string (e.g., "X > Y + log(E)").
            variables (dict): A dictionary of variable values (e.g., {"X": 5.0, "Y": 3.0}).

        Returns:
            bool: True if the condition is satisfied, False otherwise.

        Raises:
            BoolConversionError: If the result cannot be converted to a boolean value.
        """
        custom_functions = update_funcs({})
        custom_constants = update_const({})
        try:
            symbols = {var: sp.symbols(var) for var in variables}
            locals_context = {**custom_functions, **custom_constants, **symbols}
            expr = sp.sympify(expression, locals=locals_context)
            result = expr.subs(variables)
            result = bool(result)
            return result 
        except:
            raise BoolConversionError("Result cannot be converted to a boolean value.")

    @abstractmethod
    def check_evaluation(self):
        if self._selectionexp:
            try:
                temp    = np.random.rand(self._dimensions)
                param   = self.map_point_into_distribution(temp)
                self.evaluate_selection(self._selectionexp, param)
            except BoolConversionError:
                self.logger.error("Wrong selection condition in input YAML -> \n\t{}".format(self._selectionexp))
                sys.exit(2)
            except:
                self.logger.error("Random Sampler meets error when trying scan the parameter space.")
                sys.exit(2)