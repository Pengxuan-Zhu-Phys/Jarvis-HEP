#!/usr/bin/env python3 

import sympy as sp
from base import Base
from sympy.utilities.lambdify import lambdify
import os
import logging

class LogLikelihood(Base):
    def __init__(self, expressions):
        """
        Initializes the Likelihood calculator, supporting either a single expression or a list of expressions.

        Args:
            expressions (str or list): The Likelihood expression(s) to be evaluated.
                                       This can either be a single expression string or a list containing multiple expression strings.
        """
        # Convert a single expression into a list for uniform processing
        from inner_func import update_funcs, update_const
        self.named_expressions = []
        
        for expr_dict in expressions:
            if isinstance(expr_dict, dict) and 'name' in expr_dict and 'expression' in expr_dict:
                # Parse the 'expression' value to create a sympy expression, with updated functions available
                sympy_expr = sp.sympify(expr_dict['expression'], locals=update_funcs({}))
                # Append the name and sympy expression as a tuple
                self.named_expressions.append((expr_dict['name'], sympy_expr))

        self.variables = set().union(*[expr.free_symbols for _, expr in self.named_expressions])
        self.logger = None
        self.custom_functions = {}
        self.custom_functions = update_funcs(self.custom_functions)
        self.constants = {}
        self.constants = update_const(self.constants)
        self.values = {}

    def update_funcs(self, funcs):
        self.custom_functions.update(funcs)
        # self.logger.info(f"Jarvis-HEP likelihood now support the following inner functions -> \n{self.custom_functions.keys()}")

    def calculate(self, values, sample_info):
        """
        Calculates the Likelihood value based on the provided values of variables. If multiple expressions are provided, their sum is computed.

        Args:
            values (dict): A dictionary containing the values for all variables in the expressions, e.g., {'a': 1, 'b': 2, 'c': 3}.

        Returns:
            float: The calculated Likelihood value. If there are multiple expressions, the sum of their values is returned.
        """
        self.update_logger(sample_info)
        try: 
            # Ensure that values for all variables are provided
            assert {str(var) for var in self.variables}.issubset(values.keys()), "Not all variables have values provided."

            total_loglikelihood = 0. 
            self.values = {}
            for expr_dict in self.named_expressions:
                expr = expr_dict[1]
                name = expr_dict[0]
                num_expr = lambdify([str(var) for var in expr.free_symbols], expr, modules=[self.custom_functions, "numpy"])
                varis = set(expr.free_symbols)
                symbol_values_strs = {str(key): value for key, value in values.items() if key in {str(var) for var in varis}}
                likelihood = num_expr(**symbol_values_strs)
                self.childlogger.info(f"Evaluating   {name}: \n   expression \t-> {str(expr)} \n   with input \t-> {symbol_values_strs} \n   Output \t\t-> {likelihood}")
                self.values[name] = likelihood
                total_loglikelihood += likelihood

            # self.childlogger.warning(f"\t Total LogLikelihood -> \n\t LogL: {total_loglikelihood}")
            self.values['LogL'] = total_loglikelihood
            return self.values 
        except Exception as exc:
            print(exc)
            return {"LogL": float(-sp.core.numbers.Infinity())}


    def update_logger(self, sample_info):
        logger_name = f"Sample@{sample_info['uuid']} <Likelihood>"
        if not os.path.exists(sample_info['save_dir']):
            os.makedirs(sample_info['save_dir'])
        self.childlogger = logging.getLogger(logger_name)
        self.childlogger.setLevel(logging.DEBUG)
        self.childlogger.propagate = False

        # self.path['run_log_file'] = os.path.join(sample_info['save_dir'], "Sample_running.log")
        # self.path['run_log_file'] = sample_info['run_log']
        formatter = logging.Formatter("\n·•· %(name)s\n\t -> %(asctime)s - %(levelname)s >>> \n%(message)s")
        file_handler = logging.FileHandler(sample_info['run_log'], mode='a')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)

        jlog_handler = logging.FileHandler(sample_info['jarvis_log'])
        jlog_handler.setFormatter(formatter)
        jlog_handler.setLevel(logging.WARNING)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.WARNING)
        console_handler.setFormatter(formatter)

        self.childlogger.addHandler(console_handler)
        self.childlogger.addHandler(file_handler)