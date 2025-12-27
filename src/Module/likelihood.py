#!/usr/bin/env python3 

import sympy as sp
from base import Base
from sympy.utilities.lambdify import lambdify
import os
import logging
from loguru import logger
from numpy import inf
import numpy as np
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
        self.expressions = expressions
        for expr_dict in expressions:
            if isinstance(expr_dict, dict) and 'name' in expr_dict and 'expression' in expr_dict:
                # Parse the 'expression' value to create a sympy expression, with updated functions available
                sympy_expr = sp.sympify(expr_dict['expression'], locals=update_funcs({}))
                # Append the name and sympy expression as a tuple
                self.named_expressions.append((expr_dict['name'], sympy_expr))

        self.variables = set().union(*[expr.free_symbols for _, expr in self.named_expressions])
        self.logger = logger.bind(module="Jarvis-HEP.LogLikelihood", to_console=True, Jarvis=True)
        self.custom_functions = {}
        self.custom_functions = update_funcs(self.custom_functions)
        self.constants = {}
        self.constants = update_const(self.constants)
        self.values = {}

        # Cache compiled numerical callables for each expression.
        # Each item is (name, sympy_expr, var_names, var_name_set, num_expr)
        self._compiled_expressions = []
        self._compile_expressions()

    def update_funcs(self, funcs):
        self.custom_functions.update(funcs)
        # Rebuild compiled expressions so newly-added functions are available.
        self._compile_expressions()
        # self.logger.info(f"Jarvis-HEP likelihood now support the following inner functions -> \n{self.custom_functions.keys()}")

    def _compile_expressions(self):
        """Compile each sympy expression to a numerical callable once.

        We keep kwargs-style evaluation (num_expr(**dict)).
        """
        compiled = []
        for name, expr in self.named_expressions:
            var_names = [str(var) for var in expr.free_symbols]
            var_name_set = set(var_names)
            num_expr = lambdify(var_names, expr, modules=[self.custom_functions, "numpy"])
            compiled.append((name, expr, var_names, var_name_set, num_expr))
        self._compiled_expressions = compiled

    @staticmethod
    def format_summary(values: dict,
                       key_width: int | None = None,
                       value_width: int = 60,
                       float_precision: int = 6) -> str:
        """Fast, pandas-free summary formatter.

        Produces an aligned two-column listing similar to `pandas.Series(values).to_string()`.
        - Keys are left-aligned.
        - Values are right-aligned.
        - Long strings are truncated with '...'.
        - Floats are formatted with up to `float_precision` decimal places, trimming trailing zeros.
        """

        def _format_value(v):
            # Keep None explicit
            if v is None:
                s = "None"
            # Numpy scalars -> python scalars
            elif isinstance(v, (np.generic,)):
                s = str(v.item())
            elif isinstance(v, bool):
                s = "True" if v else "False"
            elif isinstance(v, (int,)):
                s = str(v)
            elif isinstance(v, (float,)):
                if np.isnan(v):
                    s = "nan"
                elif np.isposinf(v):
                    s = "inf"
                elif np.isneginf(v):
                    s = "-inf"
                else:
                    av = abs(v)
                    # Use scientific notation for very small/large magnitudes to avoid rounding to 0
                    if (av != 0.0 and av < 1e-4) or av >= 1e6:
                        s = f"{v:.{float_precision}e}"
                    else:
                        s = f"{v:.{float_precision}f}".rstrip("0").rstrip(".")
            else:
                s = str(v)

            # Truncate overly long strings similar to pandas display
            if len(s) > value_width:
                if value_width <= 3:
                    s = s[:value_width]
                else:
                    s = s[: value_width - 3] + "..."
            return s

        # Preserve insertion order of `values` (dict preserves insertion order in Python 3.7+)
        keys = list(values.keys())
        if key_width is None:
            key_width = max((len(str(k)) for k in keys), default=1)
            # Keep it reasonable; very long keys reduce readability
            key_width = min(key_width, 60)

        lines = []
        for k in keys:
            ks = str(k)
            if len(ks) > key_width:
                # Truncate long keys to keep alignment stable
                if key_width <= 3:
                    ks = ks[:key_width]
                else:
                    ks = ks[: key_width - 3] + "..."
            vs = _format_value(values[k])
            # Match the look: key left, value right with plenty of spacing
            lines.append(f"{ks:<{key_width}}  {vs:>{value_width}}")
        return "\n".join(lines)

    def calculate4dnn(self, df):
        """
        Calculate total log-likelihood for each row in a DataFrame.

        Args:
            df (pandas.DataFrame): Input data; columns should include all variables used in the expressions.

        Returns:
            numpy.ndarray: Array of total log-likelihood values, one per DataFrame row.
        """
        results = []
        for _, row in df.iterrows():
            try:
                total = 0.0
                for name, expr in self.named_expressions:
                    var_names = [str(var) for var in expr.free_symbols]
                    symbol_values = {var: row[var] for var in var_names if var in row}
                    num_expr = lambdify(var_names, expr, modules=[self.custom_functions, "numpy"])
                    total += float(num_expr(**symbol_values))
                results.append(total)
            except Exception:
                results.append(-np.inf)
        return np.array(results)

    def calculate(self, values, sample_info):
        """
        Calculates the Likelihood value based on the provided values of variables. If multiple expressions are provided, their sum is computed.

        Args:
            values (dict): A dictionary containing the values for all variables in the expressions, e.g., {'a': 1, 'b': 2, 'c': 3}.

        Returns:
            float: The calculated Likelihood value. If there are multiple expressions, the sum of their values is returned.
        """
        try: 
            slogger = self.update_logger(sample_info)
            # Ensure that values for all variables are provided
            assert {str(var) for var in self.variables}.issubset(values.keys()), "Not all variables have values provided."
            total_loglikelihood = 0. 
            self.values = {}
            for name, expr, _var_names, var_name_set, num_expr in self._compiled_expressions:
                symbol_values_strs = {str(key): value for key, value in values.items() if key in var_name_set}
                likelihood = float(num_expr(**symbol_values_strs))
                slogger.info(
                    f"Evaluating   {name}: \n   expression \t-> {str(expr)} \n   with input \t-> [{', '.join(['{} : {}'.format(kk, vv) for kk, vv in symbol_values_strs.items()])}] \n   Output \t\t-> {likelihood}"
                )
                self.values[name] = likelihood
                total_loglikelihood += likelihood

            # self.childlogger.warning(f"\t Total LogLikelihood -> \n\t LogL: {total_loglikelihood}")
            self.values['LogL'] = total_loglikelihood
            values.update(self.values)
            slogger.info(
                f"Sample SUMMARY\n============================================================================\n{LogLikelihood.format_summary(values)}\n============================================================================"
            )

            return self.values 
        except Exception as exc:
            self.logger.warning(exc)
            # return {"LogL": float(-sp.core.numbers.Infinity())}
            return {"LogL": - inf}


    def custom_format(record):
        module = record["extra"].get("module", "No module")
        if "raw" in record["extra"]:
            return "{message}"
        else:
            return f"\n·•· <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{record['message']}</level>"


    def update_logger(self, sample_info):
        logger_name = f"{sample_info['logger_name']} (Likelihood)"
        sample_logger = logger.bind(module=logger_name, to_console=True, Jarvis=True)
        return sample_logger

    def __deepcopy__(self, memo):
        # Create a new instance
        from copy import deepcopy
        copied = LogLikelihood(deepcopy(self.expressions, memo))

        # Copy deepcopy-safe attributes
        copied.variables = deepcopy(self.variables, memo)
        copied.constants = deepcopy(self.constants, memo)
        copied.values = deepcopy(self.values, memo)

        # Functions and loggers cannot be copied directly, so reinitialize them
        copied.custom_functions = self.custom_functions  # Functions remain shared
        copied.logger = None
        copied.childlogger = None
        copied.childhandler = None
        return copied