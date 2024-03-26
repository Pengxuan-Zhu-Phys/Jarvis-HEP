#!/usr/bin/env python3 

import sympy as sp
from base import Base

class Likelihood(Base):
    def __init__(self, expressions):
        """
        Initializes the Likelihood calculator, supporting either a single expression or a list of expressions.

        Args:
            expressions (str or list): The Likelihood expression(s) to be evaluated.
                                       This can either be a single expression string or a list containing multiple expression strings.
        """
        # Convert a single expression into a list for uniform processing
        if isinstance(expressions, str):
            expressions = [expressions]  
        # Parse all expressions, converting each string into a sympy expression
        self.expressions = [sp.sympify(expr) for expr in expressions]
        # Extract variables from all expressions
        self.variables = set().union(*[expr.free_symbols for expr in self.expressions])

    def calculate(self, values):
        """
        Calculates the Likelihood value based on the provided values of variables. If multiple expressions are provided, their sum is computed.

        Args:
            values (dict): A dictionary containing the values for all variables in the expressions, e.g., {'a': 1, 'b': 2, 'c': 3}.

        Returns:
            float: The calculated Likelihood value. If there are multiple expressions, the sum of their values is returned.
        """
        # Ensure that values for all variables are provided
        assert self.variables.issubset(values.keys()), "Not all variables have values provided."
        # Calculate the value for each expression and sum them up
        total_likelihood = sum(expr.evalf(subs=values) for expr in self.expressions)
        return total_likelihood

