import math
from typing import Any, Dict, Iterable, Mapping, Set, Tuple

import sympy as sp

from inner_func import update_funcs, update_const


def build_default_locals() -> Dict[str, Any]:
    d: Dict[str, Any] = {}
    d = update_funcs(d)
    d = update_const(d)
    return d


class NuisanceExpressionRegistry:
    """Single-expression compiler/evaluator.

    - set_config(name, expression) configures + compiles exactly ONE expression
    - can_eval(keys) checks whether all deps are available
    - eval(values) evaluates using dict input

    Design:
    - SymPy handles parsing.
    - deps are stored as a tuple in deterministic order (sorted by symbol name).
    - eval input is a mapping/dict; positional args are built internally from deps.
    """

    def __init__(self):
        self.name: str = ""
        self.expression: str = ""
        self._deps: Tuple[str, ...] = ()
        self._locals: Dict[str, Any] = build_default_locals()
        self._fn = None  # compiled callable

    def set_config(self, name: str, expression: str) -> None:
        """Configure and compile exactly one expression.
        User API: only requires `name` and `expression`.
        """
        self.name = str(name)
        self.expression = str(expression)
        self._compile()

    @property
    def deps(self) -> Tuple[str, ...]:
        return self._deps

    @property
    def fn(self):
        return self._fn

    def can_eval(self, available_keys: Iterable[str]) -> Tuple[bool, Set[str]]:
        avail = set(available_keys)
        missing = set(self.deps) - avail
        return (len(missing) == 0), missing

    def eval(self, values: Mapping[str, Any]) -> Any:
        args = [values[k] for k in self.deps]
        return self.fn(*args)

    # ---------- internal ----------
    def _compile(self) -> None:
        sym_expr = sp.sympify(self.expression, locals=dict(self._locals))
        self._deps = tuple(sorted(str(s) for s in sym_expr.free_symbols))

        sym_locals: Dict[str, Any] = {k: sp.Symbol(k) for k in self._deps}
        sym_locals.update(self._locals)

        expr = sp.sympify(self.expression, locals=sym_locals)
        symbols = [sym_locals[k] for k in self._deps]

        self._fn = sp.lambdify(symbols, expr, modules=[self._locals, math])