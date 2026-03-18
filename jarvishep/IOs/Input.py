#!/usr/bin/env python3

import os
import json
import numpy
import pyslha
import sympy as sp
from sympy.utilities.lambdify import lambdify

from jarvishep.IOs.IOs import InputFile
from jarvishep.inner_func import build_expression_context, update_const


def _evaluate_expression(expression, param_values, parse_locals, numeric_modules):
    expr = sp.sympify(expression, locals=parse_locals)
    symbol_names = [str(par) for par in expr.free_symbols]
    symbol_name_set = set(symbol_names)
    num_expr = lambdify(symbol_names, expr, modules=[numeric_modules, "numpy"])
    symbol_values_strs = {
        str(key): value
        for key, value in param_values.items()
        if str(key) in symbol_name_set
    }
    value = num_expr(**symbol_values_strs)
    return expr, symbol_values_strs, value

class SLHAInputFile(InputFile):
    async def write(self, param_values):
        self.path = self.decode_path(self.path)
        return await self.io_run_blocking(self._write_sync, param_values)

    def _write_sync(self, param_values):
        self.logger.info(f"Start writing the input file -> {self.path}")
        if not self.path:
            self.logger.error(f"Input file {self.path} not found")

        self.sync_make_dirs(os.path.dirname(self.path), exist_ok=True)
        observables = {}
        parse_locals, numeric_modules = build_expression_context(
            funcs=dict(self.funcs or {}),
            consts=update_const({}),
        )

        try:
            content = self.sync_read_text(self.path)
            for action in self.variables:
                if action["type"] == "Replace":
                    for var in action["variables"]:
                        if "expression" not in var:
                            value = param_values.get(var["name"], "MISSING_VALUE")
                        else:
                            expr, symbol_values_strs, value = _evaluate_expression(
                                var["expression"],
                                param_values,
                                parse_locals,
                                numeric_modules,
                            )
                            observables[var["name"]] = value
                            self.logger.info(
                                f"Evaluating: expression \n\t-> {expr} \n    with input \t -> [ {', '.join(['{}: {}, '.format(kk, vv) for kk, vv in symbol_values_strs.items()])}] \n    Output \t\t-> {value}"
                            )
                        placeholder = var["placeholder"]
                        if value != "MISSING_VALUE":
                            value = f"{float(value):.8E}"
                        content = content.replace(placeholder, value)
                elif action["type"] == "SLHA":
                    slha_content = pyslha.readSLHA(content)
                    for var in action["variables"]:
                        if "block" not in var:
                            continue
                        if "expression" not in var:
                            value = param_values.get(var["name"], "MISSING_VALUE")
                        else:
                            expr, symbol_values_strs, value = _evaluate_expression(
                                var["expression"],
                                param_values,
                                parse_locals,
                                numeric_modules,
                            )
                            observables[var["name"]] = value
                            self.logger.info(
                                f"Evaluating: expression \n\t-> {expr} \n    with input \t -> [{', '.join(['{} : {}'.format(kk, vv) for kk, vv in symbol_values_strs.items()])}] \n    Output \t\t-> {value}"
                            )
                        if value != "MISSING_VALUE":
                            value = f"{float(value):.8E}"
                        if isinstance(var["entry"], int):
                            slha_content.blocks[var["block"]][var["entry"]] = value
                        elif isinstance(var["entry"], tuple):
                            slha_content.blocks[var["block"]][tuple(var["entry"])] = value
                        else:
                            self.logger.warning(
                                f"Invalid Entry type: {type(var['entry'])}. Entry must be an integer or a tuple of integers."
                            )
                    content = pyslha.writeSLHA(slha_content, ignorenobr=True)
                elif action["type"] == "File":
                    source_path = param_values.get(action["source"], None)
                    if source_path and self.sync_exists(source_path):
                        content = self.sync_read_text(source_path)
                        break
                    self.logger.warning(
                        f"Input source file is not found: -> \n\t{action['source']} \n\t{source_path}"
                    )

            self.sync_write_text(self.path, content)
            if self.save:
                target = os.path.join(
                    self.sample_save_dir,
                    f"{os.path.basename(self.path)}@{self.module}",
                )
                self.sync_write_text(target, content)
                observables[self.name] = target

            self.logger.info(f"Finish writing the input file -> {self.path}")
            self.logger = None
            return observables
        except Exception as e:
            self.logger.error(f"Error writing SLHA input file '{self.name}': {e}")
        self.logger = None

class JsonInputFile(InputFile):
    @staticmethod
    def _to_json_compatible(value):
        if isinstance(value, numpy.generic):
            return value.item()
        if isinstance(value, numpy.ndarray):
            return value.tolist()
        if isinstance(value, dict):
            return {k: JsonInputFile._to_json_compatible(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [JsonInputFile._to_json_compatible(v) for v in value]
        return value

    async def write(self, param_values):
        self.path = self.decode_path(self.path)
        return await self.io_run_blocking(self._write_sync, param_values)

    def _write_sync(self, param_values):
        self.logger.info(f"Start writing the input file -> {self.path}")
        observables = {}
        parse_locals, numeric_modules = build_expression_context(
            funcs=dict(self.funcs or {}),
            consts=update_const({}),
        )

        try:
            content = self.sync_read_text(self.path)
            data_to_write = json.loads(content)
        except FileNotFoundError:
            self.logger.error(f"File not found: {self.path}")
            data_to_write = {}
        except json.JSONDecodeError:
            self.logger.error(f"Error decoding JSON from file: {self.path}")
            data_to_write = {}

        for action in self.variables:
            if action["type"] != "Dump":
                continue
            for var in action["variables"]:
                if "expression" not in var:
                    value = param_values.get(var["name"], "MISSING_VALUE")
                else:
                    expr, symbol_values_strs, value = _evaluate_expression(
                        var["expression"],
                        param_values,
                        parse_locals,
                        numeric_modules,
                    )
                    observables[var["name"]] = value
                    self.logger.info(
                        f"Evaluating: expression \n\t-> {expr} \n    with input \t -> [{', '.join(['{} : {}'.format(kk, vv) for kk, vv in symbol_values_strs.items()])}] \n    Output \t\t-> {value}"
                    )

                if "entry" not in var:
                    data_to_write.update({var["name"]: value})
                else:
                    self.update_json_by_entry(data_to_write, var["entry"], value)

        try:
            json_ready = self._to_json_compatible(data_to_write)
            json_str = json.dumps(json_ready, indent=4)
            self.sync_write_text(self.path, json_str)
        except Exception as e:
            self.logger.error(f"Error writing Json input file '{self.name}': {e}")

        self.logger = None
        return observables


    def update_json_by_entry(self, json_dict, entry, new_value):
        parts = entry.split('.')
        current = json_dict
        for part in parts[:-1]: 
            if part not in current:
                current[part] = {} 
            current = current[part]

        current[parts[-1]] = new_value
