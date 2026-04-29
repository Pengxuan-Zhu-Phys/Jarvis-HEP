#!/usr/bin/env python3

from copy import deepcopy
import os
import json
import xslha
import pyslha
from jarvishep.IOs.IOs import OutputFile

class SLHAOutputFile(OutputFile):
    """
        A class designed to asynchronously read, process, and optionally save SLHA files based on specified observables.

        Inherits from OutputFile and utilizes asynchronous file operations to efficiently process SLHA files without blocking the event loop. This class is particularly useful in environments where responsiveness and non-blocking I/O operations are critical.

        Attributes:
            Inherits all attributes from the OutputFile class such as path, logger, and save flag.
            variables (list of dicts): Each dictionary specifies an observable to extract from the SLHA file. The keys in the dictionary include:
                - name: A string representing the unique identifier of the observable.
                - block: The SLHA block from which the value should be extracted.
                - entry: The entry (or entries) within the block identifying the specific value or values to be extracted. Can be an integer for single values or a list for complex entries such as decay branching ratios.
    """

    async def read(self):
        self.path = self.decode_path(self.path)
        return await self.io_run_blocking(self._read_sync)

    def _read_sync(self):
        self.logger.info(f"Start reading the output file -> {self.path}")
        observables = {}
        try:
            source = self.sync_read_text(self.path)
            content = pyslha.readSLHA(source)

            for var in self.variables:
                if var["block"] == "DECAY":
                    if isinstance(var["entry"], int):
                        observables[var["name"]] = float(
                            content.decays[var["entry"]].__dict__["totalwidth"]
                        )
                    elif isinstance(var["entry"], list):
                        decays = content.decays[var["entry"][0]].__dict__["decays"]
                        value = 0.0
                        for decay in decays:
                            if set(var["entry"][1:]) == set(decay.ids):
                                value = decay.br
                                break
                        observables[var["name"]] = value
                    else:
                        self.logger.error(f"Unsupport decay entry {var['entry']}")
                elif content.blocks[var["block"]]:
                    if isinstance(var["entry"], int):
                        observables[var["name"]] = content.blocks[var["block"]][var["entry"]]
                    elif isinstance(var["entry"], list):
                        observables[var["name"]] = content.blocks[var["block"]][var["entry"]]
                    else:
                        self.logger.error(f"Unsupported block entry {var['entry']}")
                else:
                    self.logger.error(f"Unsupported SLHA read item {var}")

            if self.save:
                target = os.path.join(
                    self.sample_save_dir,
                    f"{os.path.basename(self.path)}@{self.module}",
                )
                self.sync_write_text(target, source)
                observables[self.name] = target
            else:
                target_path = os.path.join(self.sample_save_dir, ".temp")
                self.sync_make_dirs(target_path, exist_ok=True)
                target = os.path.join(
                    target_path,
                    f"{os.path.basename(self.path)}@{self.module}",
                )
                self.sync_write_text(target, source)
                observables[self.name] = target

            self.logger.info(f"Finish reading the input file -> {self.path}")
            return observables
        except Exception as e:
            self.logger.error(f"Error reading SLHA input file '{self.name}': {e}")
            return observables


class JsonOutputFile(OutputFile):
    async def read(self):
        self.path = self.decode_path(self.path)
        return await self.io_run_blocking(self._read_sync)

    def _read_sync(self):
        self.logger.info(f"Start reading the output file -> {self.path}")
        observables = {}
        content = None
        source = None

        try:
            try:
                source = self.sync_read_text(self.path)
                content = json.loads(source)
            except FileNotFoundError:
                self.logger.error(f"File not found: {self.path}")
            except json.JSONDecodeError:
                self.logger.error(f"Error decoding JSON from file: {self.path}")

            for var in self.variables:
                if "entry" in var:
                    observables[var["name"]] = self.read_json_value_by_entry(content, var["entry"])
                else:
                    observables[var["name"]] = content[var["name"]]

            if self.save:
                target = os.path.join(
                    self.sample_save_dir,
                    f"{os.path.basename(self.path)}@{self.module}",
                )
                self.sync_write_text(target, source)
                observables[self.name] = target
            else:
                target_path = os.path.join(self.sample_save_dir, ".temp")
                self.sync_make_dirs(target_path, exist_ok=True)
                target = os.path.join(
                    target_path,
                    f"{os.path.basename(self.path)}@{self.module}",
                )
                self.sync_write_text(target, source)
                observables[self.name] = target

            self.logger.info(f"Finish reading the input file -> {self.path}")
            return observables
        except Exception as e:
            self.logger.error(f"Error reading SLHA input file '{self.name}': {e}")
            return observables
        
    def read_json_value_by_entry(self, json_dict, entry):
        parts = entry.split('.')
        current = json_dict
        for part in parts:
            if part not in current:
                return None  # 如果路径的某部分不存在，则返回 None
            current = current[part]
    
        return current
    
    
class xSLHAOutputFile(OutputFile):
    """
        A class designed to asynchronously read, process, and optionally save SLHA files based on specified observables.

        Inherits from OutputFile and utilizes asynchronous file operations to efficiently process SLHA files without blocking the event loop. This class is particularly useful in environments where responsiveness and non-blocking I/O operations are critical.

        Attributes:
            Inherits all attributes from the OutputFile class such as path, logger, and save flag.
            variables (list of dicts): Each dictionary specifies an observable to extract from the SLHA file. The keys in the dictionary include:
                - name: A string representing the unique identifier of the observable.
                - block: The SLHA block from which the value should be extracted.
                - entry: The entry (or entries) within the block identifying the specific value or values to be extracted. Can be an integer for single values or a list for complex entries such as decay branching ratios.
    """

    async def read(self):
        self.path = self.decode_path(self.path)
        return await self.io_run_blocking(self._read_sync)

    def _read_sync(self):
        self.logger.info(f"Start reading the output file -> {self.path}")
        observables = {}
        try:
            source = self.sync_read_text(self.path)
            content = xslha.read(self.path)
            for var in self.variables:
                if var["block"] == "DECAY":
                    if isinstance(var["entry"], int):
                        observables[var["name"]] = content.widths[var["entry"]]
                    elif isinstance(var["entry"], list):
                        decays = content.br[var["entry"][0]]
                        dkentry = tuple(sorted(deepcopy(var["entry"][1:])))
                        observables[var["name"]] = decays.get(dkentry, 0.0)
                    else:
                        self.logger.error(f"Unsupport decay entry {var['entry']}")
                elif var["block"] == "DECAY1L":
                    if isinstance(var["entry"], int):
                        observables[var["name"]] = content.widths1L[var["entry"]]
                    elif isinstance(var["entry"], list):
                        decays = content.br1L[var["entry"][0]]
                        dkentry = tuple(sorted(deepcopy(var["entry"][1:])))
                        observables[var["name"]] = decays.get(dkentry, 0.0)
                    else:
                        self.logger.error(f"Unsupport decay entry {var['entry']}")
                elif var["block"] == "HIGGSBOUNDS":
                    if isinstance(var["entry"], list):
                        blkentry = ",".join(list(map(str, var["entry"])))
                        if blkentry in content.blocks["HIGGSBOUNDS"]:
                            observables[var["name"]] = content.blocks["HIGGSBOUNDS"][blkentry]
                        else:
                            observables[var["name"]] = 0.0
                            self.logger.info(
                                f"Entry {blkentry} in HIGGSBOUNDS unfound, setted as 0."
                            )
                    else:
                        self.logger.error(f"Unsupport entry {var['entry']} for block HIGGSBOUNDS")
                elif content.blocks[var["block"].upper()]:
                    blk = var["block"].upper()
                    if isinstance(var["entry"], int):
                        value = content.blocks[blk][str(var["entry"])]
                    elif isinstance(var["entry"], list):
                        value = content.blocks[blk][",".join(list(map(str, var["entry"])))]
                    else:
                        self.logger.error(f"Unsupported block entry {var['entry']}")
                        continue
                    if not isinstance(value, list):
                        observables[var["name"]] = value
                    else:
                        observables[f"Re[{var['name']}]"] = value[0]
                        observables[f"Im[{var['name']}]"] = value[1]
                else:
                    self.logger.error(f"Unsupported SLHA read item {var}")
            if self.save:
                target = os.path.join(
                    self.sample_save_dir,
                    f"{os.path.basename(self.path)}@{self.module}",
                )
                self.sync_write_text(target, source)
                observables[self.name] = target
            else:
                target_path = os.path.join(self.sample_save_dir, ".temp")
                self.sync_make_dirs(target_path, exist_ok=True)
                target = os.path.join(
                    target_path,
                    f"{os.path.basename(self.path)}@{self.module}",
                )
                self.sync_write_text(target, source)
                observables[self.name] = target

            self.logger.info(f"Finish reading the input file -> {self.path}")
            return observables
        except Exception as e:
            self.logger.error(f"Error reading SLHA input file '{self.name}': {e}")
            return observables


class FileOutput(OutputFile):
    async def read(self):
        self.path = self.decode_path(self.path)
        return await self.io_run_blocking(self._read_sync)

    def _read_sync(self):
        self.logger.info(f"Start reading the output file -> {self.path}")
        observables = {}
        try:
            try:
                source = self.sync_read_text(self.path)
            except FileNotFoundError:
                self.logger.error(f"File not found: {self.path}")
                source = None
            except json.JSONDecodeError:
                self.logger.error(f"Error decoding JSON from file: {self.path}")
                source = None

            if self.save:
                target = os.path.join(
                    self.sample_save_dir,
                    f"{os.path.basename(self.path)}@{self.module}",
                )
                self.sync_write_text(target, source)
                observables[self.name] = target
            return observables
        except Exception as e:
            self.logger.error(f"Error reading SLHA input file '{self.name}': {e}")
            return observables
