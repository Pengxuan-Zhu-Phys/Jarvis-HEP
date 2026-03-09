#!/usr/bin/env python3 
import asyncio
import os

import sympy as sp
from loguru import logger

from jarvishep.Module.module import Module
from jarvishep.async_subprocess import SubprocessExecutionError, SubprocessJob


class CalculatorModule(Module):
    def __init__(self, name, config):
        super().__init__(name)
        if config["modes"]:
            self.modes = True
            self.analyze_config_multi()
        else:
            self.config             = config 
            self.type               = "Calculator"
            self.required_modules   = config["required_modules"]
            self.clone_shadow       = config["clone_shadow"]
            self.installation       = config["installation"]
            self.initialization     = config["initialization"]
            self.execution          = config["execution"]
            self.input              = config["execution"].get("input", [])
            self.output             = config["execution"].get("output", [])
            self.basepath           = config['path']
            self.handlers           = {}
            self.modes              = False
            self.is_installed       = False
            self.installation_event = None
            self.is_busy            = False
            self.PackID             = None
            self.run_log_file       = None
            self.sample_info        = {}
            self._funcs             = {}
            self.subprocess_scheduler = None
            self._command_counter   = 0
            self.analyze_config()

    def assign_ID(self, PackID):
        self.PackID = PackID 

    def analyze_config(self):
        # get the variables for inputs and outputs, prepare for the workflow chart 
        from jarvishep.inner_func import build_expression_context, update_const, update_funcs
        parse_locals, _ = build_expression_context(
            funcs=update_funcs({}),
            consts=update_const({}),
        )
        for ipf in self.input:
            if ipf['type'] == "SLHA":
                ipf['variables'] = {}
                for act in ipf['actions']:
                    if act['type'] == "Replace":
                        for var in act['variables']:
                            if "expression" in var.keys():
                                expr = sp.sympify(var['expression'], locals=parse_locals)
                                varis = set(expr.free_symbols)
                                for vv in varis:
                                    self.inputs[str(vv)] = None
                                    # ipf['variables'][str(vv)] = var 
                                var.update({"inc": varis})
                                ipf['variables'][var['name']] = var 
                            else:
                                self.inputs[var['name']] = None
                                ipf['variables'][var['name']] = var
                    elif act['type'] == "SLHA":
                        for var in act['variables']:
                            if "expression" in var.keys():
                                expr = sp.sympify(var['expression'], locals=parse_locals)
                                varis = set(expr.free_symbols)
                                for vv in varis:
                                    self.inputs[str(vv)] = None
                                    # ipf['variables'][str(vv)] = var 
                                var.update({"inc": varis})
                                ipf['variables'][var['name']] = var 
                            else:
                                self.inputs[var['name']] = None
                                ipf['variables'][var['name']] = var
                    elif act['type'] == "File":
                        self.inputs[act['source']] = None
                        ipf['variables'][act['source']] = {"name": act['source']}
                # print(ipf['variables'])

            elif ipf['type'] == "Json":
                ipf['variables'] = {}
                for act in ipf['actions']:
                    if act['type'] == "Dump":
                        for var in act['variables']:
                            if "expression" in var.keys():
                                expr = sp.sympify(var['expression'], locals=parse_locals)
                                varis = set(expr.free_symbols)
                                for vv in varis:
                                    self.inputs[str(vv)] = None
                                    # ipf['variables'][str(vv)] = var 
                                var.update({"inc": varis})
                                ipf['variables'][var['name']] = var 
                            # print(var.keys())
                            else:
                                ipf['variables'][var['name']] = var 
                                self.inputs[var['name']] = None
                # pprint(ipf['variables'])

            else:
                for ipv in ipf['variables']:
                    self.inputs[ipv['name']] = None
        # pprint(self.inputs)
        for opf in self.output:
            # self.outputs[opf['name']] = None
            if opf['type'] == "File":
                pass 
            else:
                for opv in opf['variables']:
                    self.outputs[opv['name']] = None

    @property
    def funcs(self):
        return self._funcs

    def set_subprocess_scheduler(self, scheduler):
        self.subprocess_scheduler = scheduler

    def analyze_config_multi(self):
        pass 

    def custom_format(record):
        module = record["extra"].get("module", "No module")
        if "raw" in record["extra"]:
            return "{message}"
        else:
            return f"\n·•· <cyan>{module}</cyan> \n\t-> <green>{record['time']:MM-DD HH:mm:ss.SSS}</green> - [<level>{record['level']}</level>] >>> \n<level>{{message}}</level>"
    
    def create_basic_logger(self):
        logger_name = f"{self.name}-{self.PackID}"
        def filte_func(record):
            return record['extra']['module'] == logger_name
        
        self.basepath = self.decode_shadow_path(self.basepath)
        if not os.path.exists(self.basepath):
            os.makedirs(self.basepath)

        install_file_log = os.path.join(self.basepath, f"Installation_{self.name}-{self.PackID}.log")
        self.logger = logger.bind(
            module=logger_name,
            to_console=True,
            Jarvis=True,
            _log_domain="jarvis_hep",
        )
        install_handler = self.logger.add(install_file_log, format=CalculatorModule.custom_format, level="DEBUG", rotation=None, retention=None, filter=filte_func )
        self.handlers['install'] = install_handler

    def update_sample_logger(self, sample_info):
        # self.close_sample_logger()
        logger_name = f"{sample_info['logger_name']} ({self.name}-No.{self.PackID})"

        self.logger = sample_info['logger'].bind(
            module=logger_name,
            to_console=True,
            Jarvis=True,
            _log_domain="jarvis_hep",
        )
        self.logger.info("Module load instance and logger is correctly set!")

    def close_sample_logger(self):
        
        """Ensure logging processor is properly shut down after task completion"""
        # if 'sample' in self.handlers:
            # self.logger.info("Closing calculator logging handler")
            # sample_handler = self.handlers['sample']
            # self.logger.remove(sample_handler)
            # del self.handlers['sample']  
        self.logger = None 

    def _next_command_index(self) -> int:
        self._command_counter += 1
        return self._command_counter

    def _build_command_meta(self, stage: str, command_index: int, command: dict) -> dict:
        sample_uuid = None
        if isinstance(self.sample_info, dict):
            sample_uuid = self.sample_info.get("uuid")
        return {
            "module": self.name,
            "pack_id": self.PackID,
            "stage": stage,
            "command_index": int(command_index),
            "sample_uuid": sample_uuid,
            "cwd": command.get("cwd"),
        }

    def _resolve_sample_runtime_tokens(self, text: str, *, stage: str, field: str) -> str:
        """Resolve runtime-only command tokens.

        `@SampleID` is resolved from `sample_info['uuid']` only during task-run stages.
        Install stage intentionally skips this replacement.
        """
        if text is None:
            return ""
        raw = str(text)
        if stage == "install" or "@SampleID" not in raw:
            return raw

        sample_uuid = None
        if isinstance(self.sample_info, dict):
            sample_uuid = self.sample_info.get("uuid")
        if sample_uuid is None:
            raise RuntimeError(
                f"@SampleID requires sample_info['uuid'] during runtime stage '{stage}' for field '{field}'"
            )
        return raw.replace("@SampleID", str(sample_uuid))

    async def _run_command_local(self, command: dict, stage: str, command_index: int):
        process = await asyncio.create_subprocess_shell(
            command["cmd"],
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=command["cwd"],
            start_new_session=True,
        )

        async def _drain(stream, stream_name: str):
            total = 0
            text_buffer = ""
            while True:
                chunk = await stream.read(65536)
                if not chunk:
                    break
                total += len(chunk)
                if self.logger is not None:
                    text_buffer += chunk.decode(errors="replace")
                    while "\n" in text_buffer:
                        line, text_buffer = text_buffer.split("\n", 1)
                        line = line.rstrip("\r")
                        if line:
                            self.logger.info(
                                f"[{stage}#{command_index:05}][{stream_name}] {line}"
                            )
            if self.logger is not None:
                tail = text_buffer.rstrip("\r")
                if tail:
                    self.logger.info(
                        f"[{stage}#{command_index:05}][{stream_name}] {tail}"
                    )
            return total

        out_task = asyncio.create_task(_drain(process.stdout, "stdout"))
        err_task = asyncio.create_task(_drain(process.stderr, "stderr"))
        rc = await process.wait()
        stdout_bytes, stderr_bytes = await asyncio.gather(out_task, err_task)
        if self.logger is not None:
            self.logger.info(
                "Command done [{}#{:05}] rc={} out={}B err={}B".format(
                    stage,
                    command_index,
                    rc,
                    int(stdout_bytes),
                    int(stderr_bytes),
                )
            )
        if int(rc) != 0:
            raise RuntimeError(
                f"Command failed [{stage}#{command_index:05}] rc={rc} cmd={command['cmd']}"
            )


    async def install(self):
        self.create_basic_logger()
        self.logger.warning(f"Start install {self.name}-{self.PackID}")

        for cmd in self.installation:
            if self.clone_shadow:
                command = self.decode_shadow_commands(cmd)
            else:
                command = dict(cmd)

            if self.logger is not None:
                self.logger.info(
                    f" Run command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> "
                )
            command_index = self._next_command_index()
            await self.run_command(command=command, stage="install", command_index=command_index)
        
        logger.remove(self.handlers['install'])
        del self.handlers['install']
        self.logger = None
        self.is_installed = True

    async def run_command(self, command, stage="execution", command_index=0):
        cmd_text = self._resolve_sample_runtime_tokens(
            command.get("cmd", ""),
            stage=stage,
            field="cmd",
        )
        cwd_text = self._resolve_sample_runtime_tokens(
            command.get("cwd", "."),
            stage=stage,
            field="cwd",
        )
        command = {
            "cmd": cmd_text,
            "cwd": self.decode_path(cwd_text),
        }
        if self.subprocess_scheduler is None:
            await self._run_command_local(command, stage=stage, command_index=command_index)
            return

        task_id = f"{self.name}-{self.PackID or 'NA'}-{stage}-{command_index:05}"
        job = SubprocessJob(
            cmd=command["cmd"],
            cwd=command["cwd"],
            shell=True,
            log_dir=None,
            log_policy="logger",
            task_id=task_id,
            meta=self._build_command_meta(stage=stage, command_index=command_index, command=command),
        )
        try:
            result = await self.subprocess_scheduler.arun(job)
        except SubprocessExecutionError:
            raise
        except Exception as exc:
            raise RuntimeError(
                f"Subprocess scheduler failed [{stage}#{command_index:05}] -> {exc}"
            ) from exc

        if self.logger is not None:
            timeout_tag = " timeout=1" if bool(result.timed_out) else ""
            self.logger.info(
                "Command done [{}#{:05}] rc={} dur={:.3f}s out={}B err={}B{}".format(
                    stage,
                    command_index,
                    result.returncode,
                    float(result.duration_sec),
                    int(result.stdout_bytes),
                    int(result.stderr_bytes),
                    timeout_tag,
                )
            )
        if not result.ok:
            raise RuntimeError(
                "Command failed [{}#{:05}] rc={} timeout={} cmd={}".format(
                    stage,
                    command_index,
                    result.returncode,
                    result.timed_out,
                    command["cmd"],
                )
            )

    def log_stream(self, stream, level, logger):
        # Using the child logger to handle the logging 
        for line in iter(stream.readline, ''):
            logger.log(level, "\t{}".format(line.strip()))

    async def log_stream_info(self, stream):
        async for line in stream:
            self.logger.bind(raw=True).info(f"\t{line.decode()}")

    async def log_stream_error(self, stream):
        async for line in stream:
            self.logger.bind(raw=True).info(f"\t{line.decode()}")

    async def initialize(self):
        for command in self.initialization:
            if self.clone_shadow:
                command = self.decode_shadow_commands(command)
            else:
                command = dict(command)

            if self.logger is not None:
                self.logger.info(
                    f" Run initialize command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> \n"
                )
            command_index = self._next_command_index()
            await self.run_command(command=command, stage="initialize", command_index=command_index)

    async def execute_commands(self):
        for command in self.execution['commands']:
            if self.clone_shadow:
                command = self.decode_shadow_commands(command)
            else:
                command = dict(command)

            if self.logger is not None:
                self.logger.info(
                    f" Run execution command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> "
                )
            command_index = self._next_command_index()
            await self.run_command(command=command, stage="execution", command_index=command_index)

    def execute(self, input_data, sample_info):
        self.sample_info = sample_info
        self.update_sample_logger(sample_info)
        self._command_counter = 0

        result = {}
        try:
            result = asyncio.run(self._execute_async(input_data=input_data))

        except Exception as e:
            # Capture and logging the error information 
            self.logger.error(f"Error during execution: {e}")
            raise

        finally:
            # Make sure the sample logger is closed
            self.close_sample_logger()

        return result

    async def _execute_async(self, input_data):
        result = {}
        await self.initialize()

        input_obs = await self.load_input(input_data=input_data)
        if isinstance(input_obs, dict):
            result.update(input_obs)

        await self.execute_commands()

        output_obs = await self.read_output()
        if isinstance(output_obs, dict):
            result.update(output_obs)
        return result



    async def read_output(self):
        from jarvishep.IOs.IOs import IOfile
        read_coroutines = [
            IOfile.load(
                ffile['name'],
                path=ffile['path'],
                file_type=ffile['type'],
                variables=ffile.get("variables", []),
                save=ffile['save'],
                logger=self.logger,
                PackID=self.PackID,
                sample_save_dir=self.sample_info['save_dir'],
                module=self.name,
                funcs=self.funcs
            ).read()
            for ffile in self.output
        ]
        observables = await asyncio.gather(*read_coroutines) 
        # print("Line 261 ->", observables)
        try:
            merged_observables = {key: val for d in observables for key, val in d.items()}
        except:
            merged_observables = {}
        # print("Line 262->", merged_observables)
        return merged_observables

    async def load_input(self, input_data):
        """
            Asynchronously loads input data into SLHA files based on the specified configuration.

            This method reads the configuration for each file from the `self.input` list, creates 
            instances for handling the files, and then concurrently writes the input data to these 
            files using their respective `write` methods. The operation is performed asynchronously 
            to improve performance when dealing with I/O operations and multiple files.

            Args:
            input_data (dict): The input data to be written into the files. This dictionary should 
                               contain the necessary information that matches the expected structure 
                               for each file type being written.

            The method uses `asyncio.gather` to concurrently execute all write operations for the 
            files defined in `self.input`. Each file is handled based on its configuration, including 
            the path, type, actions (variables), and whether it should be saved, along with other 
            metadata like `PackID`, the directory to save the file (`sample_save_dir`), and the 
            module name (`self.name`). The logger is used for logging purposes, and it's passed to 
            each file handler instance for consistent logging throughout the operation.

            After all files have been processed and the input data written, a log message is generated 
            to indicate completion of the loading process.
        """
        from jarvishep.IOs.IOs import IOfile
        write_coroutines = [
            IOfile.create(
                ffile["name"], 
                path=ffile['path'],
                file_type=ffile["type"],
                variables=ffile["actions"],
                save=ffile['save'],
                logger=self.logger, 
                PackID=self.PackID,
                sample_save_dir=self.sample_info['save_dir'],
                module=self.name,
                funcs=self.funcs
            ).write(input_data)  # Return directly to the coroutine
            for ffile in self.input
        ]
        # Perform all write operations concurrently and wait for them all to complete
        observables = await asyncio.gather(*write_coroutines)
        merged_observables = {key: val for d in observables for key, val in d.items()}
        return merged_observables
        # self.logger.warn(self.funcs)

    def decode_shadow_commands(self, cmd):
        command = {
            "cmd":  cmd['cmd'].replace("@PackID", self.PackID),
            "cwd":  cmd['cwd'].replace("@PackID", self.PackID)
        }
        return command

    def decode_shadow_path(self, path):
        path = self.decode_path(path)
        if "@PackID" in path:
            path = path.replace("@PackID", self.PackID)
        # print(path)
        return path
