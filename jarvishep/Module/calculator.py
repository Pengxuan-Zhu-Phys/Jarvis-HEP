#!/usr/bin/env python3 
import asyncio
import os
import signal
import time

import sympy as sp
from loguru import logger

from jarvishep.Module.module import Module
from jarvishep.async_subprocess import SubprocessExecutionError, SubprocessJob


class CalculatorModule(Module):
    def __init__(self, name, config):
        super().__init__(name, selection=config.get("selection"))
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
            self.timeout            = self._normalize_timeout(config.get("timeout"))
            self.input              = config["execution"].get("input", [])
            self.output             = config["execution"].get("output", [])
            self.basepath           = config['path']
            self.handlers           = {}
            self.modes              = False
            self.is_installed       = False
            self.installation_event = None
            self.is_busy            = False
            self.PackID             = None
            self.sample_info        = {}
            self._funcs             = {}
            self.subprocess_scheduler = None
            self.io_manager         = None
            self.run_summary_collector = None
            self._command_counter   = 0
            self.analyze_config()

    @staticmethod
    def _normalize_timeout(value):
        if value is None:
            return None
        try:
            timeout = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"Unsupported Calculator timeout '{value}'. Expected a positive number or null.")
        if timeout <= 0:
            return None
        return timeout

    def assign_ID(self, PackID):
        self.PackID = PackID 

    def analyze_config(self):
        # get the variables for inputs and outputs, prepare for the workflow chart 
        from jarvishep.inner_func import build_expression_context, update_const, update_funcs
        parse_locals, _ = build_expression_context(
            funcs=update_funcs({}),
            consts=update_const({}),
        )

        def add_dump_variable(ipf, var):
            if "expression" in var.keys():
                expr = sp.sympify(var["expression"], locals=parse_locals)
                varis = set(expr.free_symbols)
                for vv in varis:
                    self.inputs[str(vv)] = None
                var.update({"inc": varis})
                ipf["variables"][var["name"]] = var
            else:
                self.inputs[var["name"]] = None
                ipf["variables"][var["name"]] = var

        for ipf in self.input:
            if ipf['type'] == "SLHA":
                ipf['variables'] = {}
                for act in ipf['actions']:
                    if act['type'] == "Replace":
                        for var in act['variables']:
                            add_dump_variable(ipf, var)
                    elif act['type'] == "SLHA":
                        for var in act['variables']:
                            add_dump_variable(ipf, var)
                    elif act['type'] == "File":
                        self.inputs[act['source']] = None
                        ipf['variables'][act['source']] = {"name": act['source']}
                # print(ipf['variables'])

            elif "actions" in ipf:
                ipf['variables'] = {}
                for act in ipf['actions']:
                    if act['type'] == "Dump":
                        for var in act['variables']:
                            add_dump_variable(ipf, var)
                # pprint(ipf['variables'])

            else:
                for ipv in ipf.get('variables', []):
                    if isinstance(ipv, dict) and ipv.get('name'):
                        self.inputs[ipv['name']] = None
        # pprint(self.inputs)
        for opf in self.output:
            # self.outputs[opf['name']] = None
            if opf['type'] == "File":
                pass 
            else:
                for opv in opf.get('variables', []):
                    if isinstance(opv, dict) and opv.get('name'):
                        self.outputs[opv['name']] = None

    @property
    def funcs(self):
        return self._funcs

    def set_funcs(self, funcs):
        super().set_funcs(funcs)

    def set_subprocess_scheduler(self, scheduler):
        self.subprocess_scheduler = scheduler

    def set_io_manager(self, io_manager):
        self.io_manager = io_manager

    def set_run_summary_collector(self, collector):
        self.run_summary_collector = collector

    def analyze_config_multi(self):
        pass 

    def custom_format(record):
        module = record["extra"].get("module", "No module")
        if "raw" in record["extra"]:
            return "{message}\n"
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
        """Detach the sample-local logger after task completion."""
        self.logger = None 

    def _next_command_index(self) -> int:
        self._command_counter += 1
        return self._command_counter

    def _resolve_sample_runtime_tokens(self, text: str, *, stage: str, field: str) -> str:
        """Resolve runtime-only command tokens.

        `@SampleID` and `@Sdir` are resolved from `sample_info` only during task-run stages.
        Install stage intentionally skips this replacement.
        """
        if text is None:
            return ""
        raw = str(text)
        if stage == "install" or ("@SampleID" not in raw and "@Sdir" not in raw):
            return raw

        if not isinstance(self.sample_info, dict):
            raise RuntimeError(
                f"Runtime token requires sample_info during runtime stage '{stage}' for field '{field}'"
            )

        resolved = raw
        if "@SampleID" in resolved:
            sample_uuid = self.sample_info.get("uuid")
            if sample_uuid is None:
                raise RuntimeError(
                    f"@SampleID requires sample_info['uuid'] during runtime stage '{stage}' for field '{field}'"
                )
            resolved = resolved.replace("@SampleID", str(sample_uuid))
        if "@Sdir" in resolved:
            sample_save_dir = self.sample_info.get("save_dir")
            if sample_save_dir is None:
                raise RuntimeError(
                    f"@Sdir requires sample_info['save_dir'] during runtime stage '{stage}' for field '{field}'"
                )
            resolved = resolved.replace("@Sdir", str(sample_save_dir))
        return resolved

    @staticmethod
    async def _terminate_process_with_fallback(process: asyncio.subprocess.Process, grace_sec: float = 5.0) -> None:
        pid = int(getattr(process, "pid", 0) or 0)
        try:
            if pid > 0:
                os.killpg(pid, signal.SIGTERM)
            else:
                process.terminate()
        except Exception:
            try:
                process.terminate()
            except Exception:
                pass

        try:
            await asyncio.wait_for(process.wait(), timeout=max(0.1, float(grace_sec)))
            return
        except Exception:
            pass

        try:
            if pid > 0:
                os.killpg(pid, signal.SIGKILL)
            else:
                process.kill()
        except Exception:
            try:
                process.kill()
            except Exception:
                pass

    async def _run_command_local(self, command: dict, stage: str, command_index: int, timeout_sec=None):
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
                            self.logger.bind(raw=True).info(line)
            if self.logger is not None:
                tail = text_buffer.rstrip("\r")
                if tail:
                    self.logger.bind(raw=True).info(tail)
            return total

        out_task = asyncio.create_task(_drain(process.stdout, "stdout"))
        err_task = asyncio.create_task(_drain(process.stderr, "stderr"))
        timed_out = False
        try:
            if timeout_sec is None:
                rc = await process.wait()
            else:
                rc = await asyncio.wait_for(process.wait(), timeout=float(timeout_sec))
        except asyncio.TimeoutError:
            timed_out = True
            await self._terminate_process_with_fallback(process)
            rc = await process.wait()
        stdout_bytes, stderr_bytes = await asyncio.gather(out_task, err_task)
        if self.logger is not None:
            timeout_tag = " timeout=1" if timed_out else ""
            done_message = "Command done [{}#{:05}] rc={} out={}B err={}B".format(
                stage,
                command_index,
                rc,
                int(stdout_bytes),
                int(stderr_bytes),
            ) + timeout_tag
            if stage == "install":
                self.logger.bind(raw=True).info(done_message)
            else:
                self.logger.info(done_message)
        if int(rc) != 0 or timed_out:
            raise RuntimeError(
                f"Command failed [{stage}#{command_index:05}] rc={rc} timeout={timed_out} cmd={command['cmd']}"
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

    async def run_command(self, command, stage="execution", command_index=0, timeout_sec=None):
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
        collector = self.run_summary_collector
        if collector is None and isinstance(self.sample_info, dict):
            collector = self.sample_info.get("run_summary_collector")

        started_monotonic = time.monotonic()
        duration_sec = None
        timed_out = False
        ok = False
        try:
            # Installation commands are written to the per-instance installation log.
            # The shared subprocess scheduler logs to its own domain, so keeping install
            # stage on the local drain path preserves the expected log file routing.
            if stage == "install" or self.subprocess_scheduler is None:
                await self._run_command_local(
                    command,
                    stage=stage,
                    command_index=command_index,
                    timeout_sec=timeout_sec,
                )
                duration_sec = max(0.0, time.monotonic() - started_monotonic)
                ok = True
                return

            task_id = f"{self.name}-{self.PackID or 'NA'}-{stage}-{command_index:05}"
            job = SubprocessJob(
                cmd=command["cmd"],
                cwd=command["cwd"],
                shell=True,
                log_dir=None,
                log_policy="logger",
                task_id=task_id,
                timeout_sec=timeout_sec,
                stream_logger=self.logger,
                meta={
                    "module": self.name,
                    "pack_id": self.PackID,
                    "stage": stage,
                    "command_index": int(command_index),
                    "sample_uuid": self.sample_info.get("uuid") if isinstance(self.sample_info, dict) else None,
                    "cwd": command.get("cwd"),
                },
            )
            try:
                result = await self.subprocess_scheduler.arun(job)
            except SubprocessExecutionError:
                duration_sec = max(0.0, time.monotonic() - started_monotonic)
                raise
            except Exception as exc:
                duration_sec = max(0.0, time.monotonic() - started_monotonic)
                raise RuntimeError(
                    f"Subprocess scheduler failed [{stage}#{command_index:05}] -> {exc}"
                ) from exc

            duration_sec = float(result.duration_sec)
            timed_out = bool(result.timed_out)

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
            ok = True
        finally:
            if collector is not None:
                if duration_sec is None:
                    duration_sec = max(0.0, time.monotonic() - started_monotonic)
                collector.record_external_command(
                    duration_sec=duration_sec,
                    ok=ok,
                    timed_out=timed_out,
                )

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

    def _execution_timeout_error(self):
        return RuntimeError(
            f"Calculator execution timed out after {self.timeout:g}s "
            f"for {self.name}-{self.PackID or 'NA'}"
        )

    def _remaining_execution_timeout(self, deadline):
        if deadline is None:
            return None
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise self._execution_timeout_error()
        return remaining

    async def _await_with_execution_timeout(self, awaitable_factory, deadline):
        timeout_sec = self._remaining_execution_timeout(deadline)
        if timeout_sec is None:
            return await awaitable_factory()
        try:
            return await asyncio.wait_for(awaitable_factory(), timeout=timeout_sec)
        except asyncio.TimeoutError as exc:
            raise self._execution_timeout_error() from exc

    async def execute_commands(self, deadline=None):
        if deadline is None and self.timeout is not None:
            deadline = time.monotonic() + float(self.timeout)

        for command in self.execution['commands']:
            if self.clone_shadow:
                command = self.decode_shadow_commands(command)
            else:
                command = dict(command)

            command_timeout = self._remaining_execution_timeout(deadline)

            if self.logger is not None:
                self.logger.info(
                    f" Run execution command -> \n\t{command['cmd']} \n in path -> \n\t{command['cwd']} \n Screen output -> "
                )
            command_index = self._next_command_index()
            await self.run_command(
                command=command,
                stage="execution",
                command_index=command_index,
                timeout_sec=command_timeout,
            )

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

        execution_deadline = None
        if self.timeout is not None:
            execution_deadline = time.monotonic() + float(self.timeout)

        input_obs = await self._await_with_execution_timeout(
            lambda: self.load_input(input_data=input_data),
            execution_deadline,
        )
        if isinstance(input_obs, dict):
            result.update(input_obs)

        await self.execute_commands(deadline=execution_deadline)

        output_obs = await self._await_with_execution_timeout(
            self.read_output,
            execution_deadline,
        )
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
                sample_uuid=self.sample_info.get("uuid"),
                sample_save_dir=self.sample_info['save_dir'],
                module=self.name,
                funcs=self.funcs,
                io_manager=self.io_manager,
                header=ffile.get("header"),
                columns=ffile.get("columns"),
                comment=ffile.get("comment"),
                spec=ffile,
            ).read()
            for ffile in self.output
        ]
        observables = await asyncio.gather(*read_coroutines)
        merged_observables = {}
        for batch in observables:
            if isinstance(batch, dict):
                merged_observables.update(batch)
        return merged_observables

    async def load_input(self, input_data):
        """Write all configured input files concurrently within the input stage."""
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
                sample_uuid=self.sample_info.get("uuid"),
                sample_save_dir=self.sample_info['save_dir'],
                module=self.name,
                funcs=self.funcs,
                io_manager=self.io_manager,
                header=ffile.get("header"),
                columns=ffile.get("columns"),
                comment=ffile.get("comment"),
                spec=ffile,
            ).write(input_data)
            for ffile in self.input
        ]
        observables = await asyncio.gather(*write_coroutines)
        merged_observables = {}
        for batch in observables:
            if isinstance(batch, dict):
                merged_observables.update(batch)
        return merged_observables

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
