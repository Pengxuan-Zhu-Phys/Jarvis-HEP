#!/usr/bin/env python3

from functools import partial
import os
from pathlib import Path
import asyncio

from jarvishep.base import Base
from jarvishep.IOs import registry as io_registry

_PROJECT_MARKERS = (".jarvis-project.json", "jarvis.project.yaml")


def _detect_project_root_from_cwd() -> str | None:
    try:
        path = Path.cwd().resolve()
    except Exception:
        return None
    for candidate in [path, *path.parents]:
        for marker in _PROJECT_MARKERS:
            if (candidate / marker).exists():
                return str(candidate)
    return None


class IOfile(Base):
    def __init__(
        self,
        name,
        path,
        file_type,
        variables,
        save,
        logger,
        PackID,
        sample_uuid,
        sample_save_dir,
        module,
        funcs,
        io_manager=None,
        header=None,
        columns=None,
        comment=None,
        spec=None,
    ):
        self.logger     = logger
        self.PackID     = PackID
        self.sample_uuid = sample_uuid
        self.name       = name
        self.path_template = None
        self.file_type = file_type
        self.variables = variables
        self.save = save
        self.path = path  
        self.sample_save_dir = sample_save_dir
        self.module = module
        self.spec = dict(spec or {})

        funcs_value = funcs
        io_manager_value = io_manager
        if (
            io_manager_value is None
            and funcs_value is not None
            and not isinstance(funcs_value, dict)
            and hasattr(funcs_value, "run_blocking")
        ):
            io_manager_value = funcs_value
            funcs_value = {}

        self.funcs = funcs_value
        self.io_manager = io_manager_value

    def sync_make_dirs(self, path, *, exist_ok=True):
        os.makedirs(str(path), exist_ok=exist_ok)

    def sync_exists(self, path):
        return bool(os.path.exists(str(path)))

    def sync_read_text(self, path, *, encoding="utf-8"):
        return Path(path).read_text(encoding=encoding)

    def sync_write_text(self, path, content, *, encoding="utf-8", ensure_parent=True):
        p = Path(path)
        if ensure_parent:
            self.sync_make_dirs(p.parent, exist_ok=True)
        return int(p.write_text(content, encoding=encoding))

    async def io_run_blocking(self, fn, *args, **kwargs):
        if self.io_manager is not None:
            return await self.io_manager.run_blocking(fn, *args, **kwargs)
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, partial(fn, *args, **kwargs))

    def _resolve_runtime_tokens(self, text):
        if text is None or not isinstance(text, str):
            return text
        resolved = str(text)
        if "@PackID" in resolved:
            if self.PackID is None:
                raise ValueError(f"Unable to resolve path: {text}")
            resolved = resolved.replace("@PackID", str(self.PackID))
        if "@SampleID" in resolved:
            if self.sample_uuid is None:
                raise ValueError(f"Unable to resolve path: {text}")
            resolved = resolved.replace("@SampleID", str(self.sample_uuid))
        if "@Sdir" in resolved:
            if self.sample_save_dir is None:
                raise ValueError(f"Unable to resolve path: {text}")
            resolved = resolved.replace("@Sdir", str(self.sample_save_dir))
        return resolved

    def decode_path(self, path) -> None:
        """
        Resolves special markers in the provided path.

        Parameters:
        - path: The path to be resolved.
        - jarvis_root: The root directory path of Jarvis.

        Returns:
        - The resolved full path.
        """
        if path is None or not isinstance(path, str):
            return path

        normalized = path.replace("\\", "/")
        if normalized.startswith("&J/") and "/src/card/" in normalized:
            raise ValueError(
                "Legacy card path prefix is no longer supported: "
                f"{path}. Use project-local packaged copies such as '&J/deps/...'."
            )

        runtime_root = None
        if isinstance(getattr(self, "path", None), dict):
            runtime_root = self.path.get("task_root") or self.path.get("jpath")
        if not runtime_root:
            runtime_root = (
                os.getenv("JARVIS_HEP_TASK_ROOT")
                or os.getenv("JHEP_TASK_ROOT")
                or _detect_project_root_from_cwd()
                or os.getcwd()
            )

        if "&J" in path:
            path = path.replace("&J", runtime_root)

        # Replace the user home directory marker ~
        if "~" in path:
            path = os.path.expanduser(path)

        path = self._resolve_runtime_tokens(path)

        if path.startswith("&"):
            raise ValueError(f"Unable to resolve path: {path}")

        if "://" not in path and not os.path.isabs(path):
            path = os.path.abspath(os.path.join(runtime_root, path))

        if self.PackID is not None:
            path = path.replace("@PackID", self.PackID)
        
        return path  

    @classmethod
    def create(
        cls,
        name,
        path,
        file_type,
        variables,
        save,
        logger,
        PackID,
        sample_uuid,
        sample_save_dir,
        module,
        funcs,
        io_manager=None,
        header=None,
        columns=None,
        comment=None,
        spec=None,
    ):
        adapter = io_registry.get(file_type, "input")
        log_label = getattr(adapter, "log_label", adapter.format_name)
        logger.debug(f"Adding the file {name} as '{log_label}' type")
        return adapter.create_handler(
            name=name,
            path=path,
            file_type=file_type,
            variables=variables,
            save=save,
            logger=logger,
            PackID=PackID,
            sample_uuid=sample_uuid,
            sample_save_dir=sample_save_dir,
            module=module,
            funcs=funcs,
            io_manager=io_manager,
            **_format_options(header=header, columns=columns, comment=comment),
            spec=spec,
        )

    @classmethod
    def load(
        cls,
        name,
        path,
        file_type,
        variables,
        save,
        logger,
        PackID,
        sample_uuid,
        sample_save_dir,
        module,
        funcs,
        io_manager=None,
        header=None,
        columns=None,
        comment=None,
        spec=None,
    ):
        adapter = io_registry.get(file_type, "output")
        log_label = getattr(adapter, "log_label", adapter.format_name)
        logger.debug(f"Loading the file {name} as '{log_label}' type")
        return adapter.create_handler(
            name=name,
            path=path,
            file_type=file_type,
            variables=variables,
            save=save,
            logger=logger,
            PackID=PackID,
            sample_uuid=sample_uuid,
            sample_save_dir=sample_save_dir,
            module=module,
            funcs=funcs,
            io_manager=io_manager,
            **_format_options(header=header, columns=columns, comment=comment),
            spec=spec,
        )

class InputFile(IOfile):
    async def write(self, param_values):
        raise NotImplementedError("This method should be implemented by subclasses.")



class OutputFile(IOfile):
    async def read(self):
        """Read the varaible values in the output files """
        raise NotImplementedError("This method should be implemented by subclasses.")


def _format_options(**values):
    return {key: value for key, value in values.items() if value is not None}




class Parameter:
    def __init__(self, name, description, distribution, ptype=None):
        self.name = name
        self.description    = description
        self.distribution   = distribution
        self.type           = "Param"
        if ptype is not None: 
            self.type           = ptype
             

    def generate_value(self):
        # 根据self.distribution的类型和参数生成参数值
        # 示例代码，具体逻辑需根据分布类型实现
        if self.distribution['type'] == 'Flat':
            return (self.distribution['parameters']['min'] + self.distribution['parameters']['max']) / 2
        elif self.distribution['type'] == 'Log':
            # Log分布的值生成逻辑，这里仅为示例
            return (self.distribution['parameters']['min'] + self.distribution['parameters']['max']) / 2
        else:
            raise ValueError("Unsupported distribution type")
