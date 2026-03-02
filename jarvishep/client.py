#!/usr/bin/env python3
from __future__ import annotations

import os
import sys

from jarvishep.project_scaffold import PROJECT_SUBDIRS, create_project_scaffold


def _mkproject_fast_path(argv) -> int | None:
    if "--mkproject" not in argv:
        return None

    conflicts = [
        flag for flag in (
            "--plot",
            "--convert",
            "--monitor",
            "--install-dependencies",
            "--check-modules",
        )
        if flag in argv
    ]
    if conflicts:
        print(f"[Jarvis-HEP] --mkproject cannot be combined with: {' '.join(conflicts)}")
        return 2

    i = argv.index("--mkproject")
    if i + 1 >= len(argv) or argv[i + 1].startswith("-"):
        print("[Jarvis-HEP] Missing project name for --mkproject")
        return 2

    project_name = argv[i + 1]
    try:
        project_root = create_project_scaffold(project_name, cwd=os.getcwd())
    except ValueError as exc:
        print(f"[Jarvis-HEP] {exc}")
        return 2
    except FileExistsError as exc:
        print(f"[Jarvis-HEP] Project directory already exists: {exc}")
        return 1

    print(f"[Jarvis-HEP] Project scaffold created at: {project_root}")
    print(f"[Jarvis-HEP] Created folders: {', '.join(PROJECT_SUBDIRS)}")
    return 0


def main(argv=None) -> int:
    argv = sys.argv if argv is None else argv

    fast_path_code = _mkproject_fast_path(argv)
    if fast_path_code is not None:
        return fast_path_code

    from jarvishep.core import Core

    jc = Core()
    jc.initialization()
    if getattr(jc.args, "mkproject", None):
        jc.mkproject()
    elif jc.scan_mode:
        jc.run_sampling()
    elif jc.args.cvtDB:
        jc.convert()
    elif jc.args.plot:
        jc.plot()
    elif jc.args.monitor:
        jc.monitor()
    elif jc.args.bdREQ:
        from jarvishep.build_requirements import install_requirements
        install_requirements()
    else:
        import argparse as _ap
        _ap.ArgumentParser(prog="Jarvis").print_help()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
