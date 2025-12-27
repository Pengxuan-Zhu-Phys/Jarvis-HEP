#!/usr/bin/env python3 
# JarvisPLOT CLI (robust + JSON-config + positional YAML)
from __future__ import annotations
import argparse, os, json

# Engine + parser
from .engine import Engine
from .parsing import load_config_spec

# Ensure FigureSpec import path is valid (optional safety)
from .figures.figure import FigureSpec  # noqa: F401

__all__ = ["main"]


def _load_json_config(path: str | None) -> dict:
    if not path:
        return {}
    p = os.path.abspath(path)
    if not os.path.exists(p):
        raise FileNotFoundError(f"Config JSON not found: {p}")
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f) or {}


def main(argv=None):
    """Entry point for JarvisPLOT CLI.
    - Direct CLI:     `python -m JarvisPLOT.cli demo.yaml` (positional YAML)
    - Via Jarvis:     `./Jarvis --jplot demo.yaml` (uses --jplot)
    - Template gen:   `python -m JarvisPLOT.cli --jplot-template --style paper_1x1 --out demo.yaml`
    - JSON config:    `python -m JarvisPLOT.cli --config cfg.json`
        where cfg.json may contain keys: {"mode":"render|template", "yaml":"...", "style":"...", "out":"..."}
    """
    parser = argparse.ArgumentParser(prog="JarvisPLOT")
    parser.add_argument("yaml", nargs="?", help="YAML path to render (positional; no --jplot needed)")
    parser.add_argument("--jplot", nargs=1, help="YAML path to render (for Jarvis launcher)")
    parser.add_argument("--jplot-template", action="store_true", help="Generate a starter YAML from a style")
    parser.add_argument("--style", default="paper_1x1")
    parser.add_argument("--out", default="jplot_template.yaml")
    parser.add_argument("--config", help="Optional JSON file to configure args (mode/yaml/style/out)")
    args = parser.parse_args(argv)

    # Load JSON config if provided
    cfg_json = {}
    try:
        cfg_json = _load_json_config(args.config)
    except Exception as e:
        print(f"[JarvisPLOT] Warning: failed to load --config: {e}")

    # Decide operation mode
    # Priority order for YAML (render): --jplot → positional yaml → cfg_json['yaml']
    yaml_path = (args.jplot[0] if args.jplot else None) or args.yaml or cfg_json.get("yaml")
    mode = (
        "template" if args.jplot_template or (cfg_json.get("mode") == "template") else
        ("render" if yaml_path else None)
    )

    if mode == "render":
        if not yaml_path:
            parser.error("Missing YAML path for render mode.")
        cfg = load_config_spec(yaml_path)
        return Engine(cfg).run()

    if mode == "template":
        style = args.style or cfg_json.get("style", "paper_1x1")
        out = args.out or cfg_json.get("out", "jplot_template.yaml")
        try:
            from .templates import generator  # lazy import
            generator.write(style, out)
        except Exception:
            # Fallback: minimal template using base.yaml if available
            import yaml, pathlib
            here = pathlib.Path(__file__).resolve().parent
            base = here / "templates" / "base.yaml"
            data = {}
            if base.exists():
                with base.open("r", encoding="utf-8") as f:
                    data = yaml.safe_load(f) or {}
            data["style"] = style
            outp = pathlib.Path(out)
            outp.parent.mkdir(parents=True, exist_ok=True)
            with outp.open("w", encoding="utf-8") as f:
                yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)
        print(f"Template written to {os.path.abspath(out)}")
        return 0

    # No mode determined → show help
    parser.print_help()
    return 2

if __name__ == "__main__":
    raise SystemExit(main())