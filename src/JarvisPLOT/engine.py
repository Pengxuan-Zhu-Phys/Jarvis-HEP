from __future__ import annotations
import os
from dataclasses import dataclass
from typing import Dict, List, Any

# Use a non-interactive backend to avoid GUI blocking
try:
    import matplotlib as _mpl
    _mpl.use("Agg")
except Exception:
    pass

from .parsing import ConfigSpec
from .figures.figure import Figure
from .datasets import DataEnv
from .styles.loader import StyleBundle


@dataclass
class RunContext:
    style: StyleBundle
    data_env: DataEnv
    namespace: Dict[str, Any]  # eval ns with named computed variables
    outputs: Dict[str, List[str]]  # figure_name -> list of saved paths
    profile: Dict[str, Any] | None = None   # resolved profile card for current figure
    profile_name: str | None = None         # name of the resolved profile


class Engine:
    def __init__(self, cfg: ConfigSpec):
        self.cfg = cfg

    def run(self) -> int:
        # Load style and apply theme
        style = StyleBundle.load(self.cfg.style_name)
        style.apply_theme()

        # Prepare data environment
        data_env = DataEnv(self.cfg.datasources, self.cfg.yaml_dir)
        
        # Resolve variables safely
        namespace: Dict[str, Any]
        try:
            namespace = self.cfg.variables.resolve_all(data_env)
        except Exception:
            namespace = data_env.namespace()

        ctx = RunContext(style=style, data_env=data_env, namespace=namespace, outputs={})

        # Execute in YAML directory so relative paths (datasources/output.dir) resolve
        cwd = os.getcwd()
        try:
            os.chdir(self.cfg.yaml_dir)
            for fig_spec in self.cfg.figures or []:
                # Resolve a profile card from the style based on layer types or explicit figure profile
                layer_types = [ls.type for ls in getattr(fig_spec, 'layers', [])]
                profile_name = None
                profile = None
                prof_fn = getattr(style, 'profile_for', None)
                if callable(prof_fn):
                    explicit = getattr(fig_spec, 'profile', None)  # allow YAML to name a profile explicitly later
                    try:
                        profile = style.profile_for(layer_types, explicit)
                        if profile:
                            # If profile was selected via layer routing, record its name for diagnostics
                            if explicit:
                                profile_name = explicit
                            else:
                                # try to infer name via profile_by_layer mapping
                                pbl = style.manifest.get('profile_by_layer', {})
                                for lt in layer_types:
                                    if lt in pbl:
                                        profile_name = pbl[lt]
                                        break
                    except Exception:
                        profile = None
                        profile_name = None

                # Attach to context so layers can introspect profile if needed
                ctx.profile = profile
                ctx.profile_name = profile_name

                # Also attach to fig_spec so Figure.build can consume defaults (axes/size) if implemented
                setattr(fig_spec, '_resolved_profile', profile)
                setattr(fig_spec, '_resolved_profile_name', profile_name)

                fig = Figure.build(fig_spec, style)        # construct runtime Figure & Axes (may consume profile)
                fig.render(ctx)                            # render all layers (can read ctx.profile)
                saved = fig.save(self.cfg.output)          # save in formats
                ctx.outputs[fig_spec.name] = saved
                fig.close()
        finally:
            os.chdir(cwd)

        # Print saved outputs for visibility
        if ctx.outputs:
            print("[JarvisPLOT] Saved outputs:")
            for name, paths in ctx.outputs.items():
                for p in paths:
                    ap = p if os.path.isabs(p) else os.path.abspath(os.path.join(self.cfg.yaml_dir, p))
                    print(f"  - {name}: {ap}")
        else:
            print("[JarvisPLOT] No figures were rendered (empty figures list?)")
        return 0
    
