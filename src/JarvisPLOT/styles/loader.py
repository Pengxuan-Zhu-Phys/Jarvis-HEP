from __future__ import annotations
import os, json
from typing import Any, Dict
from ..utils import rc_update_safely

class StyleBundle:
    """Loads a style bundle (registry -> style.json) and exposes helpers
    for figure defaults, axes presets, per-layer overrides, and theme application.
    """
    def __init__(self, root: str, manifest: Dict[str, Any]):
        self.root = root
        self.manifest = manifest or {}
        self.figure_defaults = self.manifest.get("figure_defaults", {})
        self.axes_presets = self.manifest.get("axes_presets", {})

    # ---------- loading ----------
    @classmethod
    def load(cls, name: str) -> "StyleBundle":
        here = os.path.abspath(os.path.dirname(__file__))
        reg_path = os.path.join(here, "registry.json")
        with open(reg_path, "r", encoding="utf-8") as f:
            reg = json.load(f)
        if name not in reg:
            raise KeyError(f"Unknown style bundle: {name}")
        manifest_path = os.path.join(here, reg[name])
        root = os.path.dirname(manifest_path)
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)
        return cls(root, manifest)

    # ---------- theme ----------
    def apply_theme(self) -> None:
        rel = self.manifest.get("theme")
        if not rel:
            return
        path = os.path.join(self.root, rel)
        if not os.path.exists(path):
            return
        with open(path, "r", encoding="utf-8") as f:
            theme = json.load(f)
        rc_update_safely(theme)

    # ---------- figure defaults (global and per-layer) ----------
    def figure_defaults_for(self, layer_types):
        base = dict(self.figure_defaults or {})
        by_layer = self.manifest.get("figure_defaults_by_layer", {})
        for lt in layer_types:
            spec = by_layer.get(lt)
            if spec:
                base.update(spec)
        return base

    # ---------- axes presets (named and per-layer) ----------
    def axes_preset(self, ax_name: str) -> Dict[str, Any]:
        rel = (self.axes_presets or {}).get(ax_name)
        if rel is None:
            return {}
        path = os.path.join(self.root, rel)
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}

    def axes_preset_for(self, ax_name: str, layer_types):
        if ax_name == "main":
            table = self.manifest.get("axes_presets_by_layer", {})
            for lt in layer_types:
                rel = table.get(lt)
                if rel:
                    path = os.path.join(self.root, rel)
                    try:
                        with open(path, "r", encoding="utf-8") as f:
                            return json.load(f)
                    except Exception:
                        pass
        return self.axes_preset(ax_name)

    # ---------- profiles (figure+axes+layer defaults) ----------
    def profile_for(self, layer_types, explicit: str | None = None):
        """Return a resolved profile card dict or None.
        Resolution order:
          1) explicit profile name (if provided)
          2) profile_by_layer mapping — first matching layer type wins
        """
        profs = self.manifest.get("profiles", {}) or {}
        # 1) explicit
        name = None
        if explicit:
            if explicit in profs:
                name = explicit
            else:
                # allow explicit to be a direct relative json path (advanced use)
                rel = explicit if explicit.endswith('.json') else None
                if rel:
                    path = os.path.join(self.root, rel)
                    try:
                        with open(path, "r", encoding="utf-8") as f:
                            return json.load(f)
                    except Exception:
                        return None
        # 2) by-layer routing
        if name is None:
            routing = self.manifest.get("profile_by_layer", {}) or {}
            for lt in layer_types:
                cand = routing.get(lt)
                if cand and cand in profs:
                    name = cand
                    break
        if not name:
            return None
        rel = profs.get(name)
        if not rel:
            return None
        path = os.path.join(self.root, rel)
        try:
            with open(path, "r", encoding="utf-8") as f:
                prof = json.load(f)
            return prof
        except Exception:
            return None

    def profile_name_for(self, layer_types, explicit: str | None = None):
        """Resolve and return the profile *name* according to the same
        rules as profile_for(), without opening the profile file.
        Returns a string or None.
        """
        profs = self.manifest.get("profiles", {}) or {}
        # 1) explicit
        if explicit:
            if explicit in profs:
                return explicit
            # allow explicit to be a direct relative json path; in that case
            # there is no canonical name in registry, so return None
            if explicit.endswith('.json'):
                return None
        # 2) by-layer routing
        routing = self.manifest.get("profile_by_layer", {}) or {}
        for lt in layer_types:
            cand = routing.get(lt)
            if cand and cand in profs:
                return cand
        return None

    def layer_defaults_from_profile(self, profile: dict | None, layer_type: str) -> dict:
        if not profile:
            return {}
        layers = profile.get("layers", {}) or {}
        val = layers.get(layer_type)
        return val or {}

    def expand_profile_axes(self, profile: dict | None) -> list[dict]:
        """Expand profile axes entries by loading any 'preset' JSON and merging
        with 'override'. Returns a list of ready-to-use axes cfg dicts.
        """
        if not profile:
            return []
        axes = profile.get("axes", []) or []
        out = []
        for ax in axes:
            name = ax.get("name", "main")
            preset_rel = ax.get("preset")
            override = ax.get("override", {}) or {}
            base = {}
            if preset_rel:
                path = os.path.join(self.root, preset_rel)
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        base = json.load(f)
                except Exception:
                    base = {}
            # shallow merge is fine; Figure.build may still deep-merge if needed
            merged = dict(base)
            merged.update(override)
            merged.setdefault("name", name)
            out.append(merged)
        return out

    # ---------- layer presets ----------
    def layer_preset(self, layer_type: str) -> Dict[str, Any]:
        table = self.manifest.get("layer_presets", {})
        rel = table.get(layer_type)
        if not rel:
            return {}
        path = os.path.join(self.root, rel)
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}

# Backward-compatible alias (older code referenced StyleLoader)
StyleLoader = StyleBundle
