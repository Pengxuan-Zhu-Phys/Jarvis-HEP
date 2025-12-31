# layers/base.py
class LayerBase:
    TYPE = "base"
    def __init__(self, spec): self.spec = spec
    def render(self, ax, ctx): raise NotImplementedError

# layers/__init__.py
from .base import LayerBase
from .scatter import ScatterLayer
from .contour import ContourLayer
from .line import LineLayer
from .colorbar import ColorbarLayer
from .band import BandLayer
from .ring_guides import RingGuidesLayer

registry = {
    "scatter": ScatterLayer,
    "contour": ContourLayer,
    "line": LineLayer,
    "colorbar": ColorbarLayer,
    "band": BandLayer,
    "ring_guides": RingGuidesLayer,
}