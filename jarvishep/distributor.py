#!/usr/bin/env python3 

from jarvishep.base import Base

class Distributor(Base):
    def __init__(self) -> None:
        super().__init__()

    def set_method(method) -> None:
        match method:
            case "Bridson":
                from jarvishep.Sampling.bridson import Bridson
                return Bridson()
            case "Dynesty":
                from jarvishep.Sampling.dynesty import Dynesty 
                return Dynesty()
            case "MultiNest":
                from jarvishep.Sampling.multinest import MultiNest
                return MultiNest()
            case "Grid":
                from jarvishep.Sampling.grid import Grid
                return Grid()
            case "Random":
                from jarvishep.Sampling.randoms import RandomS
                return RandomS()
            case "DNN":
                from jarvishep.Sampling.dnn import DNN
                return DNN()
            case "TPMCMC":
                from jarvishep.Sampling.tpmcmc import TPMCMC 
                return TPMCMC()
            case "MCMC":
                from jarvishep.Sampling.mcmc import MCMC 
                return MCMC()
            case "Diver":
                from jarvishep.Sampling.diver import Diver
                return Diver()
            case _:
                supported = ["Bridson", "Dynesty", "MultiNest", "Grid", "Random", "DNN", "TPMCMC", "MCMC", "Diver"]
                raise ValueError(f"Unknown Sampling.Method={method!r}. Supported: {', '.join(supported)}")
        
        
