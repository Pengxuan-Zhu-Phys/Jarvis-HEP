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
                from jarvishep.Sampling.mcmc_standard import MCMC
                return MCMC()
            case "AMMCMC":
                from jarvishep.Sampling.ammcmc import AMMCMC
                return AMMCMC()
            case "RobustAM":
                from jarvishep.Sampling.robustam import RobustAM
                return RobustAM()
            case "DRAM":
                from jarvishep.Sampling.dram import DRAM
                return DRAM()
            case "DEMCMC":
                from jarvishep.Sampling.demcmc import DEMCMC
                return DEMCMC()
            case "DREAMLite" | "DREAM-lite":
                from jarvishep.Sampling.dream_lite import DREAMLite
                return DREAMLite()
            case "EnsembleMCMC":
                from jarvishep.Sampling.ensemblemcmc import EnsembleMCMC
                return EnsembleMCMC()
            case "PTEnsemble":
                from jarvishep.Sampling.pt_ensemble import PTEnsemble
                return PTEnsemble()
            case "SliceMCMC":
                from jarvishep.Sampling.slicemcmc import SliceMCMC
                return SliceMCMC()
            case "ToyMCMC":
                from jarvishep.Sampling.mcmc_standard import MCMC
                return MCMC(method_name="ToyMCMC")
            case "Diver":
                from jarvishep.Sampling.diver import Diver
                return Diver()
            case _:
                supported = ["Bridson", "Dynesty", "MultiNest", "Grid", "Random", "DNN", "TPMCMC", "MCMC", "AMMCMC", "RobustAM", "DRAM", "DEMCMC", "DREAMLite", "EnsembleMCMC", "PTEnsemble", "SliceMCMC", "ToyMCMC", "Diver"]
                raise ValueError(f"Unknown Sampling.Method={method!r}. Supported: {', '.join(supported)}")
        
        
