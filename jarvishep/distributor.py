#!/usr/bin/env python3 

from jarvishep.base import Base

class Distributor(Base):
    RESUME_SUPPORT_STATUS = {
        "Bridson": "implemented",
        "Dynesty": "implemented",
        "MultiNest": "implemented",
        "Grid": "implemented",
        "Random": "implemented",
        "CSV": "implemented",
        "DNN": "implemented",
        "PTMCMC": "implemented",
        "RLTPMCMC": "implemented",
        "MCMC": "implemented",
        "AMMCMC": "implemented",
        "RobustAM": "implemented",
        "DRAM": "implemented",
        "DEMCMC": "implemented",
        "DREAM": "implemented",
        "DREAMLite": "implemented",
        "EnsembleMCMC": "implemented",
        "PTEnsemble": "implemented",
        "SliceMCMC": "implemented",
        "ESS": "implemented",
        "MALA": "implemented",
        "HMC": "implemented",
        "NUTS": "implemented",
        "ToyMCMC": "implemented",
        "Diver": "implemented",
    }
    _VALID_RESUME_STATUSES = {"implemented", "intentionally unsupported"}

    def __init__(self) -> None:
        super().__init__()

    @classmethod
    def get_resume_status(cls, method: str) -> str:
        return cls.RESUME_SUPPORT_STATUS.get(method, "intentionally unsupported")

    @classmethod
    def list_resume_statuses(cls) -> dict:
        return dict(cls.RESUME_SUPPORT_STATUS)

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
            case "CSV":
                from jarvishep.Sampling.csv_sampler import CSVSampler
                return CSVSampler()
            case "DNN":
                from jarvishep.Sampling.dnn import DNN
                return DNN()
            case "PTMCMC":
                from jarvishep.Sampling.tpmcmc import PTMCMC 
                return PTMCMC()
            case "RLTPMCMC":
                from jarvishep.Sampling.rltpmcmc import RLTPMCMC
                return RLTPMCMC()
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
            case "DREAM":
                from jarvishep.Sampling.dream import DREAM
                return DREAM()
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
            case "ESS":
                from jarvishep.Sampling.ess import ESS
                return ESS()
            case "MALA":
                from jarvishep.Sampling.mala import MALA
                return MALA()
            case "HMC":
                from jarvishep.Sampling.hmc import HMC
                return HMC()
            case "NUTS":
                from jarvishep.Sampling.nuts import NUTS
                return NUTS()
            case "ToyMCMC":
                from jarvishep.Sampling.mcmc_standard import MCMC
                return MCMC(method_name="ToyMCMC")
            case "Diver":
                from jarvishep.Sampling.diver import Diver
                return Diver()
            case _:
                supported = sorted(Distributor.RESUME_SUPPORT_STATUS.keys())
                raise ValueError(f"Unknown Sampling.Method={method!r}. Supported: {', '.join(supported)}")
        
        
