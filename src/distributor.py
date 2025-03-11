#!/usr/bin/env python3 

import os, sys 
from base import Base
pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(pwd, "Sampling"))) 

class Distributor(Base):
    def __init__(self) -> None:
        super().__init__()

    def set_method(method) -> None:
        match method:
            case "Bridson":
                from Sampling.bridson import Bridson
                return Bridson()
            case "Dynesty":
                from Sampling.dynesty import Dynesty 
                return Dynesty()
            case "Grid":
                from Sampling.grid import Grid
                return Grid()
            case "Random":
                from Sampling.randoms import RandomS
                return RandomS()
            case "DNN":
                from Sampling.dnn import DNN
                return DNN()
            case "TPMCMC":
                from Sampling.tpmcmc import TPMCMC 
                return TPMCMC()
        
        