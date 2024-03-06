#!/usr/bin/env python3 

import os, sys 
from base import Base
pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(pwd, "Sampling"))) 
from bridson import Bridson

class Distributor(Base):
    def __init__(self) -> None:
        super().__init__()

    def set_method(method) -> None:
        if method == "Bridson":
            return Bridson()
        
        