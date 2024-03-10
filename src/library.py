#!/usr/bin/env python3 

from copy import deepcopy
import json
import logging
import os, sys
from base import Base

jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

class Library(Base):
    def __init__(self) -> None:
        super().__init__()

    def set_logger(self, logger) -> None:
        self.logger = logger 
        self.logger.warning("Ready for preparing the Supporting Libraries ...")
        

    
