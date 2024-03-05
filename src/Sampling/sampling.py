#!/usr/bin/env python3 
from base import Base
from abc import ABCMeta, abstractmethod

class SamplingVirtial(Base):
    __metaclass__ = ABCMeta
    def __init__(self) -> None:
        super().__init__()
        self.schema         = None

    @abstractmethod
    def load_schema_file(self) -> None:
        pass 