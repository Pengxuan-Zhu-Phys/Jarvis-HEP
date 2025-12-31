#!/usr/bin/env python3 
from __future__ import annotations

from base import Base
from abc import ABCMeta, abstractmethod
from variables import Variable
import sympy as sp 
from inner_func import update_const, update_funcs
import numpy as np 
import sys, os 

class BoolConversionError(Exception):
    """Nothing but raise bool error"""
    pass
class SamplingVirtial(Base):
    __metaclass__ = ABCMeta
    def __init__(self) -> None:
        super().__init__()
        self.schema                     = None
        self.info: Dict[str, Any]       = {}
        self.vars: Tuple[Variable]      = None
        self.method                     = None
        self.max_workers                = 4
        self._selectionexp              = None
        self._loglike                   = None 
        self.nuisance_sampler           = None 
        self._with_nuisance             = False 
        self.bucket_alloc               = None 
        self.total_core                 = os.cpu_count()
        
    def load_nuisance_sampler(self): 
        if self.config['Sampling']['Nuisance']: 
            nui_config = self.config['Sampling']['Nuisance']
            if nui_config['Method'] == "Profile1D": 
                from Sampling.Source.Nuisance.profile1d import Profile1D
                self.nuisance_sampler = Profile1D()
                self.nuisance_sampler.set_logger(self.logger)
                self.nuisance_sampler.set_config(nui_config)
                self.info['sample']['nuisance'] = self.nuisance_sampler.get_info_card()
                self._with_nuisance = True

    @property
    def loglike(self):
        return self._loglike

    def set_max_workers(self, nn):
        self.max_workers = nn 

    @abstractmethod
    def load_schema_file(self) -> None:
        pass 

    @abstractmethod
    def set_config(self, config_info) -> None:
        pass 

    @abstractmethod
    def set_bucket_alloc(self) -> None: 
        limit = 200 
        width = 6 
        ba_config = self.config.get("Directory_Setting", None)
        if ba_config is not None: 
            limit = ba_config.get("limit", 200)
            width = ba_config.get("width", 6)

        from bucketallocator import BucketAllocator
        self.bucket_alloc = BucketAllocator(
            base_path=self.info['sample']['sample_dirs'],
            limit=limit,
            width=width,
            start_bucket=1,
        )
        self.bucket_alloc.check_and_update()
        
        
    @abstractmethod
    def init_generator(self) -> None:
        pass 
         

    @abstractmethod
    def set_factory(self) -> None: 
        pass 

    @abstractmethod
    def load_variable(self) -> None:
        variables = []
        for var_config in self.config['Sampling'].get("Variables", []):  # 使用get以防'Variables'键不存在
            var = Variable(
                name=var_config['name'],
                description=var_config['description'],
                distribution=var_config['distribution']['type'],
                parameters=var_config['distribution']['parameters']
            )
            variables.append(var)
        self.vars = tuple(variables)

    @abstractmethod
    def set_logger(self, logger) -> None:
        self.logger = logger 

    @abstractmethod 
    def combine_data(self, df_full) -> None:
        pass 
        

    
    @abstractmethod
    def set_likelihood(self, loglike):
        self._loglike = loglike
        self._loglike.logger = self.logger 
    
    @abstractmethod
    def logo_at_pos(self, pos, ax):
        from PIL import Image
        image_path = "src/icons/JarvisHEP.png"
        with Image.open(image_path) as image:
            image = np.array(image.convert("RGBA"))
            ax.axes("off")
            ax.imshow(image, extent=[pos[0]-0.25, pos[0]+0.25, pos[1]-0.25, pos[1]+0.25], zorder=100)
            ax.text(pos[0]+0.27, pos[1]+0.15, "Jarvis-HEP", ha="left", va='top', color="#0F66C3", fontfamily="sans-serif", fontsize="small", fontstyle="normal", fontweight="bold")
