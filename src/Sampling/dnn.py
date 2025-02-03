#!/usr/bin/env python3 
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
from mpmath.functions.functions import re
import numpy as np
import time 
from sampler import SamplingVirtial, BoolConversionError
import json
from sample import Sample
import concurrent.futures
import itertools
import sympy as sp 
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from sympy.physics.optics.utils import deviation
from _pytest.mark import param
from torch.testing._internal.common_methods_invocations import sample_inputs_grid_sampler_2d

# Add by Erdong Guo, modified by Pengxuan Zhu
class FeedforwardNN(nn.Module):
    """
    Feedforward Neural Network with ReLU activations and dropout layers.
    
    Args:
        input_size (int): Size of the input features.
        hidden_sizes (list): List of integers specifying the number of neurons in each hidden layer.
        output_size (int): Size of the output.
        dropout_prob (float): Dropout probability for dropout layers.
    """
    def __init__(self, input_size, hidden_sizes, output_size, dropout_prob=0.0):
        super(FeedforwardNN, self).__init__()
        self.layers = nn.ModuleList()  # To store all layers
        
        # Input layer
        self.layers.append(nn.Linear(input_size, hidden_sizes[0]))
        
        # Hidden layers
        for i in range(1, len(hidden_sizes)):
            self.layers.append(nn.Linear(hidden_sizes[i-1], hidden_sizes[i]))
        
        # Output layer
        self.output_layer = nn.Linear(hidden_sizes[-1], output_size)
        
        # Dropout layer
        self.dropout = nn.Dropout(dropout_prob)
    
    def forward(self, x):
        """
        Forward pass through the network.
        
        Args:
            x (torch.Tensor): Input tensor of shape (batch_size, input_size).
        
        Returns:
            torch.Tensor: Output tensor of shape (batch_size, output_size).
        """
        for layer in self.layers:
            x = F.relu(layer(x))  # Apply ReLU activation
            x = self.dropout(x)   # Apply dropout
        x = self.output_layer(x)  # Output layer (no activation)
        return x


class DNN(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "DNN"
        self._P     = None
        self._index = None
        self.tasks  = []
        self.info   = {}
        self._selectionexp = None
        self.DNN   = None 

    def load_schema_file(self):
        self.schema = self.path['DNNSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def __iter__(self):
        self.initialize()  # ensure the _P is generated before iteration 
        return self

    def __next__(self):# Stop iteration, if _P is not defined or _P is not np.array
        is_selection = True 
        while is_selection: 
            temp    = np.random.rand(self._dimensions)
            param   = self.map_point_into_distribution(temp)
            if self._selectionexp: 
                is_selection = self.evaluate_selection(self._selectionexp, param)
        # self._index += 1
        return param

    def next_sample(self):
        return self.__next__()

    def map_point_into_distribution(self, row) -> np.ndarray:
        result = {}
        for ii in range(len(row)):
            result[self.vars[ii].name] = self.vars[ii].map_standard_random_to_distribution(row[ii])
        return result

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self):
        self.load_variable()
        self._dimensions    = len(self.vars)
        self._niter         = int(self.config['Sampling']['Bounds']['Niters'])
        self._hidden_layers = list(self.config['Sampling']['Bounds']['Hidden_layers'])
        self._batch_size    = int(self.config['Sampling']['Bounds']['Batch_size'])
        self._ninit         = int(self.config['Sampling']['Bounds']['Ninit'])
        self._nepoch        = int(self.config['Sampling']['Bounds']['Nepoch'])
        self._learning_rate = float(self.config['Sampling']['Bounds']['Learning_rate'])
        self._outputs        = list(self.config['Sampling']['Bounds']['Outputs'])
        self._dimoutputs    = len(self._outputs)
        # self._maxp  = int(self.config['Sampling']['Point number'])
        if "selection" in self.config["Sampling"]:
            self._selectionexp = self.config["Sampling"]["selection"]

    def initialize(self):
        self.logger.warning("Initializing the Deep Neural Network Sampling")
        self.sampler = FeedforwardNN(
            input_size=self._dimensions, 
            hidden_sizes=self._hidden_layers, 
            output_size=self._dimoutputs)
        self.set_torch_backend()
        self.check_evaluation()
        self.info["DNN_Simu"] = os.path.join(self.info["sample"]['task_result_dir'], "SIMU")
        if not os.path.exists(self.info["DNN_Simu"]): 
            os.makedirs(self.info["DNN_Simu"])
        self.logger.warning("DNN is ready for training and sampling")

    def run_nested(self):
        self.logger.warning("Start Running DNN sampling")
        pass 

    def finalize(self):
        pass

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Random sampler")

    def evaluate_selection(self, expression, variables):
        return super().evaluate_selection(expression, variables)
  
    def check_evaluation(self):
        return super().check_evaluation()  
    
    def set_torch_backend(self): 
        if torch.backends.mps.is_available():
            self.sampler.to(torch.device("mps"))
            self.logger.warning("Using Apple MPS Backend for DNN ")
        elif torch.cuda.is_available():
            self.sampler.to(torch.device("cuda"))
            self.logger.warning("Using CUDA® Backend for DNN")
        elif hasattr(torch, 'xpu') and torch.xpu.is_available():
            self.sampler.to(torch.device("xpu"))
            self.logger.warning("Using XPU (Intel) Backend for DNN")
        elif hasattr(torch, 'xla') and torch.xla.is_available():
            self.sampler.to(torch.device("xla"))
            self.logger.warning("Using TPU (Google Cloud TPUs with XLA) Backend for DNN")
        elif hasattr(torch, 'hip') and torch.hip.is_available():
            self.sampler.to(torch.device("hip"))
            self.logger.warning("Using AMD ROCm (HIP) Backend for DNN")
        else:
            self.sampler.to(torch.device("cpu"))
            self.logger.warning("No GPU acceleration Backend found, using CPU Backend for DNN")
        
        
    def create_dataset(self, num):
        total_cores = os.cpu_count()
        from copy import deepcopy
        self._index = 0
        
        while self._index <= num: 
            while len(self.tasks) < total_cores:
                try: 
                    param = self.next_sample()
                    sample = Sample(param)
                    sample.set_config(deepcopy(self.info['sample']))
                    future = self.factory.submit_task(sample.params, sample.info)
                    self.tasks.append(future)
                    self.future_to_sample[future] = sample
                except StopIteration:
                    break 
            done, _ = concurrent.futures.wait(self.tasks, timeout=0.1, return_when=concurrent.futures.FIRST_COMPLETED)
            # Process completed futures
            for future in done:
                try:
                    result = future.result()  # Retrieve result of completed sample
                    # Retrieve the corresponding sample instance
                    sample = self.future_to_sample.pop(future, None)
                    if sample:
                        self.logger.info(f"Sample completed with result: {result}, Sample Info: {sample.info['observables']['LogL']}")
                        outputs = sample.evaluate_output(self._outputs)
                        self._index += 1
                except Exception as e:
                    self.logger.error(f"Error processing sample: {e}")