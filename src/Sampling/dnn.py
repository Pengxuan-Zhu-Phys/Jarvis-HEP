#!/usr/bin/env python3 
import os, sys
from abc import ABCMeta, abstractmethod
from mpmath.functions.functions import re
import numpy as np
import time 
from sampler import SamplingVirtial, BoolConversionError
import json
from sample import Sample
import concurrent.futures
import pandas as pd 
import sympy as sp 
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from dataconvert import DataConvert
from sympy.vector.functions import is_solenoidal
from sympy.codegen import futils
import matplotlib.pyplot as plt
from copy import deepcopy

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
            self.layers.append(nn.Linear(hidden_sizes[i - 1], hidden_sizes[i]))
        # Output layer
        self.output_layer = nn.Linear(hidden_sizes[-1], output_size)
        # Dropout layer
        self.dropout = nn.Dropout(dropout_prob)
        self.criterion = nn.MSELoss()
        self.optimizer = torch.optim.Adam(self.parameters(), lr=1.e-3)
    
    
    def set_device(self, device):
        self._device = device 
        
    @property
    def device(self):
        return self._device 
    
    def forward(self, x):
        """
        Forward pass through the network.
        
        Args:
            x (torch.Tensor): Input tensor of shape (batch_size, input_size).
        
        Returns:
            torch.Tensor: Output tensor of shape (batch_size, output_size).
        """
        for layer in self.layers:
            x = nn.functional.relu(layer(x))  # Apply ReLU activation
            x = self.dropout(x)  # Apply dropout
        x = self.output_layer(x)  # Output layer (no activation)
        return x

    def train_model(self, x, y, epochs=1000, learning_rate=0.001):
        """
        Train the model with multi-dimensional input X and output Y.
        """
        self.train()

        # Normalize Y
        self.y_mean = np.mean(y, axis=0)
        self.y_std = np.std(y, axis=0)
        if y.ndim == 1:
            y = y[:, np.newaxis]  # Ensure (num_samples, output_size)
        y = (y - self.y_mean) / self.y_std

        # Convert to PyTorch tensors
        y_tensor = torch.tensor(y, dtype=torch.float32).to(self.device)

        loss_sum = 0.
        # Training loop
        for i in range(epochs):
            self.optimizer.zero_grad()
            outputs = self(x)  # Forward pass
            loss = self.criterion(outputs, y_tensor)  # Compute loss
            loss.backward()  # Backpropagation
            self.optimizer.step()  # Update weights

            loss_sum += loss.cpu().detach().numpy()
        
        loss_ave = loss_sum / epochs
        return loss_ave

    def predict(self, x):
        yhat = self(x).cpu().detach().numpy() 
        return yhat * self.y_std + self.y_mean 


class DNN(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "DNN"
        self._P = None
        self._index = None
        self.tasks = []
        self.info = {}
        self._selectionexp = None
        self.future_to_sample = {}
        self.DNN = None 
        self._device = "cpu"
        self.dataset = pd.DataFrame()
        self._iter   = 0 
        self.plotter = concurrent.futures.ThreadPoolExecutor(max_workers=2)

    @property
    def device(self):
        return self._device 
    
    @property
    def loglike(self):
        return self._loglike

    def load_schema_file(self):
        self.schema = self.path['DNNSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def __iter__(self):
        self.initialize()  # ensure the _P is generated before iteration 
        return self

    def __next__(self):  
        is_selection = True 
        while is_selection: 
            temp = np.random.rand(self._dimensions)
            param = self.map_point_into_distribution(temp)
            if self._selectionexp: 
                is_selection = self.evaluate_selection(self._selectionexp, param)
        return temp, param

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
        self._dimensions = len(self.vars)
        self._vcolum = [var._name for var in self.vars]
        self._ucolum = ["u{}".format(ii + 1) for ii in range(len(self.vars))]
        self._niter = int(self.config['Sampling']['Bounds']['Niters'])
        self._hidden_layers = list(self.config['Sampling']['Bounds']['Hidden_layers'])
        self._batch_size = int(self.config['Sampling']['Bounds']['Batch_size'])
        self._ninit = int(self.config['Sampling']['Bounds']['Ninit'])
        self._nepoch = int(self.config['Sampling']['Bounds']['Nepoch'])
        self._learning_rate = float(self.config['Sampling']['Bounds']['Learning_rate'])
        self._outputs = list(self.config['Sampling']['Bounds']['Outputs'])
        self._ocolum = ["o{}".format(ii + 1) for ii in range(len(self._outputs))] 
        self._dnnocolum = ["dnno{}".format(ii + 1) for ii in range(len(self._outputs))] 
        self._dimoutputs = len(self._outputs)
        self._prop_new = self.config['Sampling']['Bounds'].get("Prop_new", 0.1)
        self._num_new   = int((1.0 - self._prop_new) * self._ninit)
        # self._maxp  = int(self.config['Sampling']['Point number'])
        if "selection" in self.config["Sampling"]:
            self._selectionexp = self.config["Sampling"]["selection"]

    def initialize(self):
        self.logger.warning("Initializing the Deep Neural Network Sampling")
        self.sampler = FeedforwardNN(
            input_size=self._dimensions,
            hidden_sizes=self._hidden_layers,
            output_size=self._dimoutputs)
        self.sampler.logger = self.logger
        self.set_torch_backend()
        self.check_evaluation()
        self.info['accept_avg'] = []
        self.info['logL']       = []
        self.info['sample_collect'] = pd.DataFrame()
        self.info['loss_avg']   = []
        self.info["DNN_Simu"] = os.path.join(self.info["sample"]['task_result_dir'], "SIMU")
        if not os.path.exists(self.info["DNN_Simu"]): 
            os.makedirs(self.info["DNN_Simu"])
        self.logger.warning("DNN is ready for training and sampling")

    def map_DNNoutput_to_Observable_and_logL(self, sample):
        o = self.sampler.predict(torch.tensor(sample.u, dtype=torch.float32).to(self.device)).flatten()
        sample.info['observables'] = {}
        sample.info['o'] = o
        for ii in range(len(self._ocolum)): 
            key = self._outputs[ii]
            sample.info['observables'][key] = o[ii]
        sample.info['observables']['LogL'] = self.loglike.calculate4dnn(sample.info['observables'])

    def run_nested(self):
        total_cores = os.cpu_count()
        self.logger.warning("Start Running DNN sampling")
        if self.check_output_in_likelihood():
            self.logger.error("Not all Output variables predicted from DNN entered in Likelihood!")
            sys.exit()
        dataset = self.create_dataset(self._ninit) 
        self.dataset = dataset
        print(self.dataset.data['LogL'], self.dataset.data['z'])
        #Detach to avoid tracking gradients 
        while self._iter < self._niter:
            loss_ave = self.sampler.train_model(self.dataset.u_tensor, self.dataset.o, self._nepoch, 1.e-3)
            self.info['loss_avg'].append(loss_ave)
            # dataset = self.creat_iteration_dataset(self._ninit)
            self._index = 0 
            data = []
            samples = {}
            while self._index <= self._num_new:
                temp, param = self.next_sample()
                sample = Sample(param)
                sample._u = temp 
                sample.set_config(deepcopy(self.info['sample']))
                samples[sample.uuid] = sample
                # self.loglike.calculate4dnn()
                self.map_DNNoutput_to_Observable_and_logL(sample)
                print(sample.info['observables'])
                data.append({
                    "uuid": sample.uuid, 
                    "u":    sample.u,
                    "v":    sample.params, 
                    "obs":  sample.info['observables'],
                    "o":    sample.info['o'],
                })
                self._index += 1 
                
            Ds = DataConvert(data, dtype="listdict", u_column=self._ucolum, v_column=self._vcolum, o_column=self._ocolum, device=self._device)
            Ds, accpt_ave = self.rejection_sampling(Ds)
            print(Ds)
            self._index = 0
            while self._index <= self._ninit:
                while len(self.tasks) < total_cores: 
                    if self._index < len(Ds):
                        try: 
                            sample = samples[Ds[self._index]]
                            future = self.factory.submit_task(sample.params, sample.info)
                            self.tasks.append(future)
                            self.future_to_sample[future] = sample
                        except StopIteration:
                            break                             
                    else: 
                        try: 
                            temp, param = self.next_sample()
                            sample = Sample(param)
                            sample._u = temp 
                            sample.set_config(deepcopy(self.info['sample']))
                            future = self.factory.submit_task(sample.params, sample.info)
                            self.tasks.append(future)
                            self.future_to_sample[future] = sample
                        except StopIteration:
                            break 
                    self._index += 1 
                done, _ = concurrent.futures.wait(self.tasks, timeout=0.001, return_when=concurrent.futures.FIRST_COMPLETED)
                self.tasks = [f for f in self.tasks if f not in done]
                for future in done:
                    try:
                        result = future.result()  # Retrieve result of completed sample
                        # Retrieve the corresponding sample instance
                        sample = self.future_to_sample.pop(future, None)
                        if sample:
                            # self.logger.info(f"Sample completed with result: {result}, Sample Info: {sample.info['observables']['LogL']}")
                            outputs = sample.evaluate_output(self._outputs)
                            data.append({
                                "uuid": sample.uuid,
                                "u": sample.u,
                                "v": sample.params,
                                "obs": sample.info['observables'],
                                "o": outputs,
                                "iter": self._iter
                            })
                    except Exception as e:
                        self.logger.error(f"Error processing sample: {e}")
                if not self.tasks:
                    break  
                    
            dataset = DataConvert(data, dtype="listdict", u_column=self._ucolum, v_column=self._vcolum, o_column=self._ocolum, device=self._device)
            self.dataset.update(dataset)
            self._iter += 1
            self.logger.warning("Finish {} ...".format(self._iter))      
        
    def rejection_sampling(self, ds):
        tds = deepcopy(ds.data)
        logM = tds['LogL'].max()

        acceptance_prob = np.exp(tds['LogL'].to_numpy() - logM)
        # print(acceptance_prob)
        acceptance_prob_avg = np.mean(acceptance_prob)
            
        # Perform rejection sampling
        tds['accum'] = np.cumsum(acceptance_prob)   
        query_array = np.random.random(self._num_new) * tds['accum'].max()
        indices = np.searchsorted(tds['accum'].to_numpy(), query_array, side="left") 
        tds = tds.iloc[indices]
        tds = tds.reset_index()  # Move index back to a column
        tds = tds.drop_duplicates(subset=["uuid"], keep="first")
        return list(tds['uuid']), acceptance_prob_avg
                        
    def finalize(self):
        self.plotter.shutdown()
        pass

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for DNN sampler")

    def evaluate_selection(self, expression, variables):
        return super().evaluate_selection(expression, variables)
  
    def check_evaluation(self):
        return super().check_evaluation()  
    
    def set_torch_backend(self): 
        if torch.backends.mps.is_available():
            self._device = "mps"
            self.sampler.to(torch.device("mps"))
            self.logger.warning("Using Apple MPS Backend for DNN ")
        elif torch.cuda.is_available():
            self._device = "cuda"
            self.sampler.to(torch.device("cuda"))
            self.logger.warning("Using CUDA® Backend for DNN")
        elif hasattr(torch, 'xpu') and torch.xpu.is_available():
            self._device = "xpu"
            self.sampler.to(torch.device("xpu"))
            self.logger.warning("Using XPU (Intel) Backend for DNN")
        elif hasattr(torch, 'xla') and torch.xla.is_available():
            self._device = "xla"
            self.sampler.to(torch.device("xla"))
            self.logger.warning("Using TPU (Google Cloud TPUs with XLA) Backend for DNN")
        elif hasattr(torch, 'hip') and torch.hip.is_available():
            self._device = "hip"
            self.sampler.to(torch.device("hip"))
            self.logger.warning("Using AMD ROCm (HIP) Backend for DNN")
        else:
            self._device = "cpu"
            self.sampler.to(torch.device("cpu"))
            self.logger.warning("No GPU acceleration Backend found, using CPU Backend for DNN")
        self.sampler.set_device(self._device)
        
    def create_dataset(self, num):
        total_cores = os.cpu_count()
        from copy import deepcopy
        self._index = 0
        data = []
        while self._index <= num: 
            while len(self.tasks) < total_cores:
                try: 
                    temp, param = self.next_sample()
                    sample = Sample(param)
                    sample._u = temp 
                    sample.set_config(deepcopy(self.info['sample']))
                    future = self.factory.submit_task(sample.params, sample.info)
                    self.tasks.append(future)
                    self.future_to_sample[future] = sample
                except StopIteration:
                    break 
            
                self._index += 1
            done, _ = concurrent.futures.wait(self.tasks, timeout=0.001, return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                try:
                    result = future.result()  # Retrieve result of completed sample
                    # Retrieve the corresponding sample instance
                    sample = self.future_to_sample.pop(future, None)
                    if sample:
                        # self.logger.info(f"Sample completed with result: {result}, Sample Info: {sample.info['observables']['LogL']}")
                        outputs = sample.evaluate_output(self._outputs)
                        data.append({
                            "uuid": sample.uuid,
                            "u": sample.u,
                            "v": sample.params,
                            "obs": sample.info['observables'],
                            "o": outputs,
                            "iter": self._iter
                        })
                except Exception as e:
                    self.logger.error(f"Error processing sample: {e}")
            self.tasks = [f for f in self.tasks if f not in done]
            if not self.tasks and not param:
                break  
        
        return DataConvert(data, dtype="listdict", u_column=self._ucolum, v_column=self._vcolum, o_column=self._ocolum, device=self._device)
    
    def check_output_in_likelihood(self):
        missing_vars = [var for var in self._outputs if var not in self.loglike.variables]
        return len(missing_vars) == 0
        
    def plot_status(self, Niter):
        savepath = os.path.join(self.info["DNN_Simu"], "sample.{}.png".format(Niter))
        fig = plt.figure(figsize=[10, 8])
        ax  = fig.add_axes([0.16, 0.16, 0.615, 0.82])
        axc = fig.add_axes([0.78, 0.16, 0.03, 0.82])
        axlogo = fig.add_axes([0.02, 0.016, 0.05, 0.04])
        
        df = self.dataset.loc[self.dataset['iter'] == Niter]
        a1 = ax.scatter(df['u1'], df['u2'], s=1, marker='o', c=df['LogL'], cmap="viridis", alpha=1.0)
        if Niter >= 1:
            df = self.dataset.loc[self.dataset['iter'] == Niter - 1]
            ax.scatter(df['u1'], df['u2'], s=1, marker='o', c="gray", alpha=0.7)
        if Niter >= 2:
            df = self.dataset.loc[self.dataset['iter'] == Niter - 2]
            ax.scatter(df['u1'], df['u2'], s=1, marker='o', c="gray", alpha=0.5)
        if Niter >= 3:
            df = self.dataset.loc[self.dataset['iter'] == Niter - 3]
            ax.scatter(df['u1'], df['u2'], s=1, marker='o', c="gray", alpha=0.35)
        if Niter >= 4:
            df = self.dataset.loc[self.dataset['iter'] <= Niter - 4]
            ax.scatter(df['u1'], df['u2'], s=1, marker='o', c="gray", alpha=0.25)
        
        plt.colorbar(a1, axc, ax)
                   
        plt.savefig(savepath, dpi=300)
        plt.close(fig)  # Close figure to free memory
        self.logger.info(f"Saved plot: {savepath}")
