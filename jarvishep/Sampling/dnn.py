#!/usr/bin/env python3 
import os, sys
import numpy as np
from jarvishep.Sampling.sampler import SamplingVirtial
import json
from jarvishep.sample import Sample
import concurrent.futures
import pandas as pd 
import torch
import torch.nn as nn
from jarvishep.dataconvert import DataConvert
import matplotlib.pyplot as plt
from copy import deepcopy
# from torch.autograd import Variable
# from dnn import device


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_json_safe(v) for v in value.tolist()]
    tensor_type = getattr(torch, "Tensor", None)
    if tensor_type is not None and isinstance(value, tensor_type):
        return _json_safe(value.detach().cpu().tolist())
    return value


def _format_kv_block(title, items):
    lines = [str(title)]
    for key, value in items:
        lines.append(f"  {key:<20} -> {value}")
    return "\n".join(lines)


# Add by Erdong Guo, modified by Pengxuan Zhu
class FeedforwardNN():
    """
    Feedforward Neural Network with ReLU activations and dropout layers.
    
    Args:
        input_size (int): Size of the input features.
        hidden_sizes (list): List of integers specifying the number of neurons in each hidden layer.
        output_size (int): Size of the output.
        dropout_prob (float): Dropout probability for dropout layers.
    """
    def __init__(self,
                 input_size: int,
                 hidden_sizes: list[int],
                 output_size: int,
                 otag,
                 dropout_prob: float = 0.0,
                 learning_rate: float = 1e-3):
        super().__init__()
        self.otag = otag

        # 1) Build the network
        layers = []
        in_features = input_size
        for h in hidden_sizes:
            layers.append(nn.Linear(in_features, h))
            layers.append(nn.ReLU())
            if dropout_prob > 0:
                layers.append(nn.Dropout(dropout_prob))
            in_features = h
        layers.append(nn.Linear(in_features, output_size))
        self.model = nn.Sequential(*layers)

        # 2) Loss & optimizer
        self.criterion = nn.MSELoss()
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)

        # 3) Normalization placeholders
        self.x_mean = np.zeros((1, input_size))
        self.x_std  = np.ones((1, input_size))
        self.y_mean = np.zeros((1, output_size))
        self.y_std  = np.ones((1, output_size))

        self._device = torch.device("cpu")
        self.path = ""
    
    def set_device(self, device):
        self._device = device 
        self.model.to(self.device)
        
    @property
    def device(self):
        return self._device 
    
    def forward(self, x):
        """Forward through the internal Sequential model."""
        return self.model(x)

    def train_model(self, x, y, epochs=1000, learning_rate=0.001):
        """
        Train the model with multi-dimensional input X and output Y.
        """
        self.model.train()

        self.x_mean = np.mean(x, axis=0)
        self.x_std  = np.std(x, axis=0)
        
        self.y_mean = np.mean(y, axis=0)
        self.y_std = np.std(y, axis=0)
        
        x = (x - self.x_mean) / self.x_std
        y = (y - self.y_mean) / self.y_std

        # Convert to PyTorch tensors
        x = torch.tensor(x, dtype=torch.float32).to(self.device)
        y = torch.tensor(y[:, None], dtype=torch.float32, device=self.device)
        loss_sum = 0.
        # Training loop
        # Reduce log noise: print progress every 20% of epochs.
        log_every = max(1, epochs // 5)
        for ii in range(epochs):
            yhat = self.model(x)
            loss = self.criterion(yhat, y)
            self.optimizer.zero_grad()
            loss.backward() 
            self.optimizer.step()
            # Update weights
            if ii % log_every == 0 or ii == epochs - 1:
                self.logger.warning(
                    "DNN Regressor -> {} | epoch -> {}/{} | loss -> {:.6g}".format(
                        self.otag, ii, epochs, float(loss.item())
                    )
                )
            loss_sum += float(loss.detach().cpu().item())
        
        loss_ave = float(loss_sum / epochs)
        self.save()
        return loss_ave



    def predict(self, x):
        """
        Predict with the internal nn.Sequential model, accepting numpy arrays, pandas DataFrames, lists, or torch Tensors.
        Returns denormalized numpy output of shape (n_samples, output_size).
        """
        # Move model to correct device and set to eval mode
        self.model.to(self.device)
        self.model.eval()

        # Prepare input tensor
        if isinstance(x, np.ndarray):
            arr = x
        elif isinstance(x, torch.Tensor):
            arr = x.cpu().numpy()
        else:
            # e.g., pandas DataFrame or list
            arr = np.asarray(x)
        # Normalize
        x_norm = (arr - self.x_mean) / self.x_std
        tx = torch.tensor(x_norm, dtype=torch.float32, device=self.device)

        # Forward pass without gradient tracking
        with torch.no_grad():
            yhat = self.model(tx)

        # Convert back to numpy and denormalize
        yhat_np = yhat.cpu().numpy()
        return yhat_np * self.y_std + self.y_mean

    def save(self):
        """
        Save the FeedforwardNN model checkpoint and metadata.
        """
        # Save the full nn.Sequential model
        torch.save(self.model, f"{self.path}.pth")
        # Prepare and save JSON metadata
        info = {
            "x_mean": self.x_mean.tolist(),
            "x_std": self.x_std.tolist(),
            "y_mean": self.y_mean.tolist(),
            "y_std": self.y_std.tolist(),
            "model_state": {
                k: v.cpu().numpy().tolist()
                for k, v in self.model.state_dict().items()
            }
        }
        with open(f"{self.path}.json", "w") as f:
            json.dump(info, f)

class Classifier():
    def __init__(self, model, x_dim, learning_rate=1e-2):
        self.model = model
        self.x_dim = x_dim

        # Training method
        self.criterion = nn.BCELoss()
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)

        # Data normalization parameters
        self.x_mean = np.zeros([1, x_dim])
        self.x_std = np.ones([1, x_dim])

    def set_device(self, device):
        self._device = device 
        self.model.to(device)
        
    @property
    def device(self):
        return self._device 
    

    def train_model(self, x, y, epoch=1000):
        # Normalize the data
        self.x_mean = np.mean(x, axis=0)
        self.x_std = np.std(x, axis=0)

        x = (x - self.x_mean) / self.x_std

        # Train
        x = torch.tensor(x, dtype=torch.float32, device=self.device)
        y = torch.tensor(y, dtype=torch.float32, device=self.device).unsqueeze(1)  # ← expand

        # Reduce log noise: print progress every 20% of epochs.
        log_every = max(1, epoch // 5)
        for i in range(epoch):
            yhat = self.model(x)
            loss = self.criterion(yhat, y)
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()

            if i % log_every == 0 or i == epoch - 1:
                self.logger.warning(
                    "DNN Classifier -> epoch -> {}/{} | loss -> {:.6g}".format(
                        i, epoch, float(loss.item())
                    )
                )
        self.save()

    def predict(self, x):
        x = (x - self.x_mean) / self.x_std
        x = torch.tensor(x, dtype=torch.float32, device=self.device)
        return self.model(x).cpu().detach().numpy()
    
    def save(self):
        # Save the full model checkpoint
        torch.save(self.model, f"{self.path}.pth")
        # Prepare JSON metadata
        with open(f"{self.path}.json", 'w') as f1:
            info = {
                "x_mean": self.x_mean.tolist(),
                "x_std": self.x_std.tolist(),
                "model_state": {
                    k: v.cpu().numpy().tolist()
                    for k, v in self.model.state_dict().items()
                }
            }
            json.dump(info, f1)
            



class DNN(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "DNN"
        self._P = None
        self._index = None
        self.tasks = set()
        self.info = {}
        self._selectionexp = None
        self.future_to_sample = {}
        self.DNN = None 
        self._device = "cpu"
        self.dataset = pd.DataFrame()
        self.df      = pd.DataFrame()
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
        is_selection = False 
        while not is_selection: 
            temp = np.random.rand(self._dimensions)
            param = self.map_point_into_distribution(temp)
            if self._selectionexp: 
                is_selection = self.evaluate_selection(self._selectionexp, param)
            else: 
                is_selection = True
        return param

    def next_from_df(self):
        first_row = dict(self.df.iloc[0])
        self.df = self.df.drop(self.df.index[0])
        return first_row
        

    def next_sample(self):
        if not self.df.empty:
            return self.next_from_df()
        else:
            return self.__next__()

    def map_point_into_distribution(self, row) -> np.ndarray:
        result = {}
        for ii in range(len(row)):
            result[self.vars[ii].name] = self.vars[ii].map_standard_random_to_distribution(row[ii])
        return result

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("DNN Sampler -> initializing")

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
        self._ocolum = self._outputs
        self._dnnocolum = ["dnno{}".format(ii + 1) for ii in range(len(self._outputs))] 
        self._dimoutputs = len(self._outputs)
        self._prop_new = self.config['Sampling']['Bounds'].get("Prop_new", 0.1)
        self._num_new   = int((1.0 - self._prop_new) * self._ninit)
        # self._maxp  = int(self.config['Sampling']['Point number'])
        if "selection" in self.config["Sampling"]:
            self._selectionexp = self.config["Sampling"]["selection"]

    def initialize(self):
        self.logger.warning("DNN Sampler -> initialize network and runtime")
        # self.sampler = FeedforwardNN(
        #     input_size=self._dimensions,
        #     hidden_sizes=self._hidden_layers,
        #     output_size=self._dimoutputs)
        self.set_torch_backend()        
        self.regressors = [
            FeedforwardNN(
                input_size = self._dimensions,
                hidden_sizes = self._hidden_layers,
                output_size = 1,
                otag = output_name,
                dropout_prob = 0.1,
                learning_rate = self._learning_rate
            )
                for output_name in self._outputs
            ]
        for reg in self.regressors:
            reg.logger = self.logger
            reg.set_device(self.device)
        self.check_evaluation()
        self.classifier = Classifier(
            nn.Sequential(nn.Linear(self._dimensions, 100), nn.ReLU(),
                          nn.Linear(100, 100), nn.ReLU(),
                          nn.Linear(100, 100), nn.ReLU(), 
                          nn.Linear(100, 100), nn.ReLU(),
                          nn.Linear(100, 1), nn.Sigmoid()
                          ), 
            self._dimensions,
            learning_rate = self._learning_rate
            )
        self.classifier.set_device(self.device)
        self.classifier.logger = self.logger


        self.info['accept_avg'] = []
        self.info['LogL']       = []
        self.info['sample_collect'] = pd.DataFrame()
        self.info['loss_avg']   = []
        self.info["DNN_Simu"] = os.path.join(self.info["sample"]['task_result_dir'], "SIMU")
        self.info["DNN_info"] = {
            "run":  os.path.join(self.info["DNN_Simu"], "run_info.json"),
            "iter_data":    [],
            "loss": {"model_{}".format(mod.otag): [] for mod in self.regressors},
            "model":    {
                "model_{}".format(mod.otag): os.path.join(self.info["DNN_Simu"], "Model_{}".format(mod.otag)) for mod in self.regressors
            }
        }
        for mod in self.regressors:
            mod.path = self.info["DNN_info"]['model']["model_{}".format(mod.otag)]
        self.info["DNN_info"]['model']["model_Classifier"] = os.path.join(self.info["DNN_Simu"], "Model_Classifier")
        self.classifier.path = os.path.join(self.info["DNN_Simu"], "Model_Classifier")
        if not os.path.exists(self.info["DNN_Simu"]): 
            os.makedirs(self.info["DNN_Simu"])
        self._log_dnn_runtime_parameters()
        self.logger.warning("DNN Sampler -> ready for training and sampling")

    def map_DNNoutput_to_Observable_and_logL(self, Dx):
        for reg in self.regressors:
            Dx[reg.otag] = reg.predict(Dx.to_numpy())
        Dx["LogL"] = self.loglike.calculate4dnn(Dx)
        return Dx 
        
    def GenParameters(self, num):
        idx = 0
        x = []
        while idx < num:
            para = self.__next__() 
            x.append(para)
            idx += 1
        return pd.DataFrame(x)

    def recommand(self):
        self.df = pd.DataFrame()
        while self.df.shape[0] < self._num_new: 
            x = self.GenParameters(self._ninit)
            # Filter out the non-physical and excluded points.
            y = self.classifier.predict(x.to_numpy())
            y = x[y > 0.5]
            Dfx = self.map_DNNoutput_to_Observable_and_logL(y)
            rad = np.random.rand(y.shape[0])
            ll = Dfx['LogL'] - max(Dfx['LogL'])
            cond = np.log(rad) < ll
            Dfx = Dfx[cond]
            x = x[cond]
            self.logger.info(
                "DNN Recommend -> max LogL -> {} | accepted shape -> {} | total shape -> {}".format(
                    max(Dfx['LogL']), Dfx.shape, self.df.shape
                )
            )
            self.df = pd.concat([self.df, x], ignore_index=True)

    def run_nested(self):
        self.logger.warning("DNN Sampling -> start")
        if not self.check_output_in_likelihood():
            missing = getattr(self, "_missing_likelihood_outputs", [])
            self.logger.error(
                "DNN Output Check -> failed | missing outputs in likelihood -> {}".format(missing)
            )
            sys.exit(2)
            
        # First data set 
        dataset = self.create_dataset(self._ninit) 
        dpth = os.path.join(self.info["DNN_Simu"], "iter.0.csv")
        dataset.save(dpth)
        self.info['DNN_info']['iter_data'].append(dpth)
        self.dataset = dataset
        self._log_dataset_snapshot("init", self.dataset)
        self._write_run_info()
            
        while self._iter < self._niter:
            self._iter += 1
            self.logger.warning(
                "DNN Iteration -> {}/{} | stage -> classifier training".format(
                    self._iter, self._niter
                )
            )
            self.classifier.train_model(self.dataset.v, self.dataset.valid)
            for ii, reg_model in enumerate(self.regressors):
                self.logger.warning(
                    "DNN Iteration -> {}/{} | stage -> regressor training | index -> {} | output -> {}".format(
                        self._iter, self._niter, ii, reg_model.otag
                    )
                )
                loss = reg_model.train_model(self.dataset.v, self.dataset.o(reg_model.otag))
                # record loss
                self.info["DNN_info"]["loss"][f"model_{reg_model.otag}"].append(loss)

            self.logger.warning(
                "DNN Iteration -> {}/{} | stage -> recommend points".format(
                    self._iter, self._niter
                )
            )
            self.recommand()
            
            self.logger.warning(
                "DNN Iteration -> {}/{} | stage -> exact observable evaluation".format(
                    self._iter, self._niter
                )
            )
            dataset = self.create_dataset(self._ninit)
            dpth = os.path.join(self.info["DNN_Simu"], "iter.{}.csv".format(self._iter))
            dataset.save(dpth)
            self.info['DNN_info']['iter_data'].append(dpth)
            self.dataset.update(dataset)
            self._log_dataset_snapshot(f"iter-{self._iter}", self.dataset)
            self._write_run_info()
            self.logger.warning(
                "DNN Iteration -> {}/{} | stage -> completed".format(self._iter, self._niter)
            )
        
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
        self._write_run_info()
        self.plotter.shutdown()
        pass

    def _write_run_info(self):
        with open(self.info["DNN_info"]["run"], "w") as f1:
            json.dump(_json_safe(self.info["DNN_info"]), f1)

    def _log_dnn_runtime_parameters(self):
        bounds = self.config.get("Sampling", {}).get("Bounds", {})
        items = [
            ("method", self.method),
            ("device", self.device),
            ("dimensions", self._dimensions),
            ("niter", self._niter),
            ("ninit", self._ninit),
            ("num_new", self._num_new),
            ("nepoch", self._nepoch),
            ("batch_size", self._batch_size),
            ("learning_rate", self._learning_rate),
            ("hidden_layers", self._hidden_layers),
            ("outputs", self._outputs),
            ("selection", self._selectionexp if self._selectionexp else "None"),
            ("max_workers", self.max_workers),
            ("bounds keys", sorted(list(bounds.keys()))),
        ]
        self.logger.warning(_format_kv_block("DNN Runtime Parameters", items))

    def _log_dataset_snapshot(self, tag, dataset):
        if dataset is None or getattr(dataset, "data", None) is None:
            return
        size = int(len(dataset.data))
        valid = np.asarray(dataset.valid, dtype=bool) if size > 0 else np.asarray([], dtype=bool)
        valid_ratio = float(valid.mean()) if valid.size > 0 else 0.0
        logl = np.asarray(dataset.LogL, dtype=float) if size > 0 else np.asarray([], dtype=float)
        finite = logl[np.isfinite(logl)] if logl.size > 0 else np.asarray([], dtype=float)
        if finite.size > 0:
            logl_min = float(np.min(finite))
            logl_max = float(np.max(finite))
            logl_span = "[{:.6g}, {:.6g}]".format(logl_min, logl_max)
        else:
            logl_span = "None"
        self.logger.warning(
            "DNN Dataset -> {} | size -> {} | valid ratio -> {:.4f} | finite LogL -> {}".format(
                tag, size, valid_ratio, logl_span
            )
        )

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("DNN WorkerFactory -> ready")
    
    def set_torch_backend(self): 
        if torch.backends.mps.is_available():
            self._device = "mps"
            self.logger.warning("DNN Backend -> Apple MPS")
        elif torch.cuda.is_available():
            self._device = "cuda"
            self.logger.warning("DNN Backend -> CUDA")
        elif hasattr(torch, 'xpu') and torch.xpu.is_available():
            self._device = "xpu"
            self.logger.warning("DNN Backend -> XPU (Intel)")
        elif hasattr(torch, 'xla') and torch.xla.is_available():
            self._device = "xla"
            self.logger.warning("DNN Backend -> TPU (XLA)")
        elif hasattr(torch, 'hip') and torch.hip.is_available():
            self._device = "hip"
            self.logger.warning("DNN Backend -> AMD ROCm (HIP)")
        else:
            self._device = "cpu"
            self.logger.warning("DNN Backend -> CPU")
        
    def create_dataset(self, num):
        total_cores = os.cpu_count() or 1
        from copy import deepcopy
        self._index = 0
        self.tasks = set()
        self.future_to_sample = {}
        data = []
        exhausted = False

        while (not exhausted) or self.tasks:
            while not exhausted and len(self.tasks) < total_cores and self._index < num:
                try: 
                    param = self.next_sample()
                    sample = Sample(param)
                    sample.set_config(deepcopy(self.info['sample']))
                    future = self.factory.submit_task(sample.info)
                    self.tasks.add(future)
                    self.future_to_sample[future] = sample
                except StopIteration:
                    exhausted = True
                    break

                self._index += 1
            if self._index >= num:
                exhausted = True

            if not self.tasks:
                continue

            done, _ = concurrent.futures.wait(
                self.tasks,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            for future in done:
                sample = self.future_to_sample.pop(future, None)
                try:
                    future.result()
                    if sample:
                        outputs = sample.evaluate_output(self._outputs)
                        data.append({
                            "uuid": sample.uuid,
                            "v": sample.params,
                            "obs": sample.info['observables'],
                            "o": outputs,
                            "iter": self._iter
                        })
                except Exception as exc:
                    suuid = sample.uuid if sample else "UNKNOWN"
                    self.logger.error(
                        "[WorkerFactory] future exception consumed -> uuid={} | type -> {} | detail -> {!r}".format(
                            suuid, type(exc).__name__, exc
                        )
                    )
                    raise
                finally:
                    if sample is not None:
                        sample.close()
            self.tasks.difference_update(done)
        
        return DataConvert(data, dtype="listdict", v_column=self._vcolum, o_column=self._ocolum, device=self._device)
    
    def check_output_in_likelihood(self):
        likelihood_vars = {str(var) for var in getattr(self.loglike, "variables", set())}
        self._missing_likelihood_outputs = [str(var) for var in self._outputs if str(var) not in likelihood_vars]
        return len(self._missing_likelihood_outputs) == 0
        
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
        self.logger.info(f"DNN Plot -> saved -> {savepath}")
