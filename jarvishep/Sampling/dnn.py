#!/usr/bin/env python3 
import os, sys
import numpy as np
from jarvishep.Sampling.sampler import SamplingVirtial
import json
from jarvishep.sample import Sample
import concurrent.futures
import pandas as pd 
from jarvishep.dataconvert import DataConvert
from jarvishep.log_kv import format_two_column_log
from copy import deepcopy
# from torch.autograd import Variable
# from dnn import device

torch = None
nn = None


def _ensure_torch():
    global torch, nn
    if torch is None or nn is None:
        import importlib

        try:
            torch = importlib.import_module("torch")
            nn = importlib.import_module("torch.nn")
        except Exception as exc:
            raise RuntimeError(
                "DNN sampler requires PyTorch, but PyTorch could not be imported."
            ) from exc
    return torch, nn


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_json_safe(v) for v in value.tolist()]
    if torch is not None:
        tensor_type = getattr(torch, "Tensor", None)
        if tensor_type is not None and isinstance(value, tensor_type):
            return _json_safe(value.detach().cpu().tolist())
    return value


def _format_kv_block(title, items):
    return format_two_column_log(title, items)




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
        torch, nn = _ensure_torch()
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
        torch, nn = _ensure_torch()
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

    def supports_runtime_checkpointing(self) -> bool:
        return True

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
        _ensure_torch()
        self.set_torch_backend()        
        self._build_runtime_models()
        self.check_evaluation()
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

    def _build_runtime_models(self):
        torch, nn = _ensure_torch()
        self.regressors = [
            FeedforwardNN(
                input_size=self._dimensions,
                hidden_sizes=self._hidden_layers,
                output_size=1,
                otag=output_name,
                dropout_prob=0.1,
                learning_rate=self._learning_rate,
            )
            for output_name in self._outputs
        ]
        for reg in self.regressors:
            reg.logger = self.logger
            reg.set_device(self.device)

        self.classifier = Classifier(
            nn.Sequential(
                nn.Linear(self._dimensions, 100), nn.ReLU(),
                nn.Linear(100, 100), nn.ReLU(),
                nn.Linear(100, 100), nn.ReLU(),
                nn.Linear(100, 100), nn.ReLU(),
                nn.Linear(100, 1), nn.Sigmoid(),
            ),
            self._dimensions,
            learning_rate=self._learning_rate,
        )
        self.classifier.set_device(self.device)
        self.classifier.logger = self.logger

    @staticmethod
    def _sanitize_nested(value):
        if isinstance(value, dict):
            cleaned = {}
            for key, item in value.items():
                if key in {"logger", "handlers"}:
                    continue
                cleaned[key] = DNN._sanitize_nested(item)
            return cleaned
        if isinstance(value, list):
            return [DNN._sanitize_nested(item) for item in value]
        if isinstance(value, tuple):
            return tuple(DNN._sanitize_nested(item) for item in value)
        return value

    @staticmethod
    def _tensor_state_to_cpu(state_dict):
        return {k: v.detach().cpu() for k, v in state_dict.items()}

    def _export_model_state(self, model_obj):
        return {
            "x_mean": np.asarray(model_obj.x_mean).tolist(),
            "x_std": np.asarray(model_obj.x_std).tolist(),
            "model_state": self._tensor_state_to_cpu(model_obj.model.state_dict()),
            "optimizer_state": deepcopy(model_obj.optimizer.state_dict()),
        }

    def _restore_model_state(self, model_obj, payload):
        if not payload:
            return
        model_state = payload.get("model_state")
        if model_state:
            model_obj.model.load_state_dict(model_state)
        optimizer_state = payload.get("optimizer_state")
        if optimizer_state is not None:
            model_obj.optimizer.load_state_dict(optimizer_state)
        if "x_mean" in payload:
            model_obj.x_mean = np.asarray(payload["x_mean"])
        if "x_std" in payload:
            model_obj.x_std = np.asarray(payload["x_std"])

    def _export_sampler_state(self):
        torch_rng_state = None
        torch_mod = torch
        if torch_mod is not None and hasattr(torch_mod, "get_rng_state"):
            try:
                torch_rng_state = torch_mod.get_rng_state()
            except Exception:
                torch_rng_state = None
        return {
            "iter": int(self._iter),
            "device": self._device,
            "dataset": self.dataset.copy(deep=True),
            "df": self.df.copy(deep=True),
            "info": self._sanitize_nested(deepcopy(self.info)),
            "numpy_random_state": np.random.get_state(),
            "regressors": [self._export_model_state(reg) for reg in getattr(self, "regressors", []) or []],
            "classifier": self._export_model_state(self.classifier) if getattr(self, "classifier", None) is not None else None,
            "torch_rng_state": torch_rng_state,
        }

    def _import_sampler_state(self, payload):
        self._iter = int(payload.get("iter", self._iter or 0))
        self._device = payload.get("device", self._device)
        self.dataset = payload.get("dataset", self.dataset)
        self.df = payload.get("df", self.df)
        self.info = deepcopy(payload.get("info", self.info))
        np_state = payload.get("numpy_random_state")
        if np_state is not None:
            np.random.set_state(np_state)
        if not getattr(self, "regressors", None) or not getattr(self, "classifier", None):
            self._build_runtime_models()
        for reg, reg_payload in zip(self.regressors, payload.get("regressors", [])):
            self._restore_model_state(reg, reg_payload)
        self._restore_model_state(self.classifier, payload.get("classifier"))
        torch_state = payload.get("torch_rng_state")
        torch_mod = torch
        if torch_state is not None and torch_mod is not None and hasattr(torch_mod, "set_rng_state"):
            try:
                torch_mod.set_rng_state(torch_state)
            except Exception:
                pass

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
            
        if self.dataset is None or getattr(self.dataset, "empty", True):
            dataset = self.create_dataset(self._ninit)
            dpth = os.path.join(self.info["DNN_Simu"], "iter.0.csv")
            dataset.save(dpth)
            self.info['DNN_info']['iter_data'] = [dpth]
            self.dataset = dataset
            self._log_dataset_snapshot("init", self.dataset)
            self._write_run_info()
        else:
            self._log_dataset_snapshot(f"resume-iter-{self._iter}", self.dataset)
            
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
            self.persist_runtime_checkpoint(reason="dnn_iteration")
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
        torch, _nn = _ensure_torch()
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
        self._index = 0
        self.tasks = set()
        self.future_to_sample = {}
        data = []
        exhausted = False
        base_sample_cfg = self.info['sample']
        # sample_cls = Sample
        sample_cls = getattr(sys.modules.get(self.__class__.__module__), "Sample", Sample)

        while (not exhausted) or self.tasks:
            while not exhausted and len(self.tasks) < total_cores and self._index < num:
                sample = None
                try:
                    param = self.next_sample()
                except StopIteration:
                    exhausted = True
                    break

                try:
                    sample = sample_cls(param)
                    sample.set_config(self.build_sample_config(base_sample_cfg))
                    future = self.factory.submit_task(sample.info)
                    # Keep a direct fallback binding for smoke/fake futures.
                    setattr(future, "_jarvis_sample", sample)
                    self.tasks.add(future)
                    self.future_to_sample[future] = sample
                    sample = None
                except Exception:
                    if sample is not None:
                        try:
                            sample.close()
                        except Exception:
                            pass
                    raise

                self._index += 1
            if self._index >= num:
                exhausted = True

            if not self.tasks:
                continue

            timeout = self._checkpoint_seconds_until_due()
            done, _ = concurrent.futures.wait(
                self.tasks,
                timeout=timeout,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            if not done:
                self.persist_runtime_checkpoint(reason="checkpoint_heartbeat")
                continue
            for future in done:
                sample = self.future_to_sample.pop(future, None)
                if sample is None:
                    sample = getattr(future, "_jarvis_sample", None)
                try:
                    future.result()
                    if sample:
                        self.logger.warning(f"DNN create_dataset sample type -> {sample.__class__.__module__}.{sample.__class__.__name__}")
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
                        format_two_column_log(
                            "[WorkerFactory] future exception consumed",
                            [
                                ("uuid", suuid),
                                ("type", type(exc).__name__),
                                ("detail", repr(exc)),
                            ],
                        )
                    )
                    raise
                finally:
                    if sample is not None:
                        try:
                            sample.close()
                        finally:
                            setattr(sample.__class__, "close_calls", getattr(sample.__class__, "close_calls", 0))
            self.tasks.difference_update(done)
        
        return DataConvert(data, dtype="listdict", v_column=self._vcolum, o_column=self._ocolum, device=self._device)
    
    def check_output_in_likelihood(self):
        likelihood_vars = {str(var) for var in getattr(self.loglike, "variables", set())}
        self._missing_likelihood_outputs = [str(var) for var in self._outputs if str(var) not in likelihood_vars]
        return len(self._missing_likelihood_outputs) == 0
        
    def plot_status(self, Niter):
        from matplotlib import pyplot as plt

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
