#!/usr/bin/env python3 

import numpy as np
import pandas as pd


def _load_torch():
    import importlib

    try:
        return importlib.import_module("torch")
    except Exception as exc:
        raise RuntimeError(
            "PyTorch is required for tensor conversion but could not be imported."
        ) from exc

class DataConvert:
    def __init__(self, data_list, dtype="listdict", u_column=None, v_column=None, o_column=None, device="cpu"):
        self._vcolumn       = list(v_column or [])
        self._ocolumn       = list(o_column or [])
        self._logLcolumn    = "LogL"
        self._data          = None  
        self._device        = device
        if dtype == "listdict":
            self.init_data_to_dataframe(data_list)
            
    @property
    def data(self):
        return self._data

    def _to_tensor(self, columns):
        """Generic function to convert specified DataFrame columns to a PyTorch tensor."""
        torch = _load_torch()
        return torch.tensor(self.data[columns].values, dtype=torch.float32, device=self.device)

    @property
    def device(self):
        return self._device

    @property
    def valid(self):
        return self.LogL != -np.inf


    @property
    def v(self):
        return self.data[self._vcolumn].to_numpy()

    @property
    def LogL(self):
        return self.data[self._logLcolumn].to_numpy()

    # @property
    # def o(self):
    #     return self.data[self._ocolumn].to_numpy()
    
    def o(self, idx=None):
        if idx is None:
            return [self.data[col].to_numpy() for col in self._ocolumn]

        if isinstance(idx, int):
            if idx < 0 or idx >= len(self._ocolumn):
                raise IndexError(f"Observable index {idx} out of range [0, {len(self._ocolumn)-1}].")
            col = self._ocolumn[idx]
            return self.data[col].to_numpy()

        if isinstance(idx, str):
            if idx in self._ocolumn:
                return self.data[idx].to_numpy()
            else:
                raise KeyError(f"Observable name '{idx}' not found in columns {self._ocolumn}.")

        raise TypeError(f"Unsupported type for idx: {type(idx)} (expected None, int or str).")


    @property
    def v_tensor(self):
        return self._to_tensor(self._vcolumn)

    @property
    def o_tensor(self):
        return self._to_tensor(self._ocolumn)
    
    # def update(self, idxs):
    #     self._data = self._data.loc[idxs]

    def init_data_to_dataframe(self, data):
        """convert data -> pandas.DataFrame"""
        self._data = pd.DataFrame([{
            "uuid": d["uuid"], 
            **{key: d["v"][key] for key in self._vcolumn},  
            **{key: d['o'][key] for key in d['o'].keys()},   
            **{key: d['obs'][key] for key in d['obs'].keys()},
            "it": d.get("iter", -1)
            }
            for d in data]
        )
        self._data.set_index("uuid", inplace=True)

    def _ensure_dataframe(self, data=None):
        if data is None:
            if self._data is None:
                raise ValueError("No internal dataset available.")
            return self._data
        if isinstance(data, DataConvert):
            return data.data
        if isinstance(data, pd.DataFrame):
            return data
        if isinstance(data, list):
            return pd.DataFrame(data)
        raise TypeError(f"Unsupported data type for DataConvert conversion: {type(data)}")

    def to_numpy(self, data=None):
        """Convert current or provided data into NumPy array."""
        df = self._ensure_dataframe(data)
        return df.to_numpy()

    def to_listdict(self, data=None):
        """Convert current or provided data into list[dict]."""
        df = self._ensure_dataframe(data)
        return df.to_dict(orient="records")

    def update(self, new_data_convert):
        """
        Update the current DataConvert object with new DataConvert data.

        Args:
            new_data_convert (DataConvert): A new DataConvert object to merge.

        Returns:
            None (updates the existing object in place).
        """
        if not isinstance(new_data_convert, DataConvert):
            raise TypeError("new_data_convert must be an instance of DataConvert")

        # Concatenate existing and new data
        combined_data = pd.concat([self._data, new_data_convert.data], axis=0, ignore_index=False)

        # Drop duplicates based on uuid (keep the most recent entry)
        combined_data = combined_data[~combined_data.index.duplicated(keep='last')]

        # Update the internal dataframe
        self._data = combined_data
        
        
    def save(self, pth):
        self.data.to_csv(pth)
