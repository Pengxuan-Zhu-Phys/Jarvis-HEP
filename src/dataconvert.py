#!/usr/bin/env python3 

import uuid
import numpy as np
import pandas as pd
import torch

class DataConvert:
    def __init__(self, data_list, dtype="listdict", u_column=None, v_column=None, o_column=None, device="cpu"):
        self._uuid          = "uuid"  
        self._ucolumn       = u_column  
        self._vcolumn       = v_column  
        self._ocolumn       = o_column
        self._logLcolumn    = "logL"
        self._data          = None  
        self._device        = device
        if dtype == "listdict":
            self.init_data_to_dataframe(data_list)

    @property
    def data(self):
        return self._data

    def _to_tensor(self, columns):
        """Generic function to convert specified DataFrame columns to a PyTorch tensor."""
        return torch.tensor(self.data[columns].values, dtype=torch.float32, device=self.device)

    @property
    def device(self):
        return self._device

    @property
    def u(self):
        return self.data[self._ucolumn].to_numpy()

    @property
    def v(self):
        return self.data[self._vcolumn].to_numpy()

    @property
    def logL(self):
        return self.data[self._logLcolumn].to_numpy()

    @property
    def o(self):
        return self.data[self._ocolumn].to_numpy()
    
    @property
    def u_tensor(self):
        return self._to_tensor(self._ucolumn)

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
            **dict(zip(self._ucolumn, d["u"])),  
            **{key: d["v"][key] for key in self._vcolumn},  
            **dict(zip(self._ocolumn, d["o"])),   
            **{key: d['obs'][key] for key in d['obs'].keys()},
            "it": d.get("iter", -1)
            }
            for d in data]
        )
        self._data.set_index("uuid", inplace=True)

    def to_numpy(self, data):
        """converting to numpy array"""
        df = self._ensure_dataframe(data)
        return df[[self.uuid_column] + self.input_columns + self.output_columns].values

    def to_listdict(self, data):
        """converting to list[dict] """
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