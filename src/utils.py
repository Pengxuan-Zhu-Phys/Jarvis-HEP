#!/usr/bin/env python3
import pandas as pd 
from scipy.interpolate import interp1d
import numpy as np 

def get_interpolate_1D_function_from_config(config):
    """
        Generates a 1-dimensional interpolation function based on the provided configuration.
        Author: Pengxuan Zhu

        The configuration can specify the data directly through 'x_values' and 'y_values', or
        indirectly via a 'file' path to a CSV containing the data. It allows for optional
        logarithmic scaling on the x and/or y axes ('logX', 'logY') and supports various
        interpolation methods ('kind').

        Parameters:
        - config (dict): Configuration for the interpolation function. Expected keys include:
            - 'logX' (bool, optional): Apply logarithmic scaling to x-values. Defaults to False.
            - 'logY' (bool, optional): Apply logarithmic scaling to y-values. Defaults to False.
            - 'kind' (str, optional): Type of interpolation. Default is 'cubic'.
            - 'name' (str, optional): Name of the interpolation, for identification.
            - 'x_values' (list, optional): Directly provided list of x-values.
            - 'y_values' (list, optional): Directly provided list of y-values.
            - 'file' (str, optional): Path to a CSV file containing 'x' and 'y' columns.

        Returns:
        - function: A function that takes an x-value (or array of x-values) and returns
                    the interpolated y-value(s). If both 'logX' and 'logY' are True,
                    the function performs the interpolation in log-log space. If only 'logX'
                    is True, it interpolates in log-linear space. If only 'logY' is True,
                    it interpolates in linear-log space. If neither is True, it performs
                    a standard linear interpolation.

        Note:
        - If the configuration is incorrect or the interpolation cannot be set up (e.g., due
          to missing data points), the function returns None.
        - The interpolation function returned can handle both single values and numpy arrays
          of x-values.
    """
    from copy import deepcopy
    template = {
        "logX": False,
        "logY": False,
        "kind": "cubic",
        "name": None, 
        "x_values":  [],
        "y_values":  []
    }
    try:
        temp = deepcopy(template)
        temp.update(config)
        # print(temp)
        if "file" in temp:
            data = pd.read_csv(config['file'])
            temp['x_values'] = data['x']
            temp['y_values'] = data['y']
        if not temp['logX'] and not temp['logY']: 
            interp_func = interp1d(temp['x_values'], temp['y_values'], kind=temp['kind'], fill_value="extrapolate")
            return interp_func
        elif temp['logX'] and not temp['logY']:
            xx = np.log(temp['x_values'])
            interp_func_logX = interp1d(xx, temp['y_values'], kind=temp['kind'], fill_value='extrapolate')
            def interp_func(x_new):
                log_x_new = np.log(x_new)  # Convert x_new to log scale
                return interp_func_logX(log_x_new) 
            return interp_func
        elif not temp['logX'] and temp['logY']:
            yy = np.log(temp['y_values'])
            interp_func_logY = interp1d(temp['x_values'], yy, kind=temp['kind'], fill_value="extrapolate")
            def interp_func(x_new):
                return np.exp(interp_func_logY(x_new))
            return interp_func
        else: 
            xx = np.log(temp['x_values'])
            yy = np.log(temp['y_values'])
            interp_func_logXY = interp1d(xx, yy, kind=temp['kind'], fill_value="extrapolate")
            def interp_func(x_new):
                log_x_new = np.log(x_new)
                return np.exp(interp_func_logXY(log_x_new))
            return interp_func
    except:
        return None

def convert_hdf5_to_csv(snapshot_path, csv_path):
    import h5py, csv, json
    with h5py.File(snapshot_path, 'r') as hdf5_file:
        all_data = []
        # Iterate over datasets in the HDF5 file to collect data
        for dataset_name in hdf5_file:
            data = hdf5_file[dataset_name][()]
            # Assuming 'data' is stored as a binary string of serialized JSON
            json_data = json.loads(data.decode('utf-8'))
            all_data.append(json_data)
        # Assuming all JSON objects have the same structure (same keys)
        if all_data:
            fieldnames = list(all_data[0].keys())
            for data in all_data:
                if isinstance(data, dict):   
                    if "LogL" in data.keys():                     
                        if data['LogL'] == - np.inf:
                            fieldnames = data.keys()
                            continue
                        else:
                            fieldnames = data.keys()
                            break
                    else: 
                        fieldnames = data.keys()

            with open(csv_path, 'w', newline='') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
                for data in all_data:
                    if isinstance(data, dict):
                        writer.writerow(data)
