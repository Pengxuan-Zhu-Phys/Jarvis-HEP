#!/usr/bin/env python3
import pandas as pd 
import numpy as np 
from jarvishep.observable_io import (
    csv_export_fieldnames_from_schema,
    format_csv_export_report,
    flatten_records_for_csv,
    load_schema,
    resolve_schema_path,
)

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
        from scipy.interpolate import interp1d

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

def convert_hdf5_to_csv(snapshot_path, csv_path, schema_path=None):
    import h5py, csv, json, os, warnings
    from jarvishep.observable_io import save_schema

    def _decode_one_record(raw, dataset_name):
        if isinstance(raw, (bytes, bytearray, np.bytes_)):
            text = bytes(raw).decode("utf-8")
        elif isinstance(raw, str):
            text = raw
        else:
            text = str(raw)
        try:
            return json.loads(text)
        except Exception as exc:
            raise ValueError(
                f"Failed to decode JSON record in dataset '{dataset_name}': {exc}"
            ) from exc

    def _iter_dataset_records(raw, dataset_name):
        if isinstance(raw, np.ndarray):
            if raw.ndim == 0:
                yield _decode_one_record(raw.item(), dataset_name)
                return
            for item in raw.reshape(-1):
                yield _decode_one_record(item, dataset_name)
            return
        yield _decode_one_record(raw, dataset_name)

    def _iter_records():
        with h5py.File(snapshot_path, 'r') as hdf5_file:
            for dataset_name in hdf5_file:
                data = hdf5_file[dataset_name][()]
                for record in _iter_dataset_records(data, dataset_name):
                    if isinstance(record, dict):
                        yield record

    if schema_path is None:
        schema_path = resolve_schema_path(snapshot_path)
    if str(schema_path).endswith(".schema.json"):
        pathroot = str(schema_path)[: -len(".schema.json")]
    else:
        pathroot = os.path.splitext(str(snapshot_path))[0]
    schema = load_schema(schema_path, pathroot)
    warn_msgs = schema.get("_warnings", [])
    if warn_msgs:
        for msg in warn_msgs:
            warnings.warn(f"[schema:convert_hdf5_to_csv] {msg}", RuntimeWarning, stacklevel=2)
        # Persist normalized schema so users can inspect the corrected result.
        save_schema(schema_path, schema)

    print(format_csv_export_report(schema))
    fieldnames = csv_export_fieldnames_from_schema(schema)
    row_count = 0
    for record in _iter_records():
        rows, _ = flatten_records_for_csv(
            records=[record],
            schema=schema,
            populate_name_map=False,
        )
        if not rows:
            continue
        row_count += 1

    if row_count == 0 or not fieldnames:
        return

    with open(csv_path, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for record in _iter_records():
            rows, _ = flatten_records_for_csv(
                records=[record],
                schema=schema,
                populate_name_map=False,
            )
            if rows:
                writer.writerow(rows[0])

def format_duration(seconds):
    """Formating into HH:MM:SS.msc format"""
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    millisec = seconds % 1
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}.{str(millisec)[2:5]}"
