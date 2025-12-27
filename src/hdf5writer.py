from __future__ import annotations

import threading
import time
import h5py
from queue import Queue, Empty
import csv 
import json 
from loguru import logger 
import numpy as np 
import sympy
import os 

class GlobalHDF5Writer:
    def __init__(self, dbinfo, write_interval=15):
        super().__init__()
        self._info_lock = threading.Lock()
        self.initialize(dbinfo)
        # self.filepath = dbinfo['path']
        self.write_interval = write_interval
        self.data_queue = Queue()
        self.shutdown_event = threading.Event()
        self.writer_thread = threading.Thread(target=self._write_periodically)
        # It's safer to not rely solely on daemon threads for important cleanup.
        self.writer_thread.daemon = False
        self.logger = logger.bind(module="Jarvis-HEP.hdf5-Writter", to_console=True, Jarvis=True)
        self.max_size = 1024 * 1024 * 100  # Max for 100 MB

    def initialize(self, dbinfo):
        pthroot, pthext = os.path.splitext(dbinfo['path'])
        self.infos = {
            "info":              dbinfo['info'],
            "activeNO":          0,
            "totalDB":           1,
            "converted":         [],
            "pending_converted": [],  # csv paths we failed to convert earlier
            "errors":            [],  # recent conversion/write errors
            "paths":             [],
            "pathroot":          pthroot,
            "pathext":           pthext,
            "active path":       "".join([pthroot, ".0", pthext])
        }
        self._save_infos()

    def _atomic_write_json(self, path: str, data: dict) -> None:
        tmp_path = path + ".tmp"
        with open(tmp_path, "w") as f:
            json.dump(data, f, indent=4)
            f.flush()
        os.replace(tmp_path, path)

    def _save_infos(self) -> None:
        with self._info_lock:
            self._atomic_write_json(self.infos["info"], self.infos)

    def start(self):
        """Start the writer thread."""
        self.writer_thread.start()


    def add_data(self, data):
        """Add data to the queue to be written later."""
        data = convert(data)
        # for kk, vv in data.items():
        #     print(kk, vv, type(vv))
        # print("hdf5Writer Line 30 -> ", data)
        serialized_data = json.dumps(data)
        # print("hdf5Writer Line 32 -> ", serialized_data)
        self.data_queue.put(serialized_data)

    def _write_periodically(self):
        """Periodically write accumulated data to HDF5."""
        while not self.shutdown_event.is_set():
            time.sleep(self.write_interval)
            self._write_data_to_hdf5()

        # Final write to ensure all data is flushed when shutting down
        self._write_data_to_hdf5()

    def _write_data_to_hdf5(self):
        """Write accumulated data to the HDF5 file."""
        accumulated_data = []
        while not self.data_queue.empty():
            try:
                data = self.data_queue.get_nowait()
                accumulated_data.append(data)
            except Empty:
                break
        
        if accumulated_data:
            if os.path.exists(self.infos['active path']):
                current_file_size = os.path.getsize(self.infos['active path'])
            else:
                current_file_size = 0  # New file has size 0

            if current_file_size >= self.max_size:
                # Rotate: try to convert the current active file to CSV (best-effort)
                try:
                    self.hdf5_to_csv()
                except Exception as e:
                    # hdf5_to_csv is already best-effort, but be defensive here too.
                    self.logger.warning(f"hdf5_to_csv failed during rotation: {e}")

                self.infos['paths'].append(self.infos['active path'])
                self.infos['activeNO'] += 1
                self.infos['totalDB'] += 1

                new_file_path = "".join([
                    self.infos['pathroot'], ".{}".format(self.infos['activeNO']), self.infos['pathext']
                ])
                self.infos['active path'] = new_file_path
                self.logger.warning(f"Creating new HDF5 file due to size limit: {new_file_path}")
                self._save_infos()

            with h5py.File(self.infos['active path'], 'a') as f:
                # Example: adjust dataset creation and data writing as needed.
                for data in accumulated_data:
                    # Determine how to name and store each piece of data
                    # This is an example and needs to be adapted
                    dataset_name = f"data_{time.time()}"
                    f.create_dataset(dataset_name, data=data)
            self.logger.info(f"Global writer saved {len(accumulated_data)} data points to -> {self.infos['active path']}.")
        else:
            self.logger.info("No data to write at this interval.")

    def stop(self):
        """Signal the writer thread to stop and wait for it to finish."""
        self.shutdown_event.set()
        self.writer_thread.join()
        self.logger.warning("Global HDF5 writer stopped.")

        # Optional best-effort conversion at shutdown (never raise)
        try:
            self.hdf5_to_csv()
        except Exception as e:
            self.logger.warning(f"Final hdf5_to_csv failed (ignored): {e}")

    def hdf5_to_csv(self) -> bool:
        """Convert the active HDF5 file to CSV (best-effort).

        Returns:
            bool: True if conversion succeeded, False otherwise.
        """
        csv_path = "".join([self.infos['pathroot'], ".{}".format(self.infos['activeNO']), ".csv"])

        # Already converted
        if csv_path in self.infos.get("converted", []):
            return True

        # Attempt conversion (do NOT raise)
        try:
            with h5py.File(self.infos['active path'], 'r') as hdf5_file:
                all_data = []
                for dataset_name in hdf5_file:
                    data = hdf5_file[dataset_name][()]
                    # Stored as serialized JSON string (bytes)
                    if isinstance(data, (bytes, bytearray)):
                        json_data = json.loads(data.decode('utf-8'))
                    else:
                        json_data = json.loads(str(data))
                    all_data.append(json_data)

            if not all_data:
                # Nothing to convert; treat as success but still mark as converted to avoid repeated work.
                self.infos.setdefault("converted", []).append(csv_path)
                self._save_infos()
                return True

            # Determine fieldnames robustly
            fieldnames = list(all_data[0].keys()) if isinstance(all_data[0], dict) else []
            for item in all_data:
                if isinstance(item, dict):
                    # Skip pathological -inf row if present, otherwise take the first normal row
                    if item.get('LogL', None) == -np.inf:
                        continue
                    fieldnames = list(item.keys())
                    break

            if not fieldnames:
                # Fallback: union of keys
                keys = set()
                for item in all_data:
                    if isinstance(item, dict):
                        keys.update(item.keys())
                fieldnames = sorted(keys)

            with open(csv_path, 'w', newline='') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
                for item in all_data:
                    if isinstance(item, dict):
                        writer.writerow(item)

            self.logger.warning(f"Converted HDF5 data to CSV at -> {csv_path}.")
            self.infos.setdefault("converted", []).append(csv_path)
            # Remove from pending if it was pending
            if csv_path in self.infos.get("pending_converted", []):
                self.infos["pending_converted"].remove(csv_path)
            self._save_infos()
            return True

        except (BlockingIOError, OSError) as e:
            # Most common: file locking / temporarily unavailable
            msg = f"CSV conversion skipped (file busy/locked): active={self.infos['active path']} err={e}"
            self.logger.warning(msg)
            self.infos.setdefault("pending_converted", [])
            if csv_path not in self.infos["pending_converted"]:
                self.infos["pending_converted"].append(csv_path)
            self.infos.setdefault("errors", []).append(msg)
            # Keep errors list bounded
            self.infos["errors"] = self.infos["errors"][-50:]
            self._save_infos()
            return False

        except Exception as e:
            msg = f"CSV conversion failed (ignored): active={self.infos['active path']} err={e}"
            self.logger.warning(msg)
            self.infos.setdefault("pending_converted", [])
            if csv_path not in self.infos["pending_converted"]:
                self.infos["pending_converted"].append(csv_path)
            self.infos.setdefault("errors", []).append(msg)
            self.infos["errors"] = self.infos["errors"][-50:]
            self._save_infos()
            return False


""" Custom JSON encoder, converte the numpy number format into python float. """

class NumpyEncoder(json.JSONEncoder):
    """ Custom encoder for numpy data types """
    def default(self, obj):
        if isinstance(obj, (np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)): # This includes 0-d numpy arrays.
            return obj.tolist()  # or item()
        return json.JSONEncoder.default(self, obj)
    
def convert(data):
    if isinstance(data, dict):
        return {k: convert(v) for k, v in data.items()}
    elif isinstance(data, (np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
        return int(data)
    elif isinstance(data, (np.float16, np.float32, np.float64)):
        return float(data)
    elif isinstance(data, np.str_):
        return str(data)
    elif isinstance(data, (np.ndarray,)): # This includes 0-d numpy arrays.
        return data.item()
    elif isinstance(data, (sympy.Float, sympy.Integer)):
        return float(data)
    return data
