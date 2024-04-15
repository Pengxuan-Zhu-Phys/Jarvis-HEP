import threading
import time
import h5py
from queue import Queue, Empty
import csv 
import json 
from loguru import logger 

class GlobalHDF5Writer:
    def __init__(self, filepath, write_interval=15):
        super().__init__()
        self.filepath = filepath
        self.write_interval = write_interval
        self.data_queue = Queue()
        self.shutdown_event = threading.Event()
        self.writer_thread = threading.Thread(target=self._write_periodically)
        # It's safer to not rely solely on daemon threads for important cleanup.
        self.writer_thread.daemon = False
        self.logger = logger.bind(module="Jarvis-HEP.hdf5-Writter", to_console=True, Jarvis=True)

    def start(self):
        """Start the writer thread."""
        self.writer_thread.start()

    def add_data(self, data):
        """Add data to the queue to be written later."""
        serialized_data = json.dumps(data)
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
            with h5py.File(self.filepath, 'a') as f:
                # Example: adjust dataset creation and data writing as needed.
                for data in accumulated_data:
                    # Determine how to name and store each piece of data
                    # This is an example and needs to be adapted
                    dataset_name = f"data_{time.time()}"
                    f.create_dataset(dataset_name, data=data)
            self.logger.info(f"Global writer saved {len(accumulated_data)} data points to -> {self.filepath}.")
        else:
            self.logger.info("No data to write at this interval.")

    def stop(self):
        """Signal the writer thread to stop and wait for it to finish."""
        self.shutdown_event.set()
        self.writer_thread.join()
        self.logger.warning("Global HDF5 writer stopped.")


    def hdf5_to_csv(self, csv_path):
        """Convert the HDF5 file data to a CSV file with structured columns based on JSON keys."""
        with h5py.File(self.filepath, 'r') as hdf5_file:
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
                with open(csv_path, 'w', newline='') as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                    writer.writeheader()
                    for data in all_data:
                        if isinstance(data, dict):
                            writer.writerow(data)

        self.logger.warning(f"Converted HDF5 data to CSV at -> {csv_path}.")

