#!/usr/bin/env python3 

import pickle
import pandas as pd
from datetime import datetime

class TPMCMCData:
    def __init__(self):
        # List of sample dictionaries.
        self.samples = []
        # Dictionary with chain-level status: key = chain_id, value = info dict.
        self.chains_info = {}
        # Global iteration counter.
        self.global_iter = 0

    def append_sample(self, chain_id, chain_iter, param, loglikelihood, timestamp=None):
        """Append a new sample row."""
        if timestamp is None:
            timestamp = datetime.now()
            sample_row = {
                "global_iter": self.global_iter,
                "chain_id": chain_id,
                "chain_iter": chain_iter,
                "param": param.tolist() if isinstance(param, np.ndarray) else param,
                "loglikelihood": loglikelihood,
                "timestamp": timestamp
            }
        self.samples.append(sample_row)
        self.global_iter += 1

    def update_chain(self, chain_id, current_param, iterations, last_loglikelihood, proposal_scale=None):
        """Update (or create) the status of a chain."""
        self.chains_info[chain_id] = {
            "current_param": current_param,
            "iterations": iterations,
            "last_loglikelihood": last_loglikelihood,
            "proposal_scale": proposal_scale
        }

    def save_pickle(self, filename):
        """Save the whole TPMCMCData object as a pickle file."""
        with open(filename, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load_pickle(cls, filename):
        """Load a TPMCMCData object from a pickle file."""
        with open(filename, "rb") as f:
            return pickle.load(f)

    def to_dataframe(self):
        """Convert the samples to a Pandas DataFrame for fast CSV export or analysis."""
        return pd.DataFrame(self.samples)

    def save_csv(self, filename):
        """Save all sample data as a CSV file."""
        df = self.to_dataframe()
        df.to_csv(filename, index=False)

    def save_csv_split(self, filename_prefix, rows_per_file=10000):
        """
        Optionally, split the output into multiple CSV files if the dataset is huge.
        Each file will have at most `rows_per_file` rows.
        """
        df = self.to_dataframe()
        num_files = (len(df) + rows_per_file - 1) // rows_per_file
        for i in range(num_files):
            start_row = i * rows_per_file
            end_row = start_row + rows_per_file
            sub_df = df.iloc[start_row:end_row]
            sub_df.to_csv(f"{filename_prefix}_{i+1}.csv", index=False)