#!/usr/bin/env python3

from concurrent.futures import ProcessPoolExecutor
import threading

class WorkerFactory:
    _instance = None
    _lock = threading.Lock()
    
    def __new__(cls, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super().__new__(cls)
                cls._instance.executor = ProcessPoolExecutor(max_workers=4)  # Initialize with desired number of workers
        return cls._instance
    def __init__(self):
        self.workflow = None  # Initialize with None

    def setup_workflow(self, workflow):
        self.workflow = workflow  # Assign the workflow

    def compute_likelihood(self, params):
        if not self.workflow:
            raise ValueError("Workflow is not set up.")
        # Use the workflow to compute the likelihood
        result = self.workflow.run(params)  # Assuming 'run' is the method that executes the workflow and computes the result
        return result

    def compute_likelihood(self, params):
        result = simulate_external_likelihood_calculation(params)
        return result