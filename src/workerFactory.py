
class WorkerFactory:
    # Previous implementation
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



from concurrent.futures import ProcessPoolExecutor
import threading
import dynesty  # Assuming dynesty is already installed

class WorkerFactory:
    _instance = None
    _lock = threading.Lock()
    
    def __new__(cls, *args, **kwargs):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super().__new__(cls)
                cls._instance.executor = ProcessPoolExecutor(max_workers=4)  # Initialize with desired number of workers
        return cls._instance

    def compute_likelihood(self, params):
        result = simulate_external_likelihood_calculation(params)
        return result

def simulate_external_likelihood_calculation(params):
    # Placeholder for the external likelihood computation
    # This would be where you call your simulation or external process
    return sum(params)  # Example computation

class Core:
    def __init__(self):
        self.worker_factory = WorkerFactory()
        self.sampler = dynesty.DynamicNestedSampler(
            self.likelihood,
            prior_transform,
            ndim,
            # Other dynesty sampler settings...
        )
        
    def likelihood(self, params):
        # Wrapper around the WorkerFactory's likelihood computation
        return self.worker_factory.compute_likelihood(params)

    def run(self):
        self.sampler.run_nested()
        # Additional logic for handling sampler results...

def prior_transform(u):
    # Placeholder for the prior transform function required by Dynesty
    # This transforms from unit cube to the parameter space
    return u  # Example identity transform

ndim = 2  # Number of dimensions in the parameter space

# Example usage
core = Core()
core.run()



class Core:
    def __init__(self):
        self.worker_factory = WorkerFactory()
        self.workflow = Workflow()  # Assuming Workflow class is defined elsewhere

        # Set up workflow with modules here
        # self.workflow.add_module(Module('ModuleName', inputs=..., outputs=...))

        self.worker_factory.setup_workflow(self.workflow)  # Pass the workflow to the worker factory

        self.sampler = dynesty.DynamicNestedSampler(
            self.likelihood,
            prior_transform,
            ndim,
            # Other dynesty sampler settings...
        )

    # Remainder of the Core class implementation...
