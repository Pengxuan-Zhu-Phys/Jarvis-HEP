#!/usr/bin/env python3 

# Let's translate `detypes.f90` into Python. We'll start by defining the necessary classes and structures.
# This will include the DE parameters (`deparams`), code parameters (`codeparams`), and population structures, 
# as well as the interface functions (`MinusLogLikeFunc` and `PriorFunc`).

from typing import Callable, List, Optional
import numpy as np

# Define the DE parameters class
class DEParams:
    def __init__(self, NP: int, F: Optional[np.ndarray] = None, Fsize: int = 1, lambda_: float = 0.0,
                 current: bool = False, Cr: float = 0.9, expon: bool = False, bconstrain: int = 1, jDE: bool = False,
                 lambdajDE: bool = False, removeDuplicates: bool = False):
        self.NP = NP  # Population size
        self.F = F if F is not None else np.array([0.7])  # Mutation scale factors
        self.Fsize = Fsize  # Number of entries in F
        self.lambda_ = lambda_  # Mutation scale factor for best-to-rand/current
        self.current = current  # Use current/best-to-current mutation
        self.Cr = Cr  # Crossover rate
        self.expon = expon  # Use exponential crossover (else binomial)
        self.bconstrain = bconstrain  # Boundary constraints for selection
        self.jDE = jDE  # Use self-adaptive Cr and F for DE
        self.lambdajDE = lambdajDE  # Use self-adaptive Cr, lambda, and F for DE
        self.removeDuplicates = removeDuplicates  # Weed out duplicate vectors within a generation

# Define the code parameters class
class CodeParams:
    def __init__(self, D: int, D_derived: int = 0, D_discrete: int = 0, maxNodePop: float = 1.9, tol: float = 0.01,
                 numciv: int = 1, numgen: int = 300, convthresh: float = 1e-3, convsteps: int = 10,
                 calcZ: bool = False, disableIO: bool = False, outputRaw: bool = False, outputSam: bool = False,
                 savefreq: int = 1, init_population_strategy: int = 0, discard_unfit_points: bool = False,
                 max_initialisation_attempts: int = 100, max_acceptable_value: float = 1.0, context: Optional[Callable] = None,
                 verbose: int = 1, seed: int = -1, convergence_criterion: int = 0):
        self.DE = DEParams(NP=10)  # Instantiate DE parameters with a default NP value
        self.D = D  # Dimension of parameter space
        self.D_derived = D_derived  # Dimension of derived space
        self.D_discrete = D_discrete  # Dimension of discrete parameter space
        self.discrete = np.array([])  # List of discrete dimensions
        self.partitionDiscrete = False  # Split population amongst discrete parameters
        self.numciv = numciv  # Max number of civilizations
        self.numgen = numgen  # Max number of generations
        self.convthresh = convthresh  # Convergence threshold
        self.convsteps = convsteps  # Number of steps for smoothing
        self.tol = tol  # Log-evidence tolerance
        self.maxNodePop = maxNodePop  # Max node population for posterior partitioning
        self.calcZ = calcZ  # Calculate log-evidence
        self.disableIO = disableIO  # Disable all IO
        self.outputRaw = outputRaw  # Output raw parameter samples
        self.outputSam = outputSam  # Output derived and rounded parameter samples
        self.savefreq = savefreq  # Frequency of saving progress
        self.init_population_strategy = init_population_strategy  # Initialisation strategy
        self.discard_unfit_points = discard_unfit_points  # Discard trial vectors with fitness above threshold
        self.max_initialisation_attempts = max_initialisation_attempts  # Max attempts to initialise population
        self.max_acceptable_value = max_acceptable_value  # Max acceptable fitness value
        self.context = context  # Context pointer for likelihood/prior function
        self.verbose = verbose  # Verbosity level
        self.seed = seed  # Random seed
        self.convergence_criterion = convergence_criterion  # Convergence criterion
        self.meanlike = 1e30  # Initialize large value for meanlike
        self.improvements = np.ones(convsteps)  # Array to store improvements over generations

# Define the population class
class Population:
    def __init__(self, NP: int, D: int, D_derived: int = 0):
        self.vectors = np.zeros((NP, D))  # Population vectors
        self.values = np.zeros(NP)  # Fitness values of population
        self.weights = np.ones(NP)  # Weights of population vectors
        self.multiplicities = np.ones(NP)  # Multiplicities of population vectors
        self.vectors_and_derived = np.zeros((NP, D + D_derived))  # Vectors with derived parameters
        self.FjDE = np.zeros(NP)  # Self-adaptive F parameters for jDE
        self.CrjDE = np.zeros(NP)  # Self-adaptive Cr parameters for jDE
        self.lambdajDE = np.zeros(NP)  # Self-adaptive lambda for jDE

# Define the likelihood and prior function interfaces as Python callables
MinusLogLikeFunc = Callable[[np.ndarray, int, bool, bool, Optional[Callable]], float]
PriorFunc = Callable[[np.ndarray, Optional[Callable]], float]

# Display the Python class structure
detypes_python_translation = """
Python Translation of detypes.f90:

- DEParams: Parameters related to Differential Evolution (DE), like population size, mutation factors, crossover rate.
- CodeParams: Holds all parameters for running the DE, including DEParams, bounds, initialization strategies, etc.
- Population: Contains the population vectors, fitness values, and other DE-related data.
- MinusLogLikeFunc: Abstract interface for the likelihood function, to be minimized (-ln(likelihood)).
- PriorFunc: Abstract interface for the prior function, used for Bayesian evidence calculation.
"""

detypes_python_translation
