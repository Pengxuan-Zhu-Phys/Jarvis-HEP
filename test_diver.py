#!/usr/bin/env python3  

import numpy as np
import time
from collections import deque
import uuid


# Constants for convergence criteria
MEAN_IMPROVEMENT = 0
CHECK_POP_RES = False  # Whether to check population resolution
TAU = 0.1  # jDE control parameter from Brest et al. 2006


class RunParams:
    def __init__(self, D, NP, convergence_criterion=MEAN_IMPROVEMENT, convsteps=10, convthresh=1e-6, verbose=1):
        self.D = D  # Dimensionality
        self.NP = NP  # Population size
        self.convergence_criterion = convergence_criterion
        self.convsteps = convsteps
        self.convthresh = convthresh
        self.verbose = verbose
        self.meanlike = np.finfo(float).max
        self.improvements = [1.0] * self.convsteps  # Initialize with ones

def gencrossover(X, V, n, run_params):
    """
    General crossover function that generates the trial vector U.

    Parameters:
        X : Population
            Current population.
        V : ndarray
            Donor vector.
        n : int
            Index of the current target vector.
        run_params : RunParams
            Parameters for the DE algorithm.

    Returns:
        U : ndarray
            Trial vector.
        trialCr : float
            The crossover probability used (for jDE).
    """
    if run_params.DE.jDE:
        trialCr = newCr(X, n, run_params)
        U = bincrossover(X, V, n, run_params, trialCr)
    elif run_params.DE.expon:
        U = expcrossover(X, V, n, run_params)
        trialCr = run_params.DE.Cr
    else:
        U = bincrossover(X, V, n, run_params)
        trialCr = run_params.DE.Cr

    return U, trialCr


def bincrossover(X, V, n, run_params, trialCr=None):
    """
    Binomial crossover to create trial vectors.

    Parameters:
        X : Population
            Current population.
        V : ndarray
            Donor vector.
        n : int
            Index of the current target vector.
        run_params : RunParams
            Parameters for the DE algorithm.
        trialCr : float, optional
            Self-adapting Cr for jDE.

    Returns:
        U : ndarray
            Trial vector created.
    """
    D = run_params.D
    if trialCr is not None:
        Cr = trialCr
    else:
        Cr = run_params.DE.Cr

    # Choose a random index jrand in [0, D-1]
    jrand = np.random.randint(0, D)
    # Generate random numbers for each dimension
    randj = np.random.rand(D)

    # Initialize U with target vector
    U = X.vectors[n].copy()

    # Apply crossover
    crossover_mask = (randj <= Cr)
    U[crossover_mask] = V[crossover_mask]
    # Ensure that at least one parameter is copied from V
    U[jrand] = V[jrand]

    return U


def expcrossover(X, V, n, run_params):
    """
    Exponential crossover to create trial vectors.

    Parameters:
        X : Population
            Current population.
        V : ndarray
            Donor vector.
        n : int
            Index of the current target vector.
        run_params : RunParams
            Parameters for the DE algorithm.

    Returns:
        U : ndarray
            Trial vector created.
    """
    D = run_params.D
    Cr = run_params.DE.Cr
    U = X.vectors[n].copy()

    # Choose a random starting index j in [0, D-1]
    j = np.random.randint(0, D)
    L = 0  # Length of crossover

    # Determine the length of the crossover L
    for _ in range(D):
        rand = np.random.rand()
        if rand <= Cr:
            L += 1
        else:
            break

    # Apply crossover from index j to j+L (with wrap-around)
    for k in range(L):
        idx = (j + k) % D
        U[idx] = V[idx]

    return U

def newCr(X, n, run_params):
    """
    Generates a new Cr value for self-adaptive DE (jDE).

    Parameters:
        X : Population
            Current population.
        n : int
            Index of the current target vector.
        run_params : RunParams
            Parameters for the DE algorithm.

    Returns:
        newCr : float
            The new crossover probability.
    """
    rand = np.random.rand()
    if rand < TAU:
        newCr = np.random.rand()
    else:
        newCr = X.CrjDE[n]  # Use Cr from previous generation

    return newCr

def init_CrjDE(size):
    """
    Initializes Cr values for the population in jDE.

    Parameters:
        size : int
            The size of the population.

    Returns:
        CrjDE : ndarray
            Initialized Cr values.
    """
    CrjDE = np.random.rand(size)
    return CrjDE



def converged(X, run_params):
    if run_params.verbose >= 3:
        print('  Checking convergence...')

    if run_params.convergence_criterion == MEAN_IMPROVEMENT:
        is_converged = check_SFIM(X, run_params)
    else:
        is_converged = False  # No convergence criterion used

    if is_converged:
        if run_params.verbose >= 3:
            print('  Converged.')
    else:
        if run_params.verbose >= 3:
            print('  Not converged.')

    if CHECK_POP_RES:
        is_converged = check_population_resolution(X, run_params, is_converged)

    return is_converged

def init_convergence(run_params):
    if run_params.convergence_criterion == MEAN_IMPROVEMENT:
        run_params.meanlike = np.finfo(float).max
        run_params.improvements = [1.0] * run_params.convsteps  # Initialize with ones


def check_SFIM(X, run_params):
    curval = np.mean(X.values)
    inf_threshold = 0.001 * np.finfo(float).max

    if curval > inf_threshold:
        run_params.meanlike = inf_threshold
        fracdiff = 1.0
    else:
        fracdiff = 1.0 - curval / run_params.meanlike
        run_params.meanlike = curval

    # Update the improvements list
    run_params.improvements = run_params.improvements[1:] + [fracdiff]

    sfim = sum(run_params.improvements) / run_params.convsteps

    if run_params.verbose >= 3:
        print('  Smoothed fractional improvement of the mean =', sfim)

    is_converged = sfim < run_params.convthresh
    return is_converged



class DiverIterator:
    def __init__(self, lowerbounds, upperbounds, NP=100, F=0.8, Cr=0.9, maxgen=1000, seed=None, convsteps=10, convthresh=1e-6, verbose=1):
        """
        Differential Evolution iterator that yields trial vectors with UUIDs and accepts likelihood values asynchronously.
        
        Parameters:
            lowerbounds : array_like
                Lower boundaries of the parameter space.
            upperbounds : array_like
                Upper boundaries of the parameter space.
            NP : int
                Population size.
            F : float
                Differential weight.
            Cr : float
                Crossover probability.
            maxgen : int
                Maximum number of generations.
            seed : int, optional
                Random seed.
        """
        self.lowerbounds = np.array(lowerbounds)
        self.upperbounds = np.array(upperbounds)
        self.NP = NP
        self.F = F
        self.Cr = Cr
        self.maxgen = maxgen
        self.seed = seed
        self.D = len(lowerbounds)  # Dimensionality of the problem
        self.gen = 0  # Current generation
        self.pop = None  # Current population
        self.pop_fitness = np.full(NP, np.inf)  # Initialize fitness values
        self.current_idx = 0  # Index of the individual being processed
        self.uuid_map = {}  # Map from UUIDs to indices in the population
        self.pending_uuids = set()  # Set of UUIDs awaiting likelihood values
        self.initialize_population()
        self.best_vector = None
        self.best_fitness = np.inf
        if seed is not None:
            np.random.seed(seed)

        self.run_params = RunParams(D=self.D, NP=self.NP, convsteps=convsteps, convthresh=convthresh, verbose=verbose)
        init_convergence(self.run_params)
        self.converged = False
        
    def initialize_population(self):
        # Initialize population within bounds
        self.pop = np.random.uniform(self.lowerbounds, self.upperbounds, (self.NP, self.D))
        # For initial population, we need to assign UUIDs and wait for fitness values
        self.uuid_map = {}
        self.pending_uuids = set()
        for idx in range(self.NP):
            uid = uuid.uuid4()
            self.uuid_map[uid] = idx
            self.pending_uuids.add(uid)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        # If all generations are completed, stop iteration
        if self.gen >= self.maxgen:
            raise StopIteration
        
        # If there are pending UUIDs, we cannot proceed until likelihood values are supplied
        if self.pending_uuids:
            raise StopIteration("Awaiting likelihood values for UUIDs: {}".format(self.pending_uuids))
        
        # Start processing individuals in the current generation
        if self.current_idx < self.NP:
            # Generate trial vector for individual self.current_idx
            target = self.pop[self.current_idx]
            mutant = self.mutate(self.current_idx)
            trial = self.crossover(target, mutant)
            # Ensure trial vector is within bounds
            trial = np.clip(trial, self.lowerbounds, self.upperbounds)
            # Assign a UUID to the trial vector
            uid = uuid.uuid4()
            self.uuid_map[uid] = self.current_idx
            self.pending_uuids.add(uid)
            self.current_idx += 1
            # Yield the trial vector along with its UUID
            return uid, trial
        else:
            # All trial vectors for this generation have been generated
            # Wait until all likelihood values have been received
            if self.pending_uuids:
                raise StopIteration("Awaiting likelihood values for UUIDs: {}".format(self.pending_uuids))
            else:
                # Prepare for next generation
                self.gen += 1
                self.current_idx = 0
                # Proceed to next generation
                return self.__next__()
    
    def mutate(self, idx):
        # Implement DE/rand/1 mutation strategy
        idxs = [i for i in range(self.NP) if i != idx]
        a, b, c = np.random.choice(idxs, 3, replace=False)
        mutant = self.pop[a] + self.F * (self.pop[b] - self.pop[c])
        return mutant
    
    def crossover(self, target, mutant):
        # Implement binomial crossover
        crossover_mask = np.random.rand(self.D) < self.Cr
        trial = np.where(crossover_mask, mutant, target)
        return trial
    
    def supply_likelihood(self, uid, fitness):
        """
        Supply the likelihood (fitness) value for a trial vector identified by its UUID.
        
        Parameters:
            uid : UUID
                The UUID of the trial vector.
            fitness : float
                The fitness value corresponding to the trial vector.
        """
        if uid not in self.uuid_map:
            raise ValueError("UUID {} not recognized.".format(uid))
        idx = self.uuid_map[uid]
        # Compare with target vector's fitness
        target_fitness = self.pop_fitness[idx]
        if fitness < target_fitness:
            # Update the population with the new trial vector
            # Note: Since the trial vector was yielded, we need to store it temporarily
            # For simplicity, we'll store it in a temporary map
            if not hasattr(self, 'trial_vectors'):
                self.trial_vectors = {}
            self.trial_vectors[idx] = (uid, fitness)
        # Remove the UUID from pending set
        self.pending_uuids.remove(uid)
        # If all likelihoods for the current generation are received, perform selection
        if not self.pending_uuids and self.current_idx >= self.NP:
            self.selection()
            # Reset for the next generation
            self.current_idx = 0
            self.uuid_map.clear()
            if hasattr(self, 'trial_vectors'):
                del self.trial_vectors
    
    def selection(self):
        # Perform selection based on the fitness values
        if hasattr(self, 'trial_vectors'):
            for idx, (uid, trial_fitness) in self.trial_vectors.items():
                # Since we didn't store the trial vector, we need to regenerate it
                # Alternatively, you can store the trial vectors along with their UUIDs
                target = self.pop[idx]
                mutant = self.mutate(idx)
                trial = self.crossover(target, mutant)
                trial = np.clip(trial, self.lowerbounds, self.upperbounds)
                # Update the population
                self.pop[idx] = trial
                self.pop_fitness[idx] = trial_fitness
                # Update best vector if needed
                if trial_fitness < self.best_fitness:
                    self.best_vector = trial.copy()
                    self.best_fitness = trial_fitness
        # Else, no trial vectors improved the population
    
    def get_best(self):
        """
        Get the best vector and its fitness found so far.
        
        Returns:
            best_vector : array_like
                Best parameter vector found.
            best_fitness : float
                Fitness value of the best vector.
        """
        return self.best_vector, self.best_fitness

# Define lower and upper bounds
lowerbounds = [-5, -5]
upperbounds = [5, 5]

# Create an instance of the DiverIterator with convergence parameters
diver = DiverIterator(
    lowerbounds=lowerbounds,
    upperbounds=upperbounds,
    NP=10,
    maxgen=1000,
    seed=42,
    convsteps=10,
    convthresh=1e-6,
    verbose=2
)

# Simulate asynchronous likelihood evaluation
def simulate_likelihood(trial_vector):
    return np.sum(trial_vector**2)

# Main loop to run the DE algorithm
while True:
    try:
        uid, trial_vector = next(diver)
        fitness = simulate_likelihood(trial_vector)
        diver.supply_likelihood(uid, fitness)
    except StopIteration:
        break

# Get the best solution found
best_vector, best_fitness = diver.get_best()
print("Best vector found:", best_vector)
print("Fitness of best vector:", best_fitness)

