#!/usr/bin/env python3 

import numpy as np
import time
import random

class InitModule:
    version_number = "1.1.0"

    def __init__(self):
        self.run_params = {}

    def param_assign(self, run_params, lowerbounds, upperbounds, **kwargs):
        """
        Assign parameter values (defaults if not specified) to run_params.
        """
        run_params['D'] = len(lowerbounds)

        if len(upperbounds) != run_params['D']:
            raise ValueError('ERROR: parameter space bounds do not have the same dimensionality')
        if any(lowerbounds >= upperbounds):
            raise ValueError('ERROR: invalid parameter space bounds.')

        run_params['lowerbounds'] = np.array(lowerbounds)
        run_params['upperbounds'] = np.array(upperbounds)

        # Assign optional arguments
        nDerived = kwargs.get('nDerived', 0)
        run_params['D_derived'] = nDerived

        maxgen = kwargs.get('maxgen', 300)
        run_params['numgen'] = maxgen

        convthresh = kwargs.get('convthresh', 1e-3)
        run_params['convthresh'] = convthresh

        convsteps = kwargs.get('convsteps', 10)
        run_params['convsteps'] = convsteps
        run_params['improvements'] = np.ones(convsteps)

        doBayesian = kwargs.get('doBayesian', False)
        run_params['calcZ'] = doBayesian

        if run_params['calcZ']:
            run_params['maxNodePop'] = kwargs.get('maxNodePop', 1.9)
            run_params['tol'] = kwargs.get('Ztolerance', 0.01)
            run_params['numciv'] = kwargs.get('maxciv', 2000)
        else:
            run_params['numciv'] = kwargs.get('maxciv', 1)

        savecount = kwargs.get('savecount', 1)
        run_params['savefreq'] = savecount

        # Handle other logical settings
        run_params['disableIO'] = kwargs.get('disableIO', False)
        run_params['outputRaw'] = not run_params['disableIO'] and kwargs.get('outputRaw', True)
        run_params['outputSam'] = not run_params['disableIO'] and kwargs.get('outputSam', True)

        if run_params['disableIO']:
            if run_params['outputRaw']:
                print('WARNING: disableIO = true overrides outputRaw = true. No .raw file will be output.')
            if run_params['outputSam']:
                print('WARNING: disableIO = true overrides outputSam = true. No .sam file will be output.')

        # Seed for random number generator
        seed = kwargs.get('seed', -1)
        run_params['seed'] = seed if seed > 0 else -1

        # Setup DE-specific parameters
        run_params['DE'] = {
            'jDE': kwargs.get('jDE', True),
            'lambdajDE': kwargs.get('lambdajDE', True),
            'F': np.array(kwargs.get('F', [0.7])),  # Default value
            'Cr': kwargs.get('Cr', 0.9),  # Default crossover rate
            'lambda': kwargs.get('lambda', 0.0),
            'current': kwargs.get('current', False),
            'expon': kwargs.get('expon', False),
            'NP': max(10 * run_params['D'], 2 * len(run_params['DE']['F']) + 3),
            'removeDuplicates': kwargs.get('removeDuplicates', True),
            'bconstrain': kwargs.get('bndry', 1)  # Default boundary constraint
        }

        # Adjustments for NP and constraints
        if run_params['DE']['NP'] < 2 * len(run_params['DE']['F']) + 3:
            print('WARNING: NP specified is too small. Using smallest permitted NP...')
            run_params['DE']['NP'] = 2 * len(run_params['DE']['F']) + 3

        # Setup random seed initialization strategy
        run_params['init_population_strategy'] = kwargs.get('init_population_strategy', 0)
        run_params['discard_unfit_points'] = kwargs.get('discard_unfit_points', False)
        run_params['max_initialisation_attempts'] = kwargs.get('max_initialisation_attempts', 10000)
        run_params['max_acceptable_value'] = kwargs.get('max_acceptable_value', 1e6)

    def initialize(self, X, Xnew, run_params, func, fcall, quit_flag, accept):
        """
        Initializes the first generation of target vectors.
        """
        X['multiplicities'] = np.ones(run_params['DE']['NP'])
        if run_params['DE']['jDE']:
            Xnew['FjDE'] = self.init_FjDE(run_params['DE']['NP'])
            Xnew['CrjDE'] = self.init_CrjDE(run_params['DE']['NP'])
            if run_params['DE']['lambdajDE']:
                Xnew['lambdajDE'] = self.init_lambdajDE(run_params['DE']['NP'])

        for m in range(run_params['DE']['NP']):
            n = m
            max_attempts = 1 if run_params['init_population_strategy'] == 0 else run_params['max_initialisation_attempts']
            attempt_count = 0
            while attempt_count < max_attempts:
                Xnew['vectors'][m] = np.random.uniform(run_params['lowerbounds'], run_params['upperbounds'])
                Xnew['values'][m] = func(Xnew['vectors'][m], fcall, quit_flag, True, run_params['context'])
                if Xnew['values'][m] < run_params['max_acceptable_value'] or run_params['init_population_strategy'] == 0:
                    break
                attempt_count += 1

            if run_params['init_population_strategy'] >= 2 and attempt_count == max_attempts:
                raise ValueError('ERROR: init_population_strategy = 2 but could not find valid points within max_initialisation_attempts!')

    def init_FjDE(self, NP):
        return np.random.rand(NP)

    def init_CrjDE(self, NP):
        return np.random.rand(NP)

    def init_lambdajDE(self, NP):
        return np.random.rand(NP)

    def init_all_random_seeds(self, nprocs, input_seed):
        """
        Initializes random seeds for all processes.
        """
        if input_seed > 0:
            seed = input_seed
        else:
            seed = int(time.time()) + random.randint(0, 1000)

        np.random.seed(seed)
        return seed

# Example usage
init_module = InitModule()
X = {'vectors': np.zeros((10, 3)), 'values': np.zeros(10)}
Xnew = {'vectors': np.zeros((10, 3)), 'values': np.zeros(10)}

# Example run_params
run_params = {}
init_module.param_assign(run_params, lowerbounds=[0, 0, 0], upperbounds=[1, 1, 1])
init_module.initialize(X, Xnew, run_params, func=lambda x, fcall, quit, valid, context: np.sum(x), fcall=0, quit_flag=False, accept=0)
