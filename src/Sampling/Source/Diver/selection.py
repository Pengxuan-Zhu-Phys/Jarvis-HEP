#!/usr/bin/env python3 

import numpy as np

class Selection:
    def __init__(self):
        pass

    def selector(self, X, Xnew, U, trialF, triallambda, trialCr, m, n, run_params, func, fcall, quit_flag):
        """Selector function that compares trial vectors with current population and accepts or rejects them."""
        accept = 0
        acceptable_trial_vector = False

        # Determine whether the trial vector is within bounds or not
        if any(U > run_params['upperbounds']) or any(U < run_params['lowerbounds']):
            if run_params['bconstrain'] == 1:  # 'brick wall' case
                trialvector = X['vectors'][n]
                validvector = False
            elif run_params['bconstrain'] == 2:  # Randomly re-initialize
                trialvector = np.random.random(size=len(U)) * (run_params['upperbounds'] - run_params['lowerbounds']) + run_params['lowerbounds']
                validvector = True
            elif run_params['bconstrain'] == 3:  # Reflection case
                trialvector = U
                trialvector = np.where(U > run_params['upperbounds'], run_params['upperbounds'] - (U - run_params['upperbounds']), trialvector)
                trialvector = np.where(U < run_params['lowerbounds'], run_params['lowerbounds'] + (run_params['lowerbounds'] - U), trialvector)
                validvector = True
            else:
                trialvector = U
                validvector = True
        else:
            trialvector = U
            validvector = True

        trialderived = self.roundvector(trialvector, run_params)
        trialvalue = func(trialderived, fcall, quit_flag, validvector, run_params['context'])

        if not validvector:
            trialvalue = np.inf  # If vector is not valid, assign a large value

        acceptable_fitness = (trialvalue <= run_params['max_acceptable_value'])
        acceptable_trial_vector = acceptable_fitness or not run_params.get('discard_unfit_points', True)

        if acceptable_trial_vector:
            if trialvalue <= X['values'][n]:
                Xnew['vectors'][m] = trialvector
                Xnew['vectors_and_derived'][m] = trialderived
                Xnew['values'][m] = trialvalue
                if run_params['jDE']:
                    Xnew['FjDE'][m] = trialF
                    Xnew['CrjDE'][m] = trialCr
                    if run_params['lambdajDE']:
                        Xnew['lambdajDE'][m] = triallambda
                accept += 1
            else:
                Xnew['vectors'][m] = X['vectors'][n]
                Xnew['vectors_and_derived'][m] = X['vectors_and_derived'][n]
                Xnew['values'][m] = X['values'][n]
                if run_params['jDE']:
                    Xnew['FjDE'][m] = X['FjDE'][n]
                    Xnew['CrjDE'][m] = X['CrjDE'][n]
                    if run_params['lambdajDE']:
                        Xnew['lambdajDE'][m] = X['lambdajDE'][n]
        return accept, acceptable_trial_vector

    def replace_generation(self, X, Xnew, run_params, func, fcall, quit_flag, accept, init):
        """Replace the old generation (X) with the new generation (Xnew)."""
        allvecs = Xnew['vectors']

        if run_params['removeDuplicates']:
            self.remove_duplicate_vectors(X, Xnew, run_params, allvecs, func, fcall, quit_flag, accept, init)

        X['vectors'] = allvecs
        X['values'] = Xnew['values']
        X['vectors_and_derived'] = Xnew['vectors_and_derived']

        if run_params['jDE']:
            X['FjDE'] = Xnew['FjDE']
            X['CrjDE'] = Xnew['CrjDE']
            if run_params['lambdajDE']:
                X['lambdajDE'] = Xnew['lambdajDE']

    def remove_duplicate_vectors(self, X, Xnew, run_params, allvecs, func, fcall, quit_flag, accept, init):
        """Remove duplicate vectors to maintain population diversity."""
        for k in range(run_params['NP'] - 1):
            if any(allvecs[k, 0] == allvecs[k + 1:, 0]):
                for kmatch in range(k + 1, run_params['NP']):
                    if np.all(allvecs[k] == allvecs[kmatch]):
                        self.replace_vector(Xnew, allvecs, X, run_params, func, kmatch, fcall, quit_flag, accept, revert=False)

    def replace_vector(self, Xnew, allvecs, X, run_params, func, n, fcall, quit_flag, accept, revert):
        """Replace a duplicate vector in Xnew by either reverting or generating a new random one."""
        m = n  # Index of vector in Xnew
        if revert:
            Xnew['vectors'][m] = X['vectors'][n]
            Xnew['values'][m] = X['values'][n]
            Xnew['vectors_and_derived'][m] = X['vectors_and_derived'][n]
            if run_params['jDE']:
                Xnew['FjDE'][m] = X['FjDE'][n]
                Xnew['CrjDE'][m] = X['CrjDE'][n]
                if run_params['lambdajDE']:
                    Xnew['lambdajDE'][m] = X['lambdajDE'][n]
        else:
            newvector = np.random.random(size=run_params['D']) * (run_params['upperbounds'] - run_params['lowerbounds']) + run_params['lowerbounds']
            Xnew['vectors'][m] = newvector
            Xnew['vectors_and_derived'][m] = self.roundvector(newvector, run_params)
            Xnew['values'][m] = func(Xnew['vectors_and_derived'][m], fcall, quit_flag, True, run_params['context'])
            if run_params['jDE']:
                Xnew['FjDE'][m] = self.init_FjDE(1)
                Xnew['CrjDE'][m] = self.init_CrjDE(1)
                if run_params['lambdajDE']:
                    Xnew['lambdajDE'][m] = self.init_lambdajDE(1)
        allvecs[n] = Xnew['vectors'][m]

    def roundvector(self, vector, run_params):
        """Placeholder for rounding function."""
        # Placeholder for whatever rounding logic is needed for derived vectors
        return vector

    def init_FjDE(self, size):
        """Initialize FjDE values."""
        return np.random.uniform(low=0.1, high=0.9, size=size)

    def init_CrjDE(self, size):
        """Initialize CrjDE values."""
        return np.random.uniform(low=0.1, high=0.9, size=size)

    def init_lambdajDE(self, size):
        """Initialize lambdajDE values."""
        return np.random.uniform(low=0.0, high=1.0, size=size)
