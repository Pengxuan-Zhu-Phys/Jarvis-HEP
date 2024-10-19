#!/usr/bin/env python3 

import numpy as np


class MutationModule:

    tauF = 0.1
    Fl = 0.1
    Fu = 0.9
    taulambda = 0.1

    def get_subpopulation(self, X, run_params, n):
        """
        Determine the subpopulation of points that have the same partitioned discrete parameters.
        """
        discrete_vals = np.rint(X['vectors'][n, run_params['discrete']])
        subpopulation_indices = []

        for i in range(run_params['DE']['NP']):
            same_subpop = all(np.rint(X['vectors'][i, run_params['discrete']]) == discrete_vals)
            if same_subpop:
                subpopulation_indices.append(i)

        if len(subpopulation_indices) != run_params['subpopNP']:
            raise ValueError('ERROR: subpopulation size does not match run_params["subpopNP"]!')

        Xsub = {key: X[key][subpopulation_indices] for key in X}
        nsub = subpopulation_indices.index(n)

        return Xsub, nsub

    def mutate(self, X, n, run_params):
        """
        Apply the appropriate mutation strategy depending on DE parameters.
        """
        if run_params['DE']['lambdajDE']:
            trialF = self.new_F(X, n)
            triallambda = self.new_lambda(X, n)
            V = self.lambdajDEmutation(X, n, run_params, trialF, triallambda)
        elif run_params['DE']['jDE']:
            trialF = self.new_F(X, n)
            triallambda = 0.0
            V = self.jDEmutation(X, n, run_params, trialF)
        else:
            trialF = 0.0
            triallambda = 0.0
            V = self.genmutation(X, n, run_params)

        return V, trialF, triallambda

    def genmutation(self, X, n, run_params):
        """
        General mutation strategy: V = lambda*X_best + (1-lambda)*X_I + Sum_q F(q)*(X_J(q) - X_K(q))
        """
        Fsize = run_params['DE']['Fsize']
        D = run_params['D']
        subpopNP = run_params['subpopNP']

        rbest = np.argmin(X['values']) if run_params['DE']['lambda'] > 0 else n
        ri = n if run_params['DE']['current'] or run_params['DE']['lambda'] == 1.0 else self.pick_random_ri(subpopNP, n, rbest)

        r = self.pick_unique_r_vectors(subpopNP, n, ri, rbest, Fsize)

        sumF = np.zeros(D)
        for q in range(Fsize):
            sumF += run_params['DE']['F'][q] * (X['vectors'][r[2*q-1], :] - X['vectors'][r[2*q], :])

        V = run_params['DE']['lambda'] * X['vectors'][rbest, :] + (1 - run_params['DE']['lambda']) * X['vectors'][ri, :] + sumF

        return V

    def jDEmutation(self, X, n, run_params, trialF):
        """
        jDE mutation using self-adaptive F parameter.
        """
        subpopNP = run_params['subpopNP']
        D = run_params['D']
        rbest = np.argmin(X['values']) if run_params['DE']['lambda'] > 0 else n

        r1, r2, r3 = self.pick_three_random_vectors(subpopNP, n, rbest)

        V = run_params['DE']['lambda'] * X['vectors'][rbest, :] + \
            (1 - run_params['DE']['lambda']) * X['vectors'][r1, :] + \
            trialF * (X['vectors'][r2, :] - X['vectors'][r3, :])

        return V

    def lambdajDEmutation(self, X, n, run_params, trialF, triallambda):
        """
        lambdajDE mutation using self-adaptive F and lambda parameters.
        """
        subpopNP = run_params['subpopNP']
        D = run_params['D']
        rbest = np.argmin(X['values'])

        r1, r2, r3 = self.pick_three_random_vectors(subpopNP, n, rbest)

        V = triallambda * X['vectors'][rbest, :] + \
            (1 - triallambda) * X['vectors'][r1, :] + \
            trialF * (X['vectors'][r2, :] - X['vectors'][r3, :])

        return V

    def new_F(self, X, n):
        """
        Choose a trial value for the F parameter.
        """
        rand1 = np.random.rand()
        if rand1 < self.tauF:
            rand2 = np.random.rand()
            newF = self.Fl + rand2 * self.Fu  # Select a new value between 0.1 and 1.0
        else:
            newF = X['FjDE'][n]  # Inherit F from previous generation
        return newF

    def new_lambda(self, X, n):
        """
        Choose a trial value for the lambda parameter.
        """
        rand1 = np.random.rand()
        if rand1 < self.taulambda:
            newlambda = np.random.rand()  # Select a new value between 0.0 and 1.0
        else:
            newlambda = X['lambdajDE'][n]  # Inherit lambda from previous generation
        return newlambda

    def init_FjDE(self, size):
        """
        Initialize FjDE values.
        """
        rand = np.random.rand(size)
        return self.Fl + rand * self.Fu

    def init_lambdajDE(self, size):
        """
        Initialize lambdajDE values.
        """
        return np.random.rand(size)

    def random_int(self, min_val, max_val):
        """
        Choose a random integer between min and max, inclusive.
        """
        return np.random.randint(min_val, max_val + 1)

    def pick_random_ri(self, subpopNP, n, rbest):
        """
        Pick a random index ri not equal to n or rbest.
        """
        while True:
            ri = self.random_int(1, subpopNP)
            if ri != n and ri != rbest:
                break
        return ri

    def pick_unique_r_vectors(self, subpopNP, n, ri, rbest, Fsize):
        """
        Pick unique r(q) vectors from the population.
        """
        r = []
        for _ in range(2 * Fsize):
            while True:
                rand_val = self.random_int(1, subpopNP)
                if rand_val not in [n, ri, rbest] + r:
                    r.append(rand_val)
                    break
        return r

    def pick_three_random_vectors(self, subpopNP, n, rbest):
        """
        Pick three distinct random vectors r1, r2, r3 from the population.
        """
        while True:
            r1 = self.random_int(1, subpopNP)
            if r1 not in [n, rbest]:
                break
        while True:
            r2 = self.random_int(1, subpopNP)
            if r2 not in [n, r1, rbest]:
                break
        while True:
            r3 = self.random_int(1, subpopNP)
            if r3 not in [n, r1, r2, rbest]:
                break
        return r1, r2, r3
