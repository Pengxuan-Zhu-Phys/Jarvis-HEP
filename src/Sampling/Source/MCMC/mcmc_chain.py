#!/usr/bin/env python3 

import numpy as np 


class MCMCChain:
    def __init__(self, initial_param, proposal_scale, n_iterations):
        self.param = initial_param             # Current parameter (state)
        # print("self.param -> ", self.param)
        self.proposal_scale = proposal_scale   # Standard deviation for proposal steps
        self.n_iterations = n_iterations       # Total iterations for this chain
        self.iterations = 0                    # How many iterations have been performed
        self.last_loglikelihood = None          # Latest computed log likelihood
        self.proposed_param = None 

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        if self.iterations == 0:
            # print("Using random to sampling first sample")
            proposed_param = np.random.rand(self.param.shape[0])
            self.proposed_param = proposed_param
            return proposed_param
        while self.iterations:
            # print("using step to generate the next sample")
            # Propose a new parameter by adding a Gaussian step to the current parameter.
            step = np.random.normal(0, self.proposal_scale, size=self.param.shape)
            proposed_param = self.param + step
            if np.all((proposed_param >= 0) & (proposed_param <= 1)): 
                self.proposed_param = proposed_param
                return proposed_param


    def update(self, new_loglikelihood):
        # Update the chain state with the new parameter and its log likelihood.
        if self.iterations > 0:
            from numpy import exp 
            likelihood2 = exp(new_loglikelihood)
            likelihood1 = exp(self.last_loglikelihood)
            accecptance = likelihood2 / likelihood1
            if np.random.rand() < accecptance:
                self.param = self.proposed_param
                self.last_loglikelihood = new_loglikelihood
        else: 
            self.param = self.proposed_param
            self.last_loglikelihood = new_loglikelihood
        self.iterations += 1
        
        
        

