#!/usr/bin/env python3 

import numpy as np
import scipy.special
from scipy.stats import binom, poisson, beta, expon, gamma
import scipy

class Variable:
    def __init__(self, name, description, distribution, parameters):
        self._name = name
        self._description = description
        self._distribution = distribution
        self._parameters = parameters

    @property
    def name(self):
        return self._name

    @property
    def description(self):
        return self._description

    @property
    def distribution(self):
        return self._distribution

    @property
    def parameters(self):
        return self._parameters

    def generate(self):
        """
        Generate single sample according to the defination
        """
        if self._distribution == 'Flat':
            return np.random.uniform(self._parameters['min'], self._parameters['max'])
        elif self._distribution == 'Log':
            return np.exp(np.random.uniform(np.log(self._parameters['min']), np.log(self._parameters['max'])))
        elif self._distribution == 'Normal':
            return np.random.normal(self._parameters['mean'], self._parameters['stddev'])
        elif self._distribution == 'Log-Normal':
            return np.random.lognormal(self._parameters['mean'], self._parameters['stddev'])
        elif self._distribution == "Logit":
            # Generate a sample from the Logit distribution
            mu = self._parameters['location']
            s = self._parameters['scale']
            p = np.random.uniform(0, 1)  # Generate a uniform random number between 0 and 1
            return mu + s * np.log(p / (1 - p))  # Apply the inverse CDF of the Logit distribution
        elif self._distribution == 'Binomial':
            return binom.rvs(self._parameters['n'], self._parameters['p'])
        elif self._distribution == 'Poisson':
            return poisson.rvs(self._parameters['lambda'])
        elif self._distribution == 'Beta':
            return beta.rvs(self._parameters['alpha'], self._parameters['beta'])
        elif self._distribution == 'Exponential':
            # The exponential distibution in scipy is the inverse of the probility 
            return expon.rvs(scale=1/self._parameters['rate'])  
        elif self._distribution == 'Gamma':
            return gamma.rvs(self._parameters['shape'], scale=self._parameters['scale'])
        else:
            raise ValueError(f"Unsupported distribution type: {self._distribution}")

    def map_standard_random_to_distribution(self, std_rand):
        """
        mapping the standard random number in range [0, 1] into the Variable distribution.
        """
        if self.distribution == 'Flat':
            min_val = self.parameters['min']
            max_val = self.parameters['max']
            return min_val + (max_val - min_val) * std_rand
        elif self.distribution == 'Log':
            log_min = np.log(self.parameters['min'])
            log_max = np.log(self.parameters['max'])
            return np.exp(log_min + (log_max - log_min) * std_rand)
        elif self.distribution == 'Normal':
            mean = self.parameters['mean']
            stddev = self.parameters['stddev']
            return scipy.stats.norm.ppf(std_rand, loc=mean, scale=stddev)
        elif self.distribution == 'Log-Normal':
            mean = self.parameters['mean']
            stddev = self.parameters['stddev']
            return scipy.stats.lognorm.ppf(std_rand, s=stddev, scale=np.exp(mean))
        elif self.distribution == 'Logit':
            # Using the Percent Point Function (PPF) of the logistic distribution, 
            # which is the inverse function of the logit, to map from the [0, 1] range to the real number range.
            # Note: There are no direct distribution parameters here, because the Logit distribution is usually 
            # used to describe probabilities and does not involve direct distribution parameters like mean or standard deviation.
            location = self.parameters.get('location', 0)  # The default value is 0
            scale = self.parameters.get('scale', 1)  # The default value is 1
            # The scipy's logit function is used as the inverse function of the logistic CDF
            return scipy.special.logit(std_rand) * scale + location
        # Notice: the following distribution not support the mapping function 
        elif self.distribution in ['Binomial', 'Poisson', 'Beta', 'Exponential', 'Gamma']:
            raise ValueError(f"Distribution type '{self.distribution}' does not support mapping from standard random numbers.")
        else:
            raise ValueError(f"Unsupported distribution type: {self.distribution}")