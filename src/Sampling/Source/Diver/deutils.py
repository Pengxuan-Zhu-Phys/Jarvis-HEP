#!/usr/bin/env python3 

import numpy as np

class DEUtils:
    @staticmethod
    def int_to_string(val):
        """
        Converts an integer to a string.
        """
        return str(val).strip()

    @staticmethod
    def quit_de(error_message=None):
        """
        Exits the program with an optional error message.
        """
        if error_message:
            print(error_message)
        exit()

    @staticmethod
    def quit_all_processes(error_message=None):
        """
        Exits the program for all processes (simulating MPI_ABORT behavior).
        """
        if error_message:
            print(error_message)
        exit()

    @staticmethod
    def roundvector(trialvector, run_params):
        """
        Rounds vectors to the nearest discrete values for all dimensions listed in run_params%discrete.
        All other dimensions are kept the same.
        """
        result = np.copy(trialvector)
        result[run_params['discrete']] = np.round(result[run_params['discrete']])
        return result

    @staticmethod
    def newBFs(X, BF):
        """
        Updates the current best-fit vector.
        """
        bestloc = np.argmin(X['values'])
        bestvalue = X['values'][bestloc]

        if bestvalue <= BF['values'][0]:
            BF['values'][0] = bestvalue
            BF['vectors'][0, :] = X['vectors'][bestloc, :]
            BF['vectors_and_derived'][0, :] = X['vectors_and_derived'][bestloc, :]

    @staticmethod
    def update_acceptance(accept, fcall, NP, verbose=False):
        """
        Calculates the acceptance rate and total number of function calls.
        Since we are not using MPI, this function will return the local values.
        """
        totaccept = np.sum(accept)
        totfcall = np.sum(fcall)

        if verbose:
            print(f"  Acceptance rate: {totaccept / float(NP)}")

        return totaccept, totfcall

    @staticmethod
    def sync(flag):
        """
        Syncs the value of a flag (logical OR for MPI).
        In this version, since we are not using MPI, it simply returns the flag.
        """
        return flag

# Example of how to structure parameters in Python (like run_params in Fortran)
run_params = {
    'D': 10,  # Number of dimensions
    'discrete': [0, 2, 4],  # Example discrete dimensions
}

# Example of a population structure (similar to X and BF)
X = {
    'values': np.random.rand(10),
    'vectors': np.random.rand(10, 5),
    'vectors_and_derived': np.random.rand(10, 7)
}

BF = {
    'values': np.array([float('inf')]),
    'vectors': np.zeros((1, 5)),
    'vectors_and_derived': np.zeros((1, 7))
}
