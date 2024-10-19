#!/usr/bin/env python3 

import os
import numpy as np
from typing import Optional

class IOModule:
    Ztolscale = 100.0
    Ftolscale = 100.0
    Bndtolscale = 100.0

    def __init__(self):
        self.rawlun = None
        self.samlun = None
        self.devolun = None
        self.rparamlun = None

    def io_begin(self, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF, path=None, prior=None, restart=False):
        """
        Initialize IO operations, including resuming from a previous run or starting a new run.
        """
        if restart:
            if not path:
                raise ValueError("Resuming a Diver run requires the path argument to be set to the location of the previous run files.")
            self.resume(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF, prior=prior)
        elif run_params['mpirank'] == 0 and not run_params['disableIO']:
            if not path:
                raise ValueError("The path argument must be set unless disableIO=True and not attempting to resume an old run.")
            # Create output .raw file
            if run_params['outputRaw']:
                print(f"Creating Diver .raw file at {path}.raw")
                open(f"{path}.raw", 'w').close()
            # Create output .sam file if necessary
            if run_params['outputSam'] and (run_params['D_derived'] != 0 or len(run_params['discrete']) != 0):
                print(f"Creating Diver .sam file at {path}.sam")
                open(f"{path}.sam", 'w').close()

    def save_all(self, X, BF, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, path=None, final=False):
        """
        Save all current state and samples.
        """
        if not final:
            Nsamples_saved += run_params['DE']['NP']
            self.save_samples(X, civ, gen, run_params, path)
        self.save_state(civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF, path)

    def save_samples(self, X, civ, gen, run_params, path=None):
        """
        Save samples to .raw and .sam files.
        """
        if run_params['disableIO']:
            return

        if run_params['outputRaw']:
            with open(f"{path}.raw", 'a') as raw_file:
                format_string_raw = f"(2E20.9, 2I6, {run_params['D']}E20.9)"
                for i in range(len(X['weights'])):
                    raw_file.write(f"{X['multiplicities'][i]:.9e} {X['values'][i]:.9e} {civ:6d} {gen:6d} {' '.join(map(str, X['vectors'][i, :]))}\n")

        if run_params['outputSam'] and (run_params['D_derived'] != 0 or len(run_params['discrete']) != 0):
            with open(f"{path}.sam", 'a') as sam_file:
                format_string_sam = f"(2E20.9, 2I6, {run_params['D'] + run_params['D_derived']}E20.9)"
                for i in range(len(X['weights'])):
                    sam_file.write(f"{X['multiplicities'][i]:.9e} {X['values'][i]:.9e} {civ:6d} {gen:6d} {' '.join(map(str, X['vectors_and_derived'][i, :]))}\n")

    def save_run_params(self, run_params, path=None):
        """
        Save the run parameters to a file.
        """
        if run_params['disableIO']:
            return

        with open(f"{path}.rparam", 'w') as rparam_file:
            rparam_file.write(f"{run_params['DE']['NP']:6d}\n")
            rparam_file.write(f"{run_params['DE']['jDE']}\n")
            rparam_file.write(f"{run_params['DE']['lambdajDE']}\n")
            rparam_file.write(f"{run_params['DE']['Fsize']:4d}\n")
            if run_params['DE']['Fsize'] != 0 and not run_params['DE']['jDE']:
                rparam_file.write(f"{run_params['DE']['F']}\n")
            rparam_file.write(f"{run_params['DE']['lambda']:.9e}\n")
            rparam_file.write(f"{run_params['DE']['current']}\n")
            rparam_file.write(f"{run_params['DE']['Cr']:.9e}\n")
            rparam_file.write(f"{run_params['DE']['expon']}\n")
            rparam_file.write(f"{run_params['DE']['bconstrain']:6d}\n")
            rparam_file.write(f"{run_params['D']:6d} {run_params['D_derived']:6d}\n")
            rparam_file.write(f"{run_params['lowerbounds']}\n")
            rparam_file.write(f"{run_params['upperbounds']}\n")
            rparam_file.write(f"{run_params['D_discrete']:6d}\n")
            if run_params['D_discrete'] != 0:
                rparam_file.write(f"{run_params['discrete']}\n")
                rparam_file.write(f"{run_params['partitionDiscrete']}\n")
                if run_params['partitionDiscrete']:
                    rparam_file.write(f"{run_params['repeat_scales']}\n")
                    rparam_file.write(f"{run_params['subpopNP']:6d}\n")
            rparam_file.write(f"{run_params['numciv']:6d} {run_params['numgen']:6d}\n")
            rparam_file.write(f"{run_params['convthresh']:.9e}\n")
            rparam_file.write(f"{run_params['convsteps']:6d}\n")
            rparam_file.write(f"{run_params['tol']:.9e}\n")
            rparam_file.write(f"{run_params['maxNodePop']:.9e}\n")
            rparam_file.write(f"{run_params['calcZ']}\n")
            rparam_file.write(f"{run_params['disableIO']}\n")
            rparam_file.write(f"{run_params['outputRaw']}\n")
            rparam_file.write(f"{run_params['outputSam']}\n")
            rparam_file.write(f"{run_params['savefreq']:6d}\n")
            rparam_file.write(f"{run_params['DE']['removeDuplicates']}\n")
            rparam_file.write(f"{run_params['verbose']:6d}\n")
            rparam_file.write(f"{run_params['convergence_criterion']:6d}\n")

    def save_state(self, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF, path=None):
        """
        Save the current state of the run to a file.
        """
        if run_params['disableIO']:
            return

        with open(f"{path}.devo", 'w') as devo_file:
            devo_file.write(f"{civ:10d} {gen:10d}\n")
            devo_file.write(f"{Z:.9e} {Zmsq:.9e} {Zerr:.9e} {Zold:.9e}\n")
            devo_file.write(f"{Nsamples:10d} {Nsamples_saved:10d} {fcall:10d}\n")
            devo_file.write(f"{BF['values'][0]:.9e}\n")
            devo_file.write(f"{BF['vectors'][0, :]}\n")
            devo_file.write(f"{BF['vectors_and_derived'][0, :]}\n")
            devo_file.write(f"{X['values']}\n")
            devo_file.write(f"{X['vectors']}\n")
            devo_file.write(f"{X['vectors_and_derived']}\n")


    def resume(self, path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF, prior=None):
        """
        Resume from a previous run state.
        """
        print("Restoring from previous run...")

        # Read the run state
        run_params_restored = self.read_state(path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF)

        # If restoring evidence calculation, validate consistency with prior setup
        if run_params["calcZ"]:
            if not os.path.exists(f"{path}.raw"):
                raise ValueError("Cannot resume in Bayesian mode if .raw file was not output.")
            if not run_params_restored["calcZ"]:
                raise ValueError("Cannot resume in Bayesian mode from a non-Bayesian run.")
            if prior is None:
                raise ValueError("Evidence calculation requested without specifying a prior.")

            # Check for consistency in prior bounds
            if np.any(np.abs(run_params["upperbounds"] - run_params_restored["upperbounds"]) / run_params["upperbounds"] >= self.Bndtolscale * np.finfo(float).eps) or \
               np.any(np.abs(run_params["lowerbounds"] - run_params_restored["lowerbounds"]) / run_params["lowerbounds"] >= self.Bndtolscale * np.finfo(float).eps):
                raise ValueError("Cannot resume in Bayesian mode with a modified prior box.")

            if (run_params["convthresh"] != run_params_restored["convthresh"] or 
                run_params["convsteps"] != run_params_restored["convsteps"]):
                print("WARNING: Changing the generation-level convergence parameters mid-run may make evidence inaccurate.")

            # Validate the population settings
            if run_params["DE"]["Fsize"] != run_params_restored["DE"]["Fsize"]:
                print("WARNING: Changing the number of F parameters mid-run may make evidence inaccurate.")
            elif run_params["DE"]["Fsize"] != 0:
                if np.any(np.abs(run_params["DE"]["F"] - run_params_restored["DE"]["F"]) / run_params["DE"]["F"] >= self.Ftolscale * np.finfo(float).eps):
                    print("WARNING: Changing F values mid-run may make evidence inaccurate.")

            # Validate DE algorithm consistency
            if (run_params["DE"]["lambda"] != run_params_restored["DE"]["lambda"] or 
                run_params["DE"]["current"] != run_params_restored["DE"]["current"] or 
                run_params["DE"]["Cr"] != run_params_restored["DE"]["Cr"] or 
                run_params["DE"]["expon"] != run_params_restored["DE"]["expon"] or 
                run_params["DE"]["bconstrain"] != run_params_restored["DE"]["bconstrain"] or 
                run_params["DE"]["jDE"] != run_params_restored["DE"]["jDE"] or 
                run_params["DE"]["lambdajDE"] != run_params_restored["DE"]["lambdajDE"]):
                print("WARNING: Changing DE algorithm mid-run may make evidence inaccurate!")

        print("Restored successfully.")

    def read_state(self, path, civ, gen, Z, Zmsq, Zerr, Zold, Nsamples, Nsamples_saved, fcall, run_params, X, BF):
        """
        Read the run state from a file.
        """
        # Check that .rparam exists
        if not os.path.exists(f"{path}.rparam"):
            raise ValueError(f"{path}.rparam does not exist. Cannot resume Diver.")

        # Read in the run parameters
        run_params_restored = {}
        with open(f"{path}.rparam", 'r') as rparam_file:
            run_params_restored["DE"] = {
                "NP": int(rparam_file.readline().strip()),
                "jDE": bool(int(rparam_file.readline().strip())),
                "lambdajDE": bool(int(rparam_file.readline().strip())),
                "Fsize": int(rparam_file.readline().strip())
            }
            if run_params_restored["DE"]["Fsize"] != 0 and not run_params_restored["DE"]["jDE"]:
                run_params_restored["DE"]["F"] = np.array(list(map(float, rparam_file.readline().strip().split())))

            run_params_restored["DE"]["lambda"] = float(rparam_file.readline().strip())
            run_params_restored["DE"]["current"] = bool(int(rparam_file.readline().strip()))
            run_params_restored["DE"]["Cr"] = float(rparam_file.readline().strip())
            run_params_restored["DE"]["expon"] = bool(int(rparam_file.readline().strip()))
            run_params_restored["DE"]["bconstrain"] = int(rparam_file.readline().strip())
            run_params_restored["D"] = int(rparam_file.readline().strip())
            run_params_restored["D_derived"] = int(rparam_file.readline().strip())
            run_params_restored["lowerbounds"] = np.array(list(map(float, rparam_file.readline().strip().split())))
            run_params_restored["upperbounds"] = np.array(list(map(float, rparam_file.readline().strip().split())))
            run_params_restored["D_discrete"] = int(rparam_file.readline().strip())
            if run_params_restored["D_discrete"] != 0:
                run_params_restored["discrete"] = np.array(list(map(int, rparam_file.readline().strip().split())))
                run_params_restored["partitionDiscrete"] = bool(int(rparam_file.readline().strip()))
                if run_params_restored["partitionDiscrete"]:
                    run_params_restored["repeat_scales"] = np.array(list(map(int, rparam_file.readline().strip().split())))
                    run_params_restored["subpopNP"] = int(rparam_file.readline().strip())
            else:
                run_params_restored["discrete"] = []

            run_params_restored["numciv"] = int(rparam_file.readline().strip())
            run_params_restored["numgen"] = int(rparam_file.readline().strip())
            run_params_restored["convthresh"] = float(rparam_file.readline().strip())
            run_params_restored["convsteps"] = int(rparam_file.readline().strip())
            run_params_restored["tol"] = float(rparam_file.readline().strip())
            run_params_restored["maxNodePop"] = float(rparam_file.readline().strip())
            run_params_restored["calcZ"] = bool(int(rparam_file.readline().strip()))
            run_params_restored["disableIO"] = bool(int(rparam_file.readline().strip()))
            run_params_restored["outputRaw"] = bool(int(rparam_file.readline().strip()))
            run_params_restored["outputSam"] = bool(int(rparam_file.readline().strip()))
            run_params_restored["savefreq"] = int(rparam_file.readline().strip())
            run_params_restored["DE"]["removeDuplicates"] = bool(int(rparam_file.readline().strip()))
            run_params_restored["verbose"] = int(rparam_file.readline().strip())
            run_params_restored["convergence_criterion"] = int(rparam_file.readline().strip())

        # Check that .devo exists
        if not os.path.exists(f"{path}.devo"):
            raise ValueError(f"{path}.devo does not exist. Cannot resume Diver.")

        with open(f"{path}.devo", 'r') as devo_file:
            civ, gen = map(int, devo_file.readline().split())
            Z, Zmsq, Zerr, Zold = map(float, devo_file.readline().split())
            Nsamples, Nsamples_saved, fcall = map(int, devo_file.readline().split())
            BF["values"][0] = float(devo_file.readline().strip())
            BF["vectors"] = np.array(list(map(float, devo_file.readline().split()))).reshape((1, run_params_restored["D"]))
            BF["vectors_and_derived"] = np.array(list(map(float, devo_file.readline().split()))).reshape((1, run_params_restored["D"] + run_params_restored["D_derived"]))

        return run_params_restored
