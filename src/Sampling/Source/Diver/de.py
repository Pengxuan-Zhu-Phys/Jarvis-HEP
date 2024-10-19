#!/usr/bin/env python3 

# This version of Diver is rewrite by Pengxuan Zhu, just transform the original Fortran code into Python
# Diver Github webpage: https://github.com/diveropt/Diver/tree/master 
# Any papers that use results or insights obtained with Diver should cite this paper:
# 1. Martinez, McKay, Farmer, Scott, Roebber, Putze & Conrad 2017, 
#       European Physical Journal C 77 (2017) 761, arXiv:1705.07959. 
# One can also find detailed performance comparisons of Diver with other samplers and optimisers in both the above paper and in:
# 2. DarkMachines High Dimensional Sampling Group, 
#       JHEP 05 (2021) 108, arXiv:2101.04525. 
#      In particular, this latter paper demonstrates that Diver significantly outperforms Scipy's implementation of differential evolution.

import numpy as np
import time
import uuid

class DiverVirtual:
    def __init__(self, D, NP, lowerbounds, upperbounds, max_gen=1000):
        """
        Initializes the DiverVirtual class with the given parameters.

        Parameters:
            D : int
                Dimensionality of the problem.
            NP : int
                Population size.
            lowerbounds : list or ndarray
                Lower bounds for each dimension.
            upperbounds : list or ndarray
                Upper bounds for each dimension.
            max_gen : int, optional
                Maximum number of generations (default is 1000).
        """
        # Problem parameters
        self.D = D  # Dimensionality
        self.NP = NP  # Population size
        self.lowerbounds = np.array(lowerbounds)
        self.upperbounds = np.array(upperbounds)
        self.max_gen = max_gen  # Maximum generations

        # DE parameters
        self.bconstrain = 1  # Boundary constraint method
        self.removeDuplicates = True

        # Initialize populations
        self.X = None  # Current population
        self.Xnew = None  # New population

        # Function call counters and flags
        self.fcall = 0
        self.accept = 0
        self.quit_flag = False
        self._debug = False

        # Initialize the populations
        self._init_population()

    def _init_population(self):
        """
        Initializes the populations with random vectors and evaluates them.
        """
        # Initialize the old and new populations as dictionaries
        self.X = self._create_population()
        self.Xnew = self._create_population()

        # Randomly initialize the vectors within the specified bounds
        self.X['vectors'] = np.random.uniform(
            self.lowerbounds, self.upperbounds, (self.NP, self.D)
        )
        # Evaluate the initial population using the objective function
        self.X['values'] = np.apply_along_axis(
            self.objective_function, 1, self.X['vectors']
        )
        self.fcall += self.NP

        # Assign UUIDs to each individual
        self.X['uuids'] = np.array([uuid.uuid4() for _ in range(self.NP)])

    def _create_population(self):
        """
        Creates a population data structure.

        Returns:
            population : dict
                A dictionary representing the population.
        """
        population = {
            'vectors': np.zeros((self.NP, self.D)),
            'values': np.full(self.NP, np.inf),
            'uuids': np.array([uuid.uuid4() for _ in range(self.NP)])
        }
        # print(population)
        # time.sleep(3)
        return population

    def mutate(self, n):
        """
        Performs mutation to generate a donor vector.

        Parameters:
            n : int
                Index of the target vector in the population.

        Returns:
            V : ndarray
                The donor vector.
        """
        # DE/rand/1 mutation strategy
        indices = np.delete(np.arange(self.NP), n)
        a, b, c = np.random.choice(indices, 3, replace=False)
        F = 0.5  # Scaling factor (can be customized)
        V = self.X['vectors'][a] + F * (self.X['vectors'][b] - self.X['vectors'][c])
        return V

    def crossover(self, V, n):
        """
        Performs crossover to generate a trial vector.

        Parameters:
            V : ndarray
                The donor vector.
            n : int
                Index of the target vector in the population.

        Returns:
            U : ndarray
                The trial vector.
        """
        U = np.copy(self.X['vectors'][n])
        CR = 0.9  # Crossover probability (can be customized)
        rand_vals = np.random.rand(self.D)
        jrand = np.random.randint(self.D)
        crossover_indices = rand_vals < CR
        crossover_indices[jrand] = True  # Ensure at least one parameter is from V
        U[crossover_indices] = V[crossover_indices]
        return U

    def selection(self, Us, trial_values):
        """
        Performs selection between the target and trial vectors.

        Parameters:
            Us : ndarray
                The array of trial vectors.
            trial_values : ndarray
                The array of trial values.

        Returns:
            None
        """
        improved = trial_values <= self.X['values']
        if self._debug:
            print("Line 162 -> ", improved.shape)
            print("Line 163 -> ", self.Xnew['uuids'][improved][0], self.X['uuids'][improved][0])
            print("Line 164 -> ", self.Xnew['values'][improved][0], self.X['values'][improved][0])
            print("Line 165 -> ", self.Xnew['vectors'][improved][0], self.X['vectors'][improved][0])

            print("Line 167 -> ", self.Xnew['uuids'][~improved][0], self.X['uuids'][~improved][0])
            print("Line 168 -> ", self.Xnew['values'][~improved][0], self.X['values'][~improved][0])
            print("Line 169 -> ", self.Xnew['vectors'][~improved][0], self.X['vectors'][~improved][0])
            time.sleep(10)

        self.Xnew['vectors'][improved] = Us[improved]
        self.Xnew['values'][improved] = trial_values[improved]
        self.Xnew['uuids'][improved] = self.X['uuids'][improved]
        # For non-improved individuals
        
        not_improved = ~improved
        self.Xnew['vectors'][not_improved] = self.X['vectors'][not_improved]
        self.Xnew['values'][not_improved] = self.X['values'][not_improved]
        self.Xnew['uuids'][not_improved] = self.X['uuids'][not_improved]
        self.accept += np.sum(improved)

    def _apply_boundary_constraints(self, Us):
        """
        Applies boundary constraints to an array of trial vectors Us.

        Parameters:
            Us : ndarray
                The array of trial vectors.

        Returns:
            Us : ndarray
                The array of trial vectors after applying boundary constraints.
        """
        if self.bconstrain == 1:
            # 'Brick wall' constraint: clip the values
            Us = np.clip(Us, self.lowerbounds, self.upperbounds)
        elif self.bconstrain == 2:
            # Random re-initialization
            mask = (Us < self.lowerbounds) | (Us > self.upperbounds)
            Us[mask] = np.random.uniform(self.lowerbounds, self.upperbounds, size=Us.shape)[mask]
        elif self.bconstrain == 3:
            # Reflection
            Us = np.where(Us < self.lowerbounds, self.lowerbounds + (self.lowerbounds - Us), Us)
            Us = np.where(Us > self.upperbounds, self.upperbounds - (Us - self.upperbounds), Us)
            Us = np.clip(Us, self.lowerbounds, self.upperbounds)
        else:
            # No boundary constraints enforced
            pass
        return Us

    def replace_generation(self):
        """
        Replaces the old generation with the new one after the population loop.

        Returns:
            None
        """
        # Remove duplicate vectors if desired
        if self.removeDuplicates:
            self._remove_duplicate_vectors()

        # Replace old population members with those calculated in Xnew
        self.X['vectors'] = self.Xnew['vectors'].copy()
        self.X['values'] = self.Xnew['values'].copy()
        self.X['uuids'] = self.Xnew['uuids'].copy()

    def _remove_duplicate_vectors(self):
        """
        Removes duplicate vectors from the population to maintain diversity.

        Returns:
            None
        """
        # Identify unique vectors
        unique_vectors, unique_indices = np.unique(self.Xnew['vectors'], axis=0, return_index=True)
        if len(unique_indices) < self.NP:
            # Duplicates exist
            # Identify indices of duplicates
            all_indices = np.arange(self.NP)
            duplicate_indices = np.setdiff1d(all_indices, unique_indices)
            # Replace duplicates with new random vectors
            for idx in duplicate_indices:
                self._replace_vector(idx)

    def _replace_vector(self, idx):
        """
        Replaces a duplicate vector with a new random vector.

        Parameters:
            idx : int
                Index of the vector to replace.

        Returns:
            None
        """
        # Generate a new random vector within bounds
        new_vector = np.random.uniform(self.lowerbounds, self.upperbounds)
        # Assign the new vector
        self.Xnew['vectors'][idx] = new_vector
        # Evaluate the new vector
        self.Xnew['values'][idx] = self.objective_function(new_vector)
        self.fcall += 1
        # Assign a new UUID
        self.Xnew['uuids'][idx] = uuid.uuid4()

    def run(self):
        """
        Runs the differential evolution algorithm.

        Returns:
            None
        """
        for gen in range(self.max_gen):
            self.accept = 0
            # Mutation and Crossover
            Vs = np.zeros((self.NP, self.D))
            Us = np.zeros((self.NP, self.D))
            for n in range(self.NP):
                V = self.mutate(n)
                U = self.crossover(V, n)
                Vs[n] = V
                Us[n] = U
            # Apply boundary constraints
            Us = self._apply_boundary_constraints(Us)
            # Evaluate all trial vectors
            trial_values = np.apply_along_axis(self.objective_function, 1, Us)
            self.fcall += self.NP
            # Selection
            self.selection(Us, trial_values)
            # Replace generation
            self.replace_generation()
            if self.quit_flag:
                break
            # Optional: Implement convergence criteria or logging

    def objective_function(self, x):
        """
        Evaluates the objective function at a given point x.

        Parameters:
            x : ndarray
                The point at which to evaluate the objective function.

        Returns:
            value : float
                The objective function value.
        """
        # Placeholder objective function (to be overridden in subclasses)
        raise NotImplementedError("Objective function must be defined in the subclass.")


class DiverjDE(DiverVirtual):
    def __init__(self, D, NP, lowerbounds, upperbounds, max_gen=1000):
        super().__init__(D, NP, lowerbounds, upperbounds, max_gen)
        # Set jDE parameters
        self.jDE = True
        self.FL = 0.1  # Lower bound for F
        self.FU = 0.9  # Upper bound for F
        self.tau1 = 0.1  # Parameter for adapting F
        self.tau2 = 0.1  # Parameter for adapting Cr

        # Initialize DE parameters for jDE
        self.X['FjDE'] = self._init_FjDE(self.NP)
        self.X['CrjDE'] = self._init_CrjDE(self.NP)

    def _init_FjDE(self, size):
        return self.FL + np.random.rand(size) * (self.FU - self.FL)

    def _init_CrjDE(self, size):
        return np.random.rand(size)

    def mutate(self, n):
        # Implement mutation strategy using self-adaptive F
        indices = [idx for idx in range(self.NP) if idx != n]
        a, b, c = np.random.choice(indices, 3, replace=False)
        F = self.X['FjDE'][n]
        V = self.X['vectors'][a] + F * (self.X['vectors'][b] - self.X['vectors'][c])
        return V

    def crossover(self, V, n):
        # Implement crossover using self-adaptive Cr
        U = np.zeros(self.D)
        Cr = self.X['CrjDE'][n]
        jrand = np.random.randint(self.D)
        for j in range(self.D):
            if np.random.rand() < Cr or j == jrand:
                U[j] = V[j]
            else:
                U[j] = self.X['vectors'][n][j]
        return U

    def selector(self, U, m, n):
        super().selector(U, m, n)
        # Update F and Cr for next generation
        if np.random.rand() < self.tau1:
            self.X['FjDE'][n] = self.FL + np.random.rand() * (self.FU - self.FL)
        if np.random.rand() < self.tau2:
            self.X['CrjDE'][n] = np.random.rand()


class DiverLambdajDE(DiverjDE):
    def __init__(self, D, NP, lowerbounds, upperbounds, max_gen=1000):
        super().__init__(D, NP, lowerbounds, upperbounds, max_gen)
        # Set lambdajDE parameters
        self.lambdajDE = True
        self.tau3 = 0.1  # Parameter for adapting lambda

        # Initialize lambda for each individual
        self.X['lambdajDE'] = self._init_lambdajDE(self.NP)
        self.Xnew['lambdajDE'] = np.zeros(self.NP)

    def _init_lambdajDE(self, size):
        return np.random.rand(size)

    def mutate(self, n):
        # Implement mutation strategy using self-adaptive F and lambda
        indices = [idx for idx in range(self.NP) if idx != n]
        a, b, c = np.random.choice(indices, 3, replace=False)
        F = self.X['FjDE'][n]
        lambda_param = self.X['lambdajDE'][n]
        V = self.X['vectors'][a] + F * (self.X['vectors'][b] - self.X['vectors'][c]) * lambda_param
        return V

    def selector(self, U, m, n):
        """
        Performs selection between the target and trial vectors and updates F, Cr, and lambda.
        """
        # Apply boundary constraints
        U = self._apply_boundary_constraints(U)

        # Evaluate the trial vector using the objective function
        trial_value = self.objective_function(U)
        self.fcall += 1

        # Selection process
        if trial_value <= self.X['values'][n]:
            # Accept the trial vector
            self.Xnew['vectors'][m] = U
            self.Xnew['vectors_and_derived'][m] = U
            self.Xnew['values'][m] = trial_value
            self.accept += 1
            # Keep the current F, Cr, and lambda
            self.Xnew['FjDE'][m] = self.X['FjDE'][n]
            self.Xnew['CrjDE'][m] = self.X['CrjDE'][n]
            self.Xnew['lambdajDE'][m] = self.X['lambdajDE'][n]
        else:
            # Keep the target vector
            self.Xnew['vectors'][m] = self.X['vectors'][n]
            self.Xnew['vectors_and_derived'][m] = self.X['vectors_and_derived'][n]
            self.Xnew['values'][m] = self.X['values'][n]
            self.Xnew['FjDE'][m] = self.X['FjDE'][n]
            self.Xnew['CrjDE'][m] = self.X['CrjDE'][n]
            self.Xnew['lambdajDE'][m] = self.X['lambdajDE'][n]

        # Update F, Cr, and lambda for the next generation
        if np.random.rand() < self.tau1:
            self.Xnew['FjDE'][m] = self.FL + np.random.rand() * (self.FU - self.FL)
        if np.random.rand() < self.tau2:
            self.Xnew['CrjDE'][m] = np.random.rand()
        if np.random.rand() < self.tau3:
            self.Xnew['lambdajDE'][m] = np.random.rand()

    def replace_generation(self):
        """
        Replaces the old generation with the new one after the population loop.
        """
        # Remove duplicate vectors if desired
        if self.removeDuplicates:
            self._remove_duplicate_vectors()

        # Replace old population members with those calculated in Xnew
        self.X['vectors'] = self.Xnew['vectors'].copy()
        self.X['vectors_and_derived'] = self.Xnew['vectors_and_derived'].copy()
        self.X['values'] = self.Xnew['values'].copy()
        self.X['FjDE'] = self.Xnew['FjDE'].copy()
        self.X['CrjDE'] = self.Xnew['CrjDE'].copy()
        self.X['lambdajDE'] = self.Xnew['lambdajDE'].copy()


import matplotlib.pyplot as plt

# DiverVirtual class as defined previously (make sure it's included in your code)

# Subclass with objective function and data collection
class MyDiverAlgorithm(DiverVirtual):
    def __init__(self, D, NP, lowerbounds, upperbounds, max_gen=1000):
        super().__init__(D, NP, lowerbounds, upperbounds, max_gen)
        # Initialize data collection lists
        self.best_values = []
        self.mean_values = []
        self.generations = []
        self.population_history = []  # To store population data at each generation
        self.values_history = []      # To store objective values at each generation

    def objective_function(self, x):
        # Ensure that the dimensionality is at least 2
        if len(x) < 2:
            raise ValueError("The objective function requires at least two dimensions.")

        term1 = - np.cos(x[0] * 3 * np.pi) * np.cos(x[1] * 3 * np.pi)
        term2 = 1 + 0.05 * x[0]**2 + 0.05 * x[1]**2
        return term1 + term2

    def run(self):
        """
        Runs the differential evolution algorithm and collects data for visualization.
        """
        for gen in range(self.max_gen):
            
            self.accept = 0
            # Mutation and Crossover
            Vs = np.zeros((self.NP, self.D))
            Us = np.zeros((self.NP, self.D))
            for n in range(self.NP):
                V = self.mutate(n)
                U = self.crossover(V, n)
                Vs[n] = V
                Us[n] = U
            # Apply boundary constraints
            Us = self._apply_boundary_constraints(Us)
            # Evaluate all trial vectors
            trial_values = np.apply_along_axis(self.objective_function, 1, Us)
            self.fcall += self.NP
            # Selection
            self.selection(Us, trial_values)
            # Replace generation
            self.replace_generation()
            if self.quit_flag:
                break
            # Data collection
            best_value = np.min(self.X['values'])
            mean_value = np.mean(self.X['values'])
            self.best_values.append(best_value)
            self.mean_values.append(mean_value)
            self.generations.append(gen)
            # Store population data for animation
            self.population_history.append(self.X['vectors'].copy())
            self.values_history.append(self.X['values'].copy())
            # Optional: Print progress
            # if gen % 10 == 0 or gen == self.max_gen - 1:
            print(f"Generation {gen}: Best Value = {best_value:.6f}")


def visualize_optimization(algo):
    """
    Visualizes the optimization process.

    Parameters:
        algo : MyDiverAlgorithm
            The instance of the algorithm after running optimization.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(algo.generations, algo.best_values, label='Best Objective Value')
    plt.plot(algo.generations, algo.mean_values, label='Mean Objective Value')
    plt.xlabel('Generation')
    plt.ylabel('Objective Value')
    plt.title('Optimization Progress')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def visualize_population(algo):
    """
    Visualizes the movement of the population in parameter space for 2D problems,
    with color coding based on the objective function values.

    Parameters:
        algo : MyDiverAlgorithm
            The instance of the algorithm after running optimization.
    """
    if algo.D != 2:
        print("Population visualization is only available for 2D problems.")
        return

    plt.figure(figsize=(8, 8))
    # Plot final population with color coding
    scatter = plt.scatter(
        algo.X['vectors'][:, 0],
        algo.X['vectors'][:, 1],
        c=algo.X['values'],
        cmap='viridis',
        label='Final Population',
        vmax=1.2,
        vmin=0
    )
    plt.colorbar(scatter, label='Objective Function Value')
    # Plot the best solution
    best_index = np.argmin(algo.X['values'])
    best_vector = algo.X['vectors'][best_index]
    plt.scatter(
        best_vector[0],
        best_vector[1],
        c='red',
        marker='*',
        s=200,
        edgecolors='black',
        label='Best Solution'
    )
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.title('Population in Parameter Space (Color-coded by Objective Value)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def animate_population(algo):
    """
    Creates an animation of the population over generations, retaining historical points.

    Parameters:
        algo : MyDiverAlgorithm
            The instance of the algorithm after running optimization.
    """
    if algo.D != 2:
        print("Animation is only available for 2D problems.")
        return

    fig, ax = plt.subplots(figsize=(8, 8))

    # Set up the plot limits
    x_min, x_max = algo.lowerbounds[0], algo.upperbounds[0]
    y_min, y_max = algo.lowerbounds[1], algo.upperbounds[1]
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Labels and title
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_title('Population Evolution Over Generations')

    # Initialize scatter plot
    scatter = ax.scatter([], [], c=[], cmap='viridis', s=4, vmin=min(algo.best_values), vmax=max(algo.mean_values))
    colorbar = plt.colorbar(scatter, ax=ax, label='Objective Function Value')

    # Best solution marker
    best_marker, = ax.plot([], [], 'r*', markersize=15, markeredgecolor='black', label='Best Solution')

    # Generation text
    generation_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

    # To store historical points
    historical_population = np.empty((0, 2))  # To store population points
    historical_values = np.empty(0)           # To store values of those points

    def init():
        scatter.set_offsets(np.empty((0, 2)))  # Ensure correct shape
        scatter.set_array(np.array([]))
        best_marker.set_data([], [])
        generation_text.set_text('')
        return scatter, best_marker, generation_text

    def update(frame):
        nonlocal historical_population, historical_values

        # Get population and values at the current generation
        population = algo.population_history[frame]
        values = algo.values_history[frame]

        # Append current generation points to historical data
        historical_population = np.vstack([historical_population, population])
        historical_values = np.concatenate([historical_values, values])

        # Update scatter plot with all points up to the current frame
        scatter.set_offsets(historical_population)
        scatter.set_array(historical_values)

        # Update best solution
        best_index = np.argmin(values)  # Finding the best in the current frame
        best_vector = population[best_index]
        best_marker.set_data([best_vector[0]], [best_vector[1]])  # Corrected line

        # Update generation text
        generation_text.set_text(f'Generation: {frame}')

        return scatter, best_marker, generation_text

    ani = FuncAnimation(
        fig, update, frames=len(algo.population_history),
        init_func=init, blit=False, interval=200, repeat=False
    )

    plt.legend()
    plt.show()


def animate_population_with_contour(algo):
    """
    Creates an animation of the population over generations,
    with the objective function contour as background.
    """
    if algo.D != 2:
        print("Animation is only available for 2D problems.")
        return

    fig, ax = plt.subplots(figsize=(8, 8))

    # Generate a grid of points
    x_min, x_max = algo.lowerbounds[0], algo.upperbounds[0]
    y_min, y_max = algo.lowerbounds[1], algo.upperbounds[1]
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, 200), np.linspace(y_min, y_max, 200))

    # Compute objective function values on the grid
    zz = np.array([
        algo.objective_function(np.array([x, y]))
        for x, y in zip(np.ravel(xx), np.ravel(yy))
    ])
    zz = zz.reshape(xx.shape)

    # Plot the contour map
    contour = ax.contourf(xx, yy, zz, levels=50, cmap='viridis')
    colorbar = plt.colorbar(contour, ax=ax, label='Objective Function Value')

    # Labels and title
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_title('Population Evolution Over Generations')

    # Initialize scatter plot
    scatter = ax.scatter([], [], c='white', edgecolors='black', label='Population')

    # Best solution marker
    best_marker, = ax.plot([], [], 'r*', markersize=15, markeredgecolor='black', label='Best Solution')

    # Generation text
    generation_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='white')

    def init():
        scatter.set_offsets(np.empty((0, 2)))  # Updated line
        best_marker.set_data([], [])
        generation_text.set_text('')
        return scatter, best_marker, generation_text

    def update(frame):
        # Get population at the current generation
        population = algo.population_history[frame]
        # Update scatter plot
        scatter.set_offsets(population)
        # Update best solution
        values = algo.values_history[frame]
        best_index = np.argmin(values)
        best_vector = population[best_index]
        best_marker.set_data(best_vector[0], best_vector[1])
        # Update generation text
        generation_text.set_text(f'Generation: {frame}')
        return scatter, best_marker, generation_text

    ani = FuncAnimation(
        fig, update, frames=len(algo.population_history),
        init_func=init, blit=False, interval=200, repeat=False  # Updated line
    )

    plt.legend()
    plt.show()



# Define problem parameters
D = 2  # Dimensionality (use 2D for animation)
NP = 50  # Population size
lowerbounds = [-2] * D
upperbounds = [2] * D
max_gen = 50  # Number of generations

# Create an instance of the algorithm
my_algo = MyDiverAlgorithm(D, NP, lowerbounds, upperbounds, max_gen)
my_algo._debug = True 

# Run the algorithm
my_algo.run()

# Visualize optimization progress
# visualize_optimization(my_algo)

# Animate population evolution
animate_population(my_algo)
# Or with contour
# animate_population_with_contour(my_algo)