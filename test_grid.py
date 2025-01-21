import numpy as np
import itertools

def generate_grid_samples(dimensions, num_steps):
    """
    Generate grid samples based on dimensions and dimension-specific number of steps.

    Args:
        dimensions (list of tuple): A list of tuples where each tuple specifies the (min, max) range for a dimension.
        num_steps (list of int): A list specifying the number of steps (including start and end) for each dimension.

    Returns:
        numpy.ndarray: A 2D array where each row represents a grid sample.
    """
    if len(dimensions) != len(num_steps):
        raise ValueError("The length of 'dimensions' and 'num_steps' must be the same.")

    # Generate grid points for each dimension based on the number of steps
    grid_ranges = [
        np.linspace(dim_min, dim_max, steps) 
        for (dim_min, dim_max), steps in zip(dimensions, num_steps)
    ]

    # Generate the Cartesian product of all grid points
    grid_samples = np.array(list(itertools.product(*grid_ranges)))

    return grid_samples

# Example usage
dimensions = [(0, 1), (0, 1), (0, 1)]  # Three dimensions with range [0, 1]
num_steps = [6, 4, 8]  # Number of steps for each dimension
grid_samples = generate_grid_samples(dimensions, num_steps)

print(f"Number of samples: {len(grid_samples)}")
print(grid_samples)