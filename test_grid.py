#!/usr/bin/env python3 

import numpy as np

def compute_mahalanobis_distance(samples, O_star=100, sigma_star=10):
    """
    Compute Mahalanobis distance for each sample.

    Args:
        samples (numpy.ndarray): Array of shape (num_samples, 2).
        O_star (float): Mean value for Gaussian distribution.
        sigma_star (float): Standard deviation.

    Returns:
        numpy.ndarray: Mahalanobis distance for each sample.
        float: Average Mahalanobis distance.
    """
    mean = np.array([O_star, O_star])  # Mean vector
    inv_cov = np.linalg.inv(np.eye(2) * sigma_star**2)  # Inverse covariance matrix

    # Compute Mahalanobis distance for each sample
    diff = samples - mean  # Difference from the mean
    mahalanobis_distances = np.sqrt(np.sum(diff @ inv_cov * diff, axis=1))

    # Compute average Mahalanobis distance
    avg_mahalanobis = np.mean(mahalanobis_distances)

    return mahalanobis_distances, avg_mahalanobis

def compute_average_loglikelihood(log_likelihoods):
    """
    Compute the average log-likelihood.

    Args:
        log_likelihoods (numpy.ndarray): Array of log-likelihood values.

    Returns:
        float: Average log-likelihood.
    """
    return np.mean(log_likelihoods)


def approximate_mahalanobis_from_likelihood(log_likelihoods, C=0):
    """
    Approximate Mahalanobis distance given log-likelihood values.

    Args:
        log_likelihoods (numpy array): Log-likelihood values for each sample.
        C (float): Normalization constant, can be estimated from known distributions.

    Returns:
        numpy array: Estimated Mahalanobis distances.
        float: Average estimated Mahalanobis distance.
    """
    # Compute approximate Mahalanobis distances
    mahalanobis_distances = np.sqrt(2 * (C - log_likelihoods))

    # Compute average Mahalanobis distance
    avg_mahalanobis = np.mean(mahalanobis_distances)

    return mahalanobis_distances, avg_mahalanobis

# Example usage
samples = np.random.randn(1000, 2) * 10 + 100  # Example samples around O_star
log_likelihoods = - (samples - 100)**2 / 200  # Example log-likelihood values

mahalanobis_distances, avg_mahalanobis = compute_mahalanobis_distance(samples)
avg_loglikelihood = compute_average_loglikelihood(log_likelihoods)

print("Average Mahalanobis Distance:", avg_mahalanobis)
print("Average Log-Likelihood:", avg_loglikelihood)

mahalanobis_distances, avg_mahalanobis = approximate_mahalanobis_from_likelihood(log_likelihoods, C=0)

print("Estimated Average Mahalanobis Distance:", avg_mahalanobis)