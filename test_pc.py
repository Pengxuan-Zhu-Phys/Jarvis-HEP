#!/usr/bin/env python3 

import numpy as np
import matplotlib.pyplot as plt
import cv2
import imageio.v2 as imageio
import os
from tqdm.auto import tqdm
from scipy.spatial import cKDTree
import time 
from sympy.polys.domains import old_polynomialring
POOL_SIZE = 2000
ALPHA = 0.1  # Tunable factor for likelihood-aware resampling

# Define the eggbox likelihood function
def eggbox_likelihood(points):
    return -np.log(np.abs((np.sin(points[:, 0]) * np.sin(points[:, 1])) - 0.5) )

def dynamic_sampling(pool_size=POOL_SIZE, num_iterations=50, fps=10):
    lower_bound, upper_bound = 0 * np.pi, 1 * np.pi

    points = np.random.uniform(lower_bound, upper_bound, size=(pool_size, 2))
    previous_points = []  # Store old samples

    folder = "sampling_images"
    os.makedirs(folder, exist_ok=True)
    image_files = []
    L = eggbox_likelihood(points)

    Lmax = np.max(L)
    Lmin = np.min(L)
    for it in tqdm(range(num_iterations)):
        # Compute likelihoods
        sorted_L = np.sort(L)
        # Lmax = np.max(L)
        # Lmin = np.min(L)
        L = eggbox_likelihood(points)
        # Define new Lstar
        Lstar = Lmin + 0.1 * (Lmax - Lmin)
        # Filter samples within [1/2 Lstar, 2 Lstar]
        mask = (L >= 0.5 * Lstar) & (L <= 2 * Lstar)
        filtered_points = points[mask]

        if len(filtered_points) > 1:
            tree = cKDTree(filtered_points)
            distances, _ = tree.query(filtered_points, k=2)  # k=2 for nearest neighbor
            nearest_distances = distances[:, 1]
            dmean = np.mean(distances[:, 1])
            print("Likelihood star -> {:.6f}, Lmax -> {:.6f}, dmean -> {:.6f}".format(Lstar, Lmax, np.exp(-Lstar) ))

            nfilter = filtered_points.shape[0]
            # Resample points with dmin > dmean
            # if dmean > 0.005:
            if np.any(nearest_distances > dmean):
                print(nfilter, dmean, sep="\t", end="\r")
                resample_mask = nearest_distances > dmean
                resample_points = filtered_points[resample_mask]
                # print(nearest_distances[resample_mask])
                angles = np.random.uniform(0, 2 * np.pi, size=len(resample_points))
                displacements = (np.exp(-Lstar) / 2) * np.column_stack((np.cos(angles), np.sin(angles)))
                noise = resample_points + displacements
                # print(noise)
                
                new_samples = resample_points + noise
                new_samples = np.clip(new_samples, lower_bound, upper_bound)
                filtered_points = np.vstack([filtered_points, new_samples])
                tree = cKDTree(filtered_points)
                distances, _ = tree.query(filtered_points, k=2)
                nearest_distances = distances[:, 1]
                nfilter = filtered_points.shape[0]

            # **Combine old and new points before filtering**
            # points = filtered_points
            points = np.vstack([points, filtered_points])

            # **Retain only high-likelihood points for next iteration**
            L = eggbox_likelihood(points)
            old_point = points[L <= Lstar]
            points = points[L > Lstar]  # Keep all points with L > Lstar
            L = eggbox_likelihood(points)
            previous_points.append(old_point)

            # Update Lmax and Lmin for next iteration
            Lmax = max(np.max(L), Lmax)
            Lmin = Lstar  # Set new minimum likelihood

        # **Plot**
        plt.figure(figsize=(6, 6))
        # Plot previous points in grey
        for prev_pts in previous_points:
            plt.scatter(prev_pts[:, 0], prev_pts[:, 1], color="grey", s=(it+1)/5, alpha=0.2)
        # Plot current points with likelihood-based colors
        plt.scatter(points[:, 0], points[:, 1], c=L, cmap="viridis", marker='.', s=0.1)
        plt.colorbar(label="Likelihood")
        plt.title(f"Iteration {it + 1}")
        plt.xlim(lower_bound, upper_bound)
        plt.ylim(lower_bound, upper_bound)
        filename = os.path.join(folder, f"iter_{it:03d}.png")
        plt.savefig(filename)
        plt.close()
        
        image_files.append(filename)

    # Create video
    video_filename = "dynamic_sampling.mp4"
    frame = imageio.imread(image_files[0])
    height, width, _ = frame.shape

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(video_filename, fourcc, fps, (width, height))

    for fname in image_files:
        frame = cv2.imread(fname)
        video.write(frame)

    video.release()
    print(f"Video saved as {video_filename}")

if __name__ == "__main__":
    dynamic_sampling()
