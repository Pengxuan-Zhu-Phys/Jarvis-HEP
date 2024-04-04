#!/usr/bin/env python3 
import multiprocessing

def process_function(item, shared_dict):
    # Perform computation
    result = item * 2

    # Update the shared dictionary
    shared_dict["{}".format(item)] = result

if __name__ == '__main__':
    # Create a Manager object
    manager = multiprocessing.Manager()

    # Create a shared dictionary using the Manager
    shared_dict = manager.dict()

    # Create a Pool of worker processes
    pool = multiprocessing.Pool()

    # Call the process_function in parallel using map
    items = [1, 2, 3, 4, 5]
    for item in items:
        process_function(item, shared_dict)

    # Close the Pool to prevent any more tasks from being submitted
    pool.close()

    # Wait for all the worker processes to finish
    pool.join()

    # Print the shared dictionary
    print(shared_dict)
