import numpy as np

from .utility import arg_sort, multicore
import config

def calculate_manhattan_distance(v1, v2):
    return np.sum(np.abs(np.subtract(v1, v2)) )

def calculate_hamming_distance(v1, v2):
    return np.sum(v1 != v2)

def calculate_distance_multicore(sub_batch):
    '''
    Calculates the pairwise Manhattan or Hamming distance for a sub_batch of vectors in parallel
    '''

    tasks, common_data = sub_batch
    metric = common_data[0]
    if metric == "manhattan":
        distance_function = calculate_manhattan_distance
    elif metric == "hamming":
        distance_function = calculate_hamming_distance
    else:
        raise NotImplementedError(f"There is no distance function for distance metric: ({metric})")
        
    results = list()
    for task in tasks:
        results.append(distance_function(*task))
    return results


def calculate_distance_matrix(vectors, metric):
    '''
    Calculates the pairwise Manhattan or Hamming distance for a group of vectors
    
    Inputs:
        vectors (A list of N vectors)
        metric (string): The metric used to calculate the distance between pairs of vectors
    
    Returns
        distance_matrix (2d numpy array of size N by N). Entry (i,j) encodes the phenotypic distance between vectors i and j
    '''
    
    batch = list()
    num_vectors = len(vectors)
    # This batch preparation is inefficient (nested loop, list appending) but it is not currently a bottleneck
    for i in range(num_vectors):
        for j in range(len(vectors)):
            batch.append((vectors[i], vectors[j]))
    results = multicore(batch, calculate_distance_multicore, config.num_processes, [metric])
    #results = calculate_distance_batch((batch, [metric])) #Single core

    task_id = 0
    distance_matrix = np.zeros((num_vectors, num_vectors))    
    # Retrieve batch results. This is also inefficient, but not a bottleneck.
    for i in range(num_vectors):
        for j in range(num_vectors):
            distance_matrix[i, j] = results[task_id]
            task_id += 1  
            
    return distance_matrix

def calculate_novelty_scores(individuals, neighborhood_size, distance_matrix_func, distance_matrix_func_params):
    '''
    Calculates the novelty score for each individual in a provided population
    
    Inputs:
        individuals (list of SelectiveRibozyme Objects): The population for which we wish to calculate novelty
        neighborhood_size (int): The number of nearest neighbors that are summed together to calculate novelty
        distance_matrix_func (function): The function that calculates the phenotypic distance between each pair of individuals
        distance_matrix_func_params (dict): Parameters used by the above distance_matrix_func
    
    Returns:
        individuals: The same list of SelectiveRibozyme objects as the input, but with the novelty term calculated 
    '''
    
    print("Starting novelty")
    distance_matrix = distance_matrix_func(individuals, distance_matrix_func_params)
        
    num_items = distance_matrix.shape[0]
    novelty_scores = np.zeros(num_items)
    for i in range(num_items):
        distances = list()
        for j in range(num_items):
            distances.append(distance_matrix[i, j])                      

        neighbor_idxs = arg_sort(distances)[0:neighborhood_size] # Find the nearest neighbors
        novelty_scores[i] = sum([distances[n_idx] for n_idx in neighbor_idxs])
        
    for i, ind in enumerate(individuals):
        ind.potential_fitnesses['novelty'] = novelty_scores[i]  
        
    return individuals