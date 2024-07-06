import numpy as np

from ..common.fitness_evaluation import calculate_on_off_viability_scores, coarse_grain_structures
from ..common.folding_preparation import generate_segments, generate_strands, build_base_maps
from ..common.folding  import fold
from ..common.utility import multicore, flatten
from ..common.novelty import calculate_distance_matrix, calculate_novelty_scores
from .folding_preparation import generate_folding_tasks
import config

def predict_structures_evaluate_performance(individuals, folding_package_names, constant_data):
    '''
    Prepares the individuals for multicore processing (structural prediction and performance evaluation), calls the multicore function, and retrieves the results of the multicore fucntion
    '''
    
    extension_region_seg_ids, truth_vector, folding_temperature, are_segments_low_level, sensor_mismatches, folding_package_paths, ribozyme_type, template, wt_len, stride_size = constant_data

    batch = prepare_inds_batch(individuals, template, sensor_mismatches, wt_len, stride_size)
    common_data = (extension_region_seg_ids, truth_vector, folding_package_names, folding_temperature, are_segments_low_level, sensor_mismatches, folding_package_paths, ribozyme_type)
    results = multicore(batch, predict_structures_evaluate_performance_multicore, config.num_processes, common_data)

    # Store results back in individuals
    for i, ind in enumerate(individuals):
        result = results[i]
        ind.potential_fitnesses = result[0]
        ind.spmms = result[1]
        ind.mfe_dot_bracket_structures = result[2]   

def prepare_inds_batch(individuals, template, sensor_mismatches, wt_len, stride_size ):
    '''
    Packages data for each individual into a batch to be used by multiprocess fitness evaluation

    Parameters:
        inds (list of SwitchableRibozyme objects): The list of individuals whose fitness we want to evaluate
        template (RibozymeTemplate object)
        sensor_mismatches 
        wt_len (int)
        stride_size (int): The number of nucleotides by which consuective sRz + Mutant configurations shift

    Returns:
        batch (complex data structure): Contains batched data required to evaluate each individual's fitness in parallel
    '''

    batch = list()
    for ind in individuals:
        ind.segments = generate_segments(ind.sccs, template.num_strands)
        ind.strands = generate_strands(ind.segments, template.num_strands)
        ind.base_to_seg, ind.name_to_bases = build_base_maps(ind.segments)
        ind.tasks = generate_folding_tasks(sensor_mismatches, ind.strands, ind.name_to_bases, wt_len, stride_size)
        batch.append({'folding_tasks':ind.tasks, 'base_to_seg':ind.base_to_seg, 'name_to_bases':ind.name_to_bases}) # base to seg and name_to_bases will only be different between individuals during post-processing with the longer RzBS extension
    return batch

def predict_structures_evaluate_performance_multicore(input_data): 
    '''
    The core of the EA. Operates on a batch of selected individual properties.
    Predicts secondary structures for each state of every individual in the subbatch
    
    Parameters:
        Input data is split into sub_batch and common_data
            sub_batch pertains to a certain subset of the individuals in the population
            see prepare_inds_batch to see the info contained in the sub batch 
            common_data contains information shared across all individuals in the sub batch

    '''

    sub_batch, common_data = input_data
    extension_region_seg_ids, truth_vector, folding_package_names, folding_temperature, are_segments_low_level, mismatches, folding_package_paths, ribozyme_type = common_data
    num_states = len(truth_vector)

    folding_tasks = [sub_batch[i]['folding_tasks'] for i in range(len(sub_batch))] # This is a list of len(sub_batch) sublists where each sublist contains one task per state 
    folding_tasks = flatten(folding_tasks) # fold expects a single flat list of folding tasks

    # Folds the entire subbatch
    folding_results = fold(folding_tasks, folding_package_names, folding_temperature, 2, folding_package_paths) # 2 is for 2 strands so it uses cofold
   
    results = list()   
    for i in range(len(sub_batch)):    
        # Retrieve the folding results for the appropriate individual and states
        bppms = list()
        mfe_dot_bracket_structures = list()    
        for state_id in range(num_states):
            task_id = i*num_states + state_id
            (mfe_db, bppm, _, _) = folding_results[task_id] # Cofold has no ensemble diversity           
            mfe_dot_bracket_structures.append(mfe_db)           
            bppms.append(bppm)

        # Calculate the fitness scores
        base_to_seg, name_to_bases = sub_batch[i]['base_to_seg'], sub_batch[i]['name_to_bases']
        potential_fitnesses = dict()
        potential_fitnesses['ON'], potential_fitnesses['OFF'], potential_fitnesses['viability'] = calculate_on_off_viability_scores(bppms, truth_vector, name_to_bases, ribozyme_type)
        potential_fitnesses['mismatch'] = calculate_mismatch_score(bppms, mismatches, name_to_bases)

        # Coarse grain the bppms into spmms
        spmms = coarse_grain_structures(bppms, base_to_seg, extension_region_seg_ids, truth_vector, are_segments_low_level)

        results.append( ( potential_fitnesses, spmms, mfe_dot_bracket_structures) )

    return results

def calculate_mismatch_score(bppms, mismatches, name_to_bases):
    '''
    Measures the extent to which the sensor nucleotides that are mismatched with the sensor bind elsewhere
    I don' think this is a criticial score, and in future versions of TC we should explore removing it
    '''

    sensor_base_ids = name_to_bases['SENSOR'] # These are 1-indexed and right inclusive
    sensor_lims = (sensor_base_ids[0] , sensor_base_ids[-1] + 1) # These are 1-indexed and right exclusive
    # mismatches gives the indices of the sensor mismatches w.r.t. to their local position in sensor whereas mismatches_global gives their position w.r.t. the entire strand
    mismatches_global = [o + sensor_lims[0] for o in mismatches] 

    mismatch_score = 0      
    for bppm in bppms:
        for m in mismatches_global:
            mismatch_score += bppm[(m, None)]
    mismatch_score = mismatch_score / float(len(mismatches_global) )

    return mismatch_score

def calculate_distance_matrix_tc(individuals, params):
    '''
    Determines the phenotypic distance between each pair of individuals based on two metrics
    1) The manhattan distance between their corresponding segment-pair magnitude matrices (SPMMs) 
    2) Whether or not they bind to the same RzBS
        
    Parameters:
        individuals (list of SwitchableRibozyme objects)
        params['rzbs_diversity_weight'] (float) How much importance we give to two individuals having different rzbs.
            
    Returns:
        distance_matrix (num_inds X num_inds np float array) The pairwise distance between each individual
    '''

    if len(individuals) == 0:
        return None
        
    spmms_list = [ind.spmms for ind in individuals]
    rzbs_sites = [np.array([ind.rzbs_position]) for ind in individuals]
    
    # The more distinct the coarse grained secondary structures are between individuals i and j, the higher the value of entry (i,j)
    distance_matrix_spmms = calculate_distance_matrix(spmms_list, "manhattan")
    
    # Entry (i,j) is 0 if i and j have same RzBS; otherwise is 1
    distance_matrix_rzbs_sites = calculate_distance_matrix(rzbs_sites, "hamming")
    
    # Combine the two distance matrices into a single one
    # If params['rzbs_diversity_weight'] is set to 0, we ignore the rzbs matrix
    # If params['rzbs_diversity_weight'] is set to 1, entries in the sppm matrix are multiplied by 2 if they have distinct RzBSs. 
    # Entries that have the same RzBS are unmodified

    distance_matrix = distance_matrix_spmms * (1 + params['rzbs_diversity_weight']*distance_matrix_rzbs_sites)   
    return distance_matrix       
       
 
##########################################################################################################
## TO DO: Find a better place to put this function
def fold_unhybridized_srz(input_data):
    '''
    Folds the SRzs on their own. Not used in EA loop; just used for visualization purposes
    '''

    tasks, common_data = input_data  
    folding_temperature, folding_package_names, folding_package_paths = common_data
    results = fold(tasks, folding_package_names, folding_temperature, 1, folding_package_paths)
    return results
            


